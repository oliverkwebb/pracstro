use crate::{coord, time};

pub struct Cart(f64, f64, f64);
impl Cart {
    pub fn celestial(self) -> coord::Coord {
        let Cart(tx, ty, tz) = self;
        let r = (tx * tx + ty * ty + tz * tz).sqrt();
        let l = time::Period::from_radians(ty.atan2(tx));
        let t2 = time::Period::from_radians(0.5 * std::f64::consts::PI - (tz / r).acos());

        coord::Coord::from_celestial(l, t2)
    }
}

/// Calculate the coordinates of the sun at a given time
///
/// From https://www.celestialprogramming.com/snippets/sunPositionVsop.html
pub fn whereis_sun(d: time::Date) -> coord::Coord {
    let t = (d.julian() - 2451545.0) / 365250.0;
    let mut y = 0.00010466965 * (0.09641690558 + 18849.22754997420 * t).cos();
    y = y + 0.00835292314 * (0.13952878991 + 12566.15169998280 * t).cos() - 0.02442699036;
    y = y + 0.99989211030 * (0.18265890456 + 6283.07584999140 * t).cos();
    y = y + (0.00093046324 + 0.00051506609 * (4.43180499286 + 12566.15169998280 * t).cos()) * t;

    let mut x = 0.00561144206 + 0.00010466628 * (1.66722645223 + 18849.22754997420 * t).cos();
    x = x + 0.00835257300 * (1.71034539450 + 12566.15169998280 * t).cos();
    x = x + 0.99982928844 * (1.75348568475 + 6283.07584999140 * t).cos();
    x = x + (0.00123403056 + 0.00051500156 * (6.00266267204 + 12566.15169998280 * t).cos()) * t;

    let z = 0.00227822442 * (3.41372504278 + 6283.07584999140 * t).cos() * t;

    let tx = -(x + y * 0.000000440360 + z * -0.000000190919);
    let ty = -(x * -0.000000479966 + y * 0.917482137087 + z * -0.397776982902);
    let tz = -(y * 0.397776982902 + z * 0.917482137087);

    let r = (tx * tx + ty * ty + tz * tz).sqrt();
    let l = time::Period::from_radians(ty.atan2(tx));
    let t2 = time::Period::from_radians(0.5 * std::f64::consts::PI - (tz / r).acos());

    coord::Coord::from_celestial(l, t2)
}

// Ephemeris for planets uses Keplerian motion with correction for perturbations of other planets
// Error is at most 10' for most use, well within range of wanted accuracy.
#[derive(PartialEq, Copy, Clone)]
pub struct Planet {
    number: u8, // Planet Number
    a: f64,
    e: f64,
    i: f64,
    l: f64,
    w: f64,
    o: f64,
    rates: [f64; 6],
}
impl Planet {
    pub fn locationcart(self, d: time::Date) -> Cart {
        fn kepler(m: f64, e: f64, ee: f64) -> f64 {
            let dM = m - (ee - e.to_degrees() * (ee.to_radians().sin()));
            dM / (1.0 - e * (ee.to_radians()).cos())
        }
        let T = (d.julian() - 2451545.0) / 36525.0;
        let a = self.a + self.rates[0] * T;
        let e = self.e + self.rates[1] * T;
        let l = self.l + self.rates[3] * T;
        let w = self.w + self.rates[4] * T;
        let o = self.o + self.rates[5] * T;
        let i = (self.i + self.rates[2] * T).to_radians();

        let ww = (w - o).to_radians();
        let mut m = l - w;

        while m > 180.0 {
            m = m - 360.0;
        }

        let mut E = m + 57.29578 * e * (m.to_radians().sin());
        let mut de: f64 = 1.0;
        while de.abs() > 1e-7 {
            de = kepler(m, e, E);
            E = E + de;
        }

        let xp = a * ((E.to_radians()).cos() - e);
        let yp = a * (1.0 - e * e).sqrt() * (E.to_radians().sin());

        let O = o.to_radians();
        let xecl = (ww.cos() * O.cos() - ww.sin() * O.sin() * i.cos()) * xp
            + (-ww.sin() * O.cos() - ww.cos() * O.sin() * i.cos()) * yp;
        let yecl = (ww.cos() * O.sin() + ww.sin() * O.cos() * i.cos()) * xp
            + (-ww.sin() * O.sin() + ww.cos() * O.cos() * i.cos()) * yp;
        let zecl = (ww.sin() * i.sin()) * xp + (ww.cos() * i.sin()) * yp;

        let eps = 23.43928_f64.to_radians();
        let tx = xecl;
        let ty = eps.cos() * yecl - eps.sin() * zecl;
        let tz = eps.sin() * yecl + eps.cos() * zecl;

        Cart(tx, ty, tz)
    }

    pub fn location(self, d: time::Date) -> coord::Coord {
        let c = self.locationcart(d);
        if self.number == 2 {
            // If we aren't the earth, get the coords of the earth
            return c.celestial();
        }
        let e = earth.locationcart(d);

        Cart(c.0 - e.0, c.1 - e.1, c.2 - e.2).celestial()
    }
}

const venus: Planet = Planet {
    number: 1,
    a: 0.72333566,
    e: 0.00677672,
    i: 3.39467605,
    l: 181.97909950,
    w: 131.60246718,
    o: 76.67984255,
    rates: [
        0.00000390,
        -0.00004107,
        -0.00078890,
        58517.81538729,
        0.00268329,
        -0.27769418,
    ],
};
const earth: Planet = Planet {
    number: 2,
    a: 1.00000261,
    e: 0.01671123,
    i: -0.00001531,
    l: 100.46457166,
    w: 102.93768193,
    o: 0.0,
    rates: [
        0.00000562,
        -0.00004392,
        -0.01294668,
        35999.37244981,
        0.32327364,
        0.0,
    ],
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sunpos() {
        assert_eq!(
            whereis_sun(time::Date::from_julian(2268932.541667)),
            coord::Coord::from_celestial(
                time::Period::from_degrees(298.49306),
                time::Period::from_degrees(-20.91664)
            )
        );
    }

    // "Is this a reliable way of getting the ecliptic longitude of the sun?"
    #[test]
    fn test_lambdasun() {
        assert_eq!(
            whereis_sun(time::Date::from_calendar(1980, 7, 27.0))
                .ecliptic(time::Date::from_calendar(1980, 7, 27.0))
                .0,
            time::Period::from_degminsec(124, 23, 40.8)
        )
    }

    #[test]
    fn test_planet() {
        assert_eq!(
            venus.location(time::Date::from_calendar(2025, 3, 12.0)),
            coord::Coord::from_celestial(
                time::Period::from_clock(0, 16, 11.5),
                time::Period::from_degminsec(10, 52, 50.7)
            )
        );
    }
}
