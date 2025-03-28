/*! Solar Dynamics and handling of planets orbits

This function has one main type, [`Planet`]. With two main methods, [`Planet::distance()`] and [`Planet::location()`].

```
use pracstro::{time, sol};

let now = time::Date::now();
for p in sol::PLANETS {
    let (ra, de) =  p.location(now).equatorial();
    let ((rah, ram, _), (ded, dem, _)) = (ra.clock(), de.to_latitude().degminsec());
    println!("{:<10} {:>2}h{:02} RA {:>3}°{:02}' De {:.2} AU", p.name, rah, ram, ded, dem, p.distance(now));
}
```

Orbital property and correction numbers from <https://ssd.jpl.nasa.gov/planets/approx_pos.html>
*/

use crate::{coord, time};

/// The location, distance, and so forth of the sun
pub mod sun {
    use crate::{coord, time};
    /// Calculate the coordinates of the sun at a given time
    ///
    /// From <https://www.celestialprogramming.com/snippets/sunPositionVsop.html>
    pub fn locationcart(d: time::Date) -> (f64, f64, f64) {
        let t = d.centuries() / 10.0;
        let mut y = 0.00010466965 * (0.09641690558 + 18_849.227_549_974_2 * t).cos();
        y += 0.00835292314 * (0.13952878991 + 12_566.151_699_982_8 * t).cos() - 0.02442699036;
        y += 0.99989211030 * (0.18265890456 + 6283.07584999140 * t).cos();
        y += (0.00093046324 + 0.00051506609 * (4.43180499286 + 12_566.151_699_982_8 * t).cos()) * t;

        let mut x =
            0.00561144206 + 0.00010466628 * (1.66722645223 + 18_849.227_549_974_2 * t).cos();
        x += 0.00835257300 * (1.71034539450 + 12_566.151_699_982_8 * t).cos();
        x += 0.99982928844 * (1.75348568475 + 6283.07584999140 * t).cos();
        x += (0.00123403056 + 0.00051500156 * (6.00266267204 + 12_566.151_699_982_8 * t).cos()) * t;

        let z = 0.00227822442 * (3.41372504278 + 6283.07584999140 * t).cos() * t;

        let tx = -(x + y * 0.000000440360 + z * -0.000000190919);
        let ty = -(x * -0.000000479966 + y * 0.917482137087 + z * -0.397776982902);
        let tz = -(y * 0.397776982902 + z * 0.917482137087);

        (tx, ty, tz)
    }

    /// Calculate the coordinates of the sun at a given time
    ///
    /// From <https://www.celestialprogramming.com/snippets/sunPositionVsop.html>
    pub fn location(d: time::Date) -> coord::Coord {
        let (x, y, z) = locationcart(d);
        coord::Coord::from_cartesian(x, y, z)
    }

    /// Calculate the distance to the sun, in AU
    pub fn distance(d: time::Date) -> f64 {
        let (tx, ty, tz) = locationcart(d);
        (tx * tx + ty * ty + tz * tz).sqrt()
    }
}

/// Generalized Planet Structure containing keplerian orbital properties and corrections.
///
/// Ephemeris for planets uses Keplerian motion with correction for perturbations of other planets
/// Error is at most 10' for most use, well within range of wanted accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct Planet {
    /// Planet Number
    pub number: u8,
    /// Planet Name
    pub name: &'static str,
    /// Semi-Major Axis (AU)
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination (Degrees)
    pub i: f64,
    /// Mean longitude (Degrees)
    pub l: f64,
    /// Longitude of the Periapsis (Degrees)
    pub w: f64,
    /// Longitude of the ascending node (Degrees)
    pub o: f64,
    /// Correction rates for all 6 preceding properties.
    pub rates: [f64; 6],
    /// Correction values for the mean anomaly, needed in larger planets
    pub extra: Option<(f64, f64, f64, f64)>,
    // Physical Properties
    /// Angular Diameter at 1AU (Degrees)
    pub theta0: f64,
    /// Visual Magnitude at 1AU
    pub v0: f64,
}
impl Planet {
    /// Returns the location of the planets as rectangular coordinates as relative to the Sun
    ///
    /// From <>
    pub fn locationcart(&self, d: time::Date) -> (f64, f64, f64) {
        fn kepler(m: f64, e: f64, ee: f64) -> f64 {
            let dm = m - (ee - e.to_degrees() * (ee.to_radians().sin()));
            dm / (1.0 - e * (ee.to_radians()).cos())
        }
        let t = d.centuries();
        let a = self.a + self.rates[0] * t;
        let e = self.e + self.rates[1] * t;
        let l = self.l + self.rates[3] * t;
        let w = self.w + self.rates[4] * t;
        let o = self.o + self.rates[5] * t;
        let i = (self.i + self.rates[2] * t).to_radians();

        let ww = (w - o).to_radians();
        let mut m = l - w;
        if let Some((b, c, s, f)) = self.extra {
            m = m + b * t * t + c * ((f * t).to_radians().cos()) + s * ((f * t).to_radians().sin());
        }
        while m > 180.0 {
            m -= 360.0;
        }

        let mut ee = m + 57.29578 * e * (m.to_radians().sin());
        let mut de: f64 = 1.0;
        while de.abs() > 1e-7 {
            de = kepler(m, e, ee);
            ee += de;
        }

        let xp = a * ((ee.to_radians()).cos() - e);
        let yp = a * (1.0 - e * e).sqrt() * (ee.to_radians().sin());

        let o = o.to_radians();
        let xecl = (ww.cos() * o.cos() - ww.sin() * o.sin() * i.cos()) * xp
            + (-ww.sin() * o.cos() - ww.cos() * o.sin() * i.cos()) * yp;
        let yecl = (ww.cos() * o.sin() + ww.sin() * o.cos() * i.cos()) * xp
            + (-ww.sin() * o.sin() + ww.cos() * o.cos() * i.cos()) * yp;
        let zecl = (ww.sin() * i.sin()) * xp + (ww.cos() * i.sin()) * yp;

        let eps = 23.43928_f64.to_radians();
        let tx = xecl;
        let ty = eps.cos() * yecl - eps.sin() * zecl;
        let tz = eps.sin() * yecl + eps.cos() * zecl;

        (tx, ty, tz)
    }

    /// Returns coordinates as subtracted from the earths coordinates
    pub fn location(&self, d: time::Date) -> coord::Coord {
        let c = self.locationcart(d);
        if self.number == 2 {
            // If we aren't the earth, get the coords of the earth
            return coord::Coord::default();
        }
        let e = EARTH.locationcart(d);

        coord::Coord::from_cartesian(c.0 - e.0, c.1 - e.1, c.2 - e.2)
    }

    /// Returns distance in AU
    pub fn distance(&self, d: time::Date) -> f64 {
        if self.number == 2 {
            // If we aren't the earth, get the coords of the earth
            return 0.0;
        }
        let c = self.locationcart(d);
        let e = EARTH.locationcart(d);
        let (tx, ty, tz) = (c.0 - e.0, c.1 - e.1, c.2 - e.2);

        (tx * tx + ty * ty + tz * tz).sqrt()
    }

    /// Returns angular diameter of the planet at current time
    pub fn angdia(&self, d: time::Date) -> time::Period {
        time::Period::from_degrees(self.theta0 / self.distance(d))
    }

    fn sun_distance(&self, d: time::Date) -> f64 {
        let (tx, ty, tz) = self.locationcart(d);
        (tx * tx + ty * ty + tz * tz).sqrt()
    }

    /// Get apparent magnitude of a planet
    pub fn magnitude(&self, d: time::Date) -> f64 {
        5.0 * ((self.distance(d) * self.sun_distance(d)) / self.illumfrac(d).sqrt()).log10()
            + self.v0
    }

    /// Gets the phase angle of a planet
    ///
    /// This is simple trig work with the triangle between the planet, earth, and sun.
    pub fn phaseangle(&self, d: time::Date) -> time::Period {
        let sep = sun::location(d).dist(self.location(d));
        let (tx, ty, tz) = self.locationcart(d);
        let sp = (tx * tx + ty * ty + tz * tz).sqrt();
        // Todo: Replace one with distance to the sun in AU
        time::Period::asin(1.0 * (sep.sin() / sp))
    }

    /// Gets the illuminated fraction of the planets surface
    pub fn illumfrac(&self, d: time::Date) -> f64 {
        // Todo: Replace one with distance to the sun in AU
        if self.number > 2 {
            0.5 * (1.0 + self.phaseangle(d).cos())
        } else {
            0.5 * (1.0 - self.phaseangle(d).cos())
        }
    }
}

/// Mercury
pub const MERCURY: Planet = Planet {
    number: 0,
    name: "Mercury",
    a: 0.38709927,
    e: 0.20563593,
    i: 7.00497902,
    l: 252.25032350,
    w: 77.45779628,
    o: 48.33076593,
    rates: [
        0.00000037,
        0.00001906,
        -0.00594749,
        149472.67411175,
        0.16047689,
        -0.12534081,
    ],
    extra: None,
    theta0: 0.0017972222,
    v0: -0.42,
};
/// Venus
pub const VENUS: Planet = Planet {
    number: 1,
    name: "Venus",
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
    extra: None,
    theta0: 0.0047,
    v0: -4.4,
};
/// Earth (Technically the Earth-Moon Barycenter)
pub const EARTH: Planet = Planet {
    number: 2,
    name: "Earth",
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
    extra: None,
    theta0: 180.0,
    v0: -12.0,
};
/// Mars
pub const MARS: Planet = Planet {
    number: 3,
    name: "Mars",
    a: 1.52371034,
    e: 0.09339410,
    i: 1.84969142,
    l: -4.55343205,
    w: -23.94362959,
    o: 49.55953891,
    rates: [
        0.00001847,
        0.00007882,
        -0.00813131,
        19140.30268499,
        0.44441088,
        -0.29257343,
    ],
    extra: None,
    theta0: 0.0026,
    v0: -1.52,
};
/// Jupiter
pub const JUPITER: Planet = Planet {
    number: 4,
    name: "Jupiter",
    a: 5.20248019,
    e: 0.04853590,
    i: 1.29861416,
    l: 34.33479152,
    w: 14.27495244,
    o: 100.29282654,
    rates: [
        -0.00002864,
        0.00018026,
        -0.00322699,
        3034.90371757,
        0.18199196,
        0.13024619,
    ],
    extra: Some((-0.00012452, 0.06064060, -0.35635438, 38.35125000)),
    theta0: 0.05465,
    v0: -9.4,
};
/// Saturn
pub const SATURN: Planet = Planet {
    number: 5,
    name: "Saturn",
    a: 9.54149883,
    e: 0.05550825,
    i: 2.49424102,
    l: 50.07571329,
    w: 92.86136063,
    o: 113.63998702,
    rates: [
        -0.00003065,
        -0.00032044,
        0.00451969,
        1222.11494724,
        0.54179478,
        -0.25015002,
    ],
    extra: Some((0.00025899, -0.13434469, 0.87320147, 38.35125000)),
    theta0: 0.046,
    v0: -8.9,
};
/// Uranus
pub const URANUS: Planet = Planet {
    number: 6,
    name: "Uranus",
    a: 19.18797948,
    e: 0.04685740,
    i: 0.77298127,
    l: 314.20276625,
    w: 172.43404441,
    o: 73.96250215,
    rates: [
        -0.00020455,
        -0.00001550,
        -0.00180155,
        428.49512595,
        0.09266985,
        0.05739699,
    ],
    extra: Some((0.00058331, -0.97731848, 0.17689245, 7.67025000)),
    theta0: 0.0182777777,
    v0: -7.19,
};
/// Neptune
pub const NEPTUNE: Planet = Planet {
    number: 7,
    name: "Neptune",
    a: 30.06952752,
    e: 0.00895439,
    i: 1.77005520,
    l: 304.22289287,
    w: 46.68158724,
    o: 131.78635853,
    rates: [
        0.00006447,
        0.00000818,
        0.00022400,
        218.46515314,
        0.01009938,
        -0.00606302,
    ],
    extra: Some((-0.00041348, 0.68346318, -0.10162547, 7.67025000)),
    theta0: 0.0172777777,
    v0: -6.87,
};
/// Pluto
pub const PLUTO: Planet = Planet {
    number: 8,
    name: "Pluto",
    a: 39.48686035,
    e: 0.24885238,
    i: 17.14104260,
    l: 238.96535011,
    w: 224.09702598,
    o: 110.30167986,
    rates: [
        0.00449751,
        0.00006016,
        0.00000501,
        145.18042903,
        -0.00968827,
        -0.00809981,
    ],
    extra: None,
    theta0: 0.0022777777,
    v0: -1.0,
};

/// Defines the planets in order
///
/// Can be used in a iterator to loop over planets
///
/// ```
/// use pracstro::*;
/// for p in sol::PLANETS {
///   println!("{} AU", p.distance(time::Date::now()));
/// }
/// ```
pub const PLANETS: [&Planet; 9] = [
    &MERCURY, &VENUS, &EARTH, &MARS, &JUPITER, &SATURN, &URANUS, &NEPTUNE, &PLUTO,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sunpos() {
        assert_eq!(
            sun::location(time::Date::from_julian(2268932.541667)),
            coord::Coord::from_equatorial(
                time::Period::from_degrees(298.49306),
                time::Period::from_degrees(-20.91664)
            )
        );
    }

    // "Is this a reliable way of getting the ecliptic longitude of the sun?"
    #[test]
    fn test_lambdasun() {
        assert_eq!(
            sun::location(time::Date::from_calendar(1980, 7, 27, time::Period::default()))
                .ecliptic(time::Date::from_calendar(1980, 7, 27, time::Period::default()))
                .0,
            time::Period::from_degminsec(124, 23, 40.8)
        )
    }

    #[test]
    fn test_planet() {
        assert_eq!(
            VENUS.location(time::Date::from_calendar(2025, 3, 12, time::Period::default())),
            coord::Coord::from_equatorial(
                time::Period::from_clock(0, 17, 44.5),
                time::Period::from_degminsec(10, 54, 50.7)
            )
        );
        assert_eq!(
            JUPITER.location(time::Date::from_julian(2460748.41871)),
            coord::Coord::from_equatorial(
                time::Period::from_clock(4, 47, 10.5),
                time::Period::from_degminsec(22, 01, 7.7)
            )
        );
        assert_eq!(
            MARS.distance(time::Date::from_julian(2460748.41871)),
            0.9721731869765928
        );
        assert_eq!(
            JUPITER.distance(time::Date::from_julian(2460748.41871)),
            5.1839327273287825
        );
    }

    #[test]
    fn test_phase() {
        assert_eq!(
            VENUS.illumfrac(time::Date::from_calendar(2025, 03, 24, time::Period::default())),
            0.010588440315090597
        );
        assert_eq!(
            MARS.illumfrac(time::Date::from_calendar(2025, 03, 24, time::Period::default())),
            0.9097044984005929
        );
        assert_eq!(
            VENUS.illumfrac(time::Date::from_calendar(1996, 07, 22, time::Period::default())),
            0.2930458224882831
        );
        assert_eq!(
            VENUS.illumfrac(time::Date::from_calendar(1996, 06, 22, time::Period::default())),
            0.04343209659342473
        );
        assert_eq!(
            VENUS.illumfrac(time::Date::from_calendar(1996, 06, 22, time::Period::default())),
            0.04343209659342473
        );
    }
}
