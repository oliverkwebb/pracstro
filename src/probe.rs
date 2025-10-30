/*! Segmented Solar Dynamics and handling of planets orbits

This function has one main type, [`SegmentedPlanet`] With methods for:

* Cartesian Coordinates
* Distance from earth

Orbital property and correction numbers from <https://ssd.jpl.nasa.gov/planets/approx_pos.html>
and JPL Horizons <https://ssd.jpl.nasa.gov/horizons/>
*/

use crate::{coord, sol::EARTH, time};

/// Generalized Planet Structure containing keplerian orbital properties and corrections.
///
/// Ephemeris for planets uses Keplerian motion with correction for perturbations of other planets
/// Error is at most 10' for most use, well within range of wanted accuracy.
#[derive(Clone, Debug, PartialEq)]
pub struct SegmentedPlanet {
    /// Planet Name
    pub name: &'static str,
    /// Semi-Major Axis (AU)
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination (Degrees)
    pub i: f64,
    /// Longitude of the Periapsis (Degrees)
    pub w: f64,
    /// Longitude of the ascending node (Degrees)
    pub o: f64,
    /// Mean longitude (Degrees)
    pub l: f64,
    /// Correction rates of the mean longitude per Julian Century (degrees)
    pub l_delta_century: f64,
    /// Epoch of the Mean Longitude
    pub l_epoch: time::Date,
}
impl SegmentedPlanet {
    /// Returns the location of the planets as rectangular coordinates as relative to the Sun, in AU
    ///
    /// From <https://ssd.jpl.nasa.gov/planets/approx_pos.html>
    pub fn locationcart(&self, d: time::Date) -> (f64, f64, f64) {
        let t = (d.julian() - self.l_epoch.julian()) / 36525.0;
        let a = self.a;
        let e = self.e;
        let i = time::Angle::from_degrees(self.i);
        let o = time::Angle::from_degrees(self.o);
        let w = time::Angle::from_degrees(self.w);
        let l = time::Angle::from_degrees(self.l + (self.l_delta_century * t));
        let ww = w - o;
        let mut m = (l - w).degrees();
        m = time::Angle::from_degrees(m).to_latitude().degrees();

        fn kepler(m: f64, e: f64, ee: f64) -> f64 {
            let dm = m - (ee - e.to_degrees() * (ee.to_radians().sin()));
            dm / (1.0 - e * (ee.to_radians()).cos())
        }
        let mut ee = m + 57.29578 * e * (m.to_radians().sin());
        let mut de: f64 = 1.0;
        while de.abs() > 1e-7 {
            de = kepler(m, e, ee);
            ee += de;
        }

        if e < 1.0 {
            let v = 2.0
                * ((1.0 + e).sqrt() * (ee.to_radians() / 2.0).sin())
                    .atan2((1.0 - e).sqrt() * (ee.to_radians() / 2.0).cos());
            eprintln!("{}", l.degrees());
            let xp = a * ((ee.to_radians()).cos() - e);
            let yp = a * (1.0 - e * e).sqrt() * (ee.to_radians().sin());

            let xecl = (ww.cos() * o.cos() - ww.sin() * o.sin() * i.cos()) * xp
                + (-ww.sin() * o.cos() - ww.cos() * o.sin() * i.cos()) * yp;
            let yecl = (ww.cos() * o.sin() + ww.sin() * o.cos() * i.cos()) * xp
                + (-ww.sin() * o.sin() + ww.cos() * o.cos() * i.cos()) * yp;
            let zecl = (ww.sin() * i.sin()) * xp + (ww.cos() * i.sin()) * yp;

            let eps = 23.43928_f64.to_radians();
            let tx = xecl;
            let ty = eps.cos() * yecl - eps.sin() * zecl;
            let tz = eps.sin() * yecl + eps.cos() * zecl;

            return (tx, ty, tz);
        } else {
            todo!();
        }
    }

    /// Returns coordinates as subtracted from the earths coordinates
    pub fn location(&self, d: time::Date) -> coord::Coord {
        let c = self.locationcart(d);
        let e = EARTH.locationcart(d);

        coord::Coord::from_cartesian(c.0 - e.0, c.1 - e.1, c.2 - e.2)
    }

    /// Returns distance in AU
    pub fn distance(&self, d: time::Date) -> f64 {
        let c = self.locationcart(d);
        let e = EARTH.locationcart(d);
        let (tx, ty, tz) = (c.0 - e.0, c.1 - e.1, c.2 - e.2);

        (tx * tx + ty * ty + tz * tz).sqrt()
    }

    fn sun_distance(&self, d: time::Date) -> f64 {
        let (tx, ty, tz) = self.locationcart(d);
        (tx * tx + ty * ty + tz * tz).sqrt()
    }
}

/// Voyager 2 Test Object
pub const VOYAGER2TEST: SegmentedPlanet = SegmentedPlanet {
    name: "Voyager 2 Test of 1983-11-30",
    a: -3.922739981,
    e: 3.44861399869959,
    i: 2.66291222081523,
    w: 77.45779628,
    o: 48.33076593,
    l: 298.15,
    l_delta_century: 55390.810728252,
    l_epoch: time::Date::from_julian(2445668.5),
};
/// Mars
pub const MARS: SegmentedPlanet = SegmentedPlanet {
    name: "Mars Control Test Object",
    a: 1.52371034,
    e: 0.09339410,
    i: 1.84969142,
    l: -4.55343205,
    w: -23.94362959,
    o: 49.55953891,
    l_delta_century: 19140.30268499,
    l_epoch: time::Date::from_julian(2451545.0),
};
/// Parker Solar Probe Test Object
pub const PARKERTEST: SegmentedPlanet = SegmentedPlanet {
    name: "Parker Solar Probe Test of 2025-09-09",
    a: 0.3884911788,
    e: 0.881936788230005,
    i: 3.39525331656616,
    w: 68.4610434250691,
    o: 76.6442910488463,
    l: 477.32,
    l_delta_century: 148669.7918,
    l_epoch: time::Date::from_julian(2460928.5),
};
/// Halleys Commet
pub const HALLEY: SegmentedPlanet = SegmentedPlanet {
    name: "Halley's Commet",
    a: 17.8591256074516,
    e: 0.967894644637869,
    i: 162.151115110664,
    o: 59.5723978313194,
    w: 112.470381868335,
    l: 360.1,
    l_delta_century: 365.25,
    l_epoch: time::Date::from_julian(2460928.5),
};
/// Mars, Again
pub const SUPERSURE: SegmentedPlanet = SegmentedPlanet {
    name: "Mars Again",
    a: 1.52362660139838,
    e: 0.0934994380486591,
    i: 1.84752995631806,
    o: 49.483749371897,
    w: 286.651495403059,
    l: 592.71,
    l_delta_century: 19141.02,
    l_epoch: time::Date::from_julian(2460927.5),
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_voy2831130() {
        return;
        let x = SUPERSURE.location(time::Date::from_calendar(
            2025,
            9,
            9,
            time::Angle::from_degrees(0.0),
        ));
        eprintln!(
            "{}, {}",
            x.equatorial().0.decimal(),
            x.equatorial().1.to_latitude().degrees()
        );
        todo!();
    }
}
