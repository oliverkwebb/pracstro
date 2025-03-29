//! Lunar Dynamics
use crate::{coord, sol, time};

/// Structure for the moons orbital properties at an epoch.
///
/// There's only one moon, but having the data and routines all in one type is cleaner.
/// The data contained in this is floating point instead of the Time/Coordinate types provided because static
/// Time/Angle/Coordinate data using those types is unpleasant to construct. And because the code I've written before
/// works with floats.
pub struct Moon {
    /// The epoch on which this data is based off
    pub epoch: f64,
    /// Mean Longitude
    pub l0: f64,
    /// Mean Longitude of the node
    pub n0: f64,
    /// Inclination
    pub i: time::Period,
    /// Eccentricity
    pub e: f64,
    /// Semi-major axis
    pub a: f64,
    /// Angular size at `a` distance
    pub theta0: time::Period,
    /// Parallax at `a` distance
    pub pi0: time::Period,
}

/// The moon at epoch January 1980 0.0
pub const MOON: Moon = Moon {
    epoch: 2444238.5, // January 1980 0.0
    l0: 64.975464,
    n0: 151.950429,
    i: time::Period::from_degrees(5.145396),
    e: 0.054900,
    a: 0.002569562, // AU
    theta0: time::Period::from_degrees(0.0013312900722),
    pi0: time::Period::from_degrees(0.0024428825934),
};

impl Moon {
    /// Gets a ton of information about the moon that is used by other functions
    ///
    /// From moontool.c by John Walker
    fn mooninfo(self, d: time::Date) -> (time::Period, coord::Coord, f64) {
        /* Calculation of the Sun's position */
        let day = d.julian() - self.epoch; /* Date within epoch */
        let m = time::Period::from_degrees(((360.0 / 365.2422) * day) + 278.833540 - 282.596403); /* Convert from perigee co-ordinates to epoch 1980.0 */
        let lambdasun = sol::SUN.location(d).ecliptic(d).0;

        // Moon's mean longitude
        let ml = time::Period::from_degrees(13.1763966 * day + self.l0);

        // Moon's mean anomaly
        let mm = time::Period::from_degrees(ml.degrees() - 0.1114041 * day - 349.383063); /* 349:  Mean longitude of the perigee at the epoch */

        // Evection
        let ev = 1.2739 * ((ml - lambdasun) * 2.0 - mm).sin();

        // Annual equation
        let ae = 0.1858 * m.sin();

        // Corrected anomaly
        let mmp = time::Period::from_degrees(mm.degrees() + ev - ae - (0.37 * m.sin()));

        // Correction for the equation of the centre
        let mec = time::Period::from_degrees(6.2886 * mmp.sin());

        // Corrected longitude
        let lp = time::Period::from_degrees(
            ml.degrees() + ev + (6.2886 * mmp.sin()) - ae + (0.214 * (mmp * 2.0).sin()),
        );

        // True longitude
        let lpp = time::Period::from_degrees(lp.degrees() + (1.3166 * (lp - lambdasun).sin()));

        // Corrected longitude of the node
        let np = time::Period::from_degrees(
            time::Period::from_degrees(self.n0 - 0.0529539 * day).degrees() - 0.16 * m.sin(),
        );

        let lambdamoon =
            time::Period::from_radians(((lpp - np).sin() * self.i.cos()).atan2((lpp - np).cos()))
                + np;
        let betamoon = time::Period::asin((lpp - np).sin() * self.i.sin());

        let dist = (self.a * (1.0 - self.e * self.e)) / (1.0 + self.e * (mmp + mec).cos());

        (
            // Age of the Moon in degrees
            lpp - lambdasun,
            // Coordinates of the moon
            coord::Coord::from_ecliptic(lambdamoon, betamoon, d),
            // Distance from earth
            dist,
        )
    }

    /// Returns age of phase in Days
    pub fn phaseage(self, d: time::Date) -> f64 {
        29.53058868 * self.mooninfo(d).0.turns()
    }

    /// Returns the illuminated fraction of the Moons surface
    pub fn illumfrac(self, d: time::Date) -> f64 {
        (1.0 - self.mooninfo(d).0.cos()) / 2.0
    }

    /// Returns the phase angle of the moon
    pub fn phaseangle(self, d: time::Date) -> time::Period {
        self.mooninfo(d).0
    }
    /// Returns the coordinates of the moon
    ///
    /// This code has a low accuracy of around 5 degrees
    pub fn location(self, d: time::Date) -> coord::Coord {
        self.mooninfo(d).1
    }

    /// Returns the distance to the moon in AU
    pub fn distance(self, d: time::Date) -> f64 {
        self.mooninfo(d).2
    }

    /// Returns angular diameter of the planet at current time
    pub fn angdia(self, d: time::Date) -> time::Period {
        self.theta0 / self.distance(d)
    }

    /// Calculates the moons horizontal parallax
    pub fn parallax(self, d: time::Date) -> time::Period {
        self.pi0 / self.distance(d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_moonlocation() {
        assert_eq!(
            MOON.location(time::Date::from_julian(2460748.554861)),
            coord::Coord::from_equatorial(
                time::Period::from_degminsec(172, 11, 15.7),
                time::Period::from_degminsec(3, 59, 15.2)
            )
        );
    }

    #[test]
    fn test_moonphase() {
        assert_eq!(
            MOON.illumfrac(time::Date::from_calendar(
                2025,
                03,
                29,
                time::Period::default()
            )),
            0.002799062630499616
        );
        assert_eq!(
            MOON.illumfrac(time::Date::from_calendar(
                2025,
                04,
                09,
                time::Period::default()
            )),
            0.8694887493109439
        );
    }

    #[test]
    fn test_moondist() {
        assert_eq!(
            MOON.distance(time::Date::from_julian(2460748.467894)),
            0.0026765709280575905
        );
        assert_eq!(
            MOON.angdia(time::Date::from_julian(2460748.467894)),
            time::Period::from_degrees(0.499999999)
        );
    }
}
