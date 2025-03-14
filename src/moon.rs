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
    pub i: f64,
    /// Eccentricity
    pub e: f64,
    /// Semi-major axis
    pub a: f64,
    /// Angular size at `a` distance
    pub theta0: f64,
    /// Parallax at `a` distance
    pub pi0: f64,
}

/// The moon at epoch January 1980 0.0
pub const MOON: Moon = Moon {
    epoch: 2444238.5, // January 1980 0.0
    l0: 64.975464,
    n0: 151.950429,
    i: 5.145396,
    e: 0.054900,
    a: 384401.0, // KM
    theta0: 0.5181,
    pi0: 0.9507,
};

impl Moon {
    /// Gets a ton of information about the moon that is used by other functions
    ///
    /// From moontool.c by John Walker
    pub fn mooninfo(self, d: time::Date) -> (f64, coord::Coord, f64) {
        fn fixangle(a: f64) -> f64 {
            a - 360.0 * (a / 360.0).floor()
        }
        /* Calculation of the Sun's position */
        let day = d.julian() - self.epoch; /* Date within epoch */
        let m = fixangle(fixangle((360.0 / 365.2422) * day) + 278.833540 - 282.596403); /* Convert from perigee co-ordinates to epoch 1980.0 */
        let lambdasun = sol::where_is_sun(d).ecliptic(d).0.degrees();

        // Moon's mean longitude
        let ml = time::Period::from_degrees(13.1763966 * day + self.l0);

        // Moon's mean anomaly
        let mm = time::Period::from_degrees(ml.degrees() - 0.1114041 * day - 349.383063); /* 349:  Mean longitude of the perigee at the epoch */

        // Evection
        let ev = 1.2739
            * (2.0 * (ml.degrees() - lambdasun) - mm.degrees())
                .to_radians()
                .sin();

        // Annual equation
        let ae = 0.1858 * m.to_radians().sin();

        // Corrected anomaly
        let mmp = mm.degrees() + ev - ae - (0.37 * m.to_radians().sin());

        // Correction for the equation of the centre
        let mec = 6.2886 * mmp.to_radians().sin();

        // Corrected longitude
        let lp = ml.degrees() + ev + (6.2886 * mmp.to_radians().sin()) - ae
            + (0.214 * (2.0 * mmp).to_radians().sin());

        // True longitude
        let lpp = lp + (0.6583 * (2.0 * (lp - lambdasun).to_radians().sin()));

        // Corrected longitude of the node
        let np = fixangle(self.n0 - 0.0529539 * day) - 0.16 * m.to_radians().sin();

        let lambdamoon = (((lpp - np).to_radians().sin() * self.i.to_radians().cos())
            .atan2((lpp - np).to_radians().cos()))
        .to_degrees()
            + np;
        let betamoon =
            (lpp - np).to_radians().sin().asin().to_degrees() * self.i.to_radians().sin();

        let dist =
            (self.a * (1.0 - self.e * self.e)) / (1.0 + self.e * (mmp + mec).to_radians().cos());

        (
            // Age of the Moon in degrees
            lpp - lambdasun,
            // Coordinates of the moon
            coord::Coord::from_ecliptic(
                time::Period::from_degrees(lambdamoon),
                time::Period::from_degrees(betamoon),
                d,
            ),
            // Distance from earth
            dist,
        )
    }

    /// Returns (Illuminated Fraction, Moon Age)
    pub fn phase(self, d: time::Date) -> (f64, f64) {
        let (age, _, _) = self.mooninfo(d);
        fn fixangle(a: f64) -> f64 {
            a - 360.0 * (a / 360.0).floor()
        }
        /// Synodic month (new Moon to new Moon)
        const SYNMONTH: f64 = 29.53058868;
        (
            (1.0 - age.to_radians().cos()) / 2.0,
            SYNMONTH * (fixangle(age) / 360.0),
        )
    }

    /// Returns the coordinates of the moon
    ///
    /// Note that this code has a low accuracy of around 5 degrees
    pub fn location(self, d: time::Date) -> coord::Coord {
        self.mooninfo(d).1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_moonlocation() {
        /*
        assert_eq!(
            MOON.location(time::Date::from_julian(2446982.804861111)),
            coord::Coord::from_celestial(
                time::Period::from_degminsec(227, 1, 31.55),
                time::Period::from_degminsec(-17, 35, 39.2)
            )
        );
        */
        /*
        assert_eq!(
            MOON.location(time::Date::from_julian(2453922.554861111)),
            coord::Coord::from_celestial(
                time::Period::from_clock(14, 26, 57.7),
                time::Period::from_degminsec(-18, 05, 39.2)
            )
        );
        */
        assert_eq!(
            MOON.location(time::Date::from_julian(2460748.554861)),
            coord::Coord::from_celestial(
                time::Period::from_degminsec(172, 11, 15.7),
                time::Period::from_degminsec(3, 59, 15.2)
            )
        );
    }

    #[test]
    fn test_moonphase() {
        assert_eq!(
            MOON.phase(time::Date::from_julian(2460748.467894)),
            (0.9992826878175012, 14.513510258299618)
        );
    }
}
