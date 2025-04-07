//! Coordinate handling
//!
//! This module contains one type, [`Coord`]. That has methods to convert two and from several
//! different coordinate systems. Mainly:
//! - Equatorial (Hour Angle, Declination)
//! - Horizon (Azimuth, Altitude)
//! - Ecliptic (Beta, Lambda)
//!
//! This type also contains algorithms for converting from Cartesian (rectangular) coordinates, rise and set times, distance between angles, etc.

use crate::time::*;

/// Gets the mean obliquity of the ecliptic at a certain date
pub fn mean_obliquity_ecl(d: Date) -> Angle {
    let t = d.centuries();
    Angle::from_degrees(
        23.439_292 - ((46.815 * t + 0.0006 * (t * t) - 0.00181 * (t * t * t)) / 3600.0),
    )
}

/**
Pair of angles, Representing "How far up" and "How far round"

| Property          | Latitude          | Longitude           | Depends On                      | To Method              | From Method                 |
|-------------------|-------------------|---------------------|---------------------------------|------------------------|-----------------------------|
| Equatorial        | Declination (δ)   | Right Ascension (α) |                                 | [`Coord::equatorial()`]| [`Coord::from_equatorial()`]|
| Horizontal        | Altitude (a)      | Azimuth (A)         | Date, Time, Latitude, Longitude | [`Coord::horizon()`]   | [`Coord::from_horizon()`]   |
| Ecliptic          | Ecl. Latitude (β) | Ecl. Longitude (λ)  | Date[^1]                        | [`Coord::ecliptic()`]  | [`Coord::from_ecliptic()`]  |
| Cartesian         | N/A (3D system)   | N/A (3D system)     | Distance                        |                        | [`Coord::from_cartesian()`] |

Additional Methods:
* Distance between coordinates: [`Coord::dist()`]
* Rise and set times of a coordinate in the sky [`Coord::riseset()`]
* Precession [`Coord::precess()`]

[^1]: The plane of the ecliptic varies slightly with perturbations in the orbit and inclination of the earth.
*/
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct Coord(Angle, Angle);
impl Coord {
    /// Right Ascension and Declination
    pub const fn equatorial(self) -> (Angle, Angle) {
        (self.0, self.1)
    }
    /// Right Ascension and Declination
    pub const fn from_equatorial(x: Angle, y: Angle) -> Self {
        Coord(x, y)
    }

    /// Azimuth and Altitude, dependent on longitude, Latitude and time
    ///
    /// From Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn horizon(self, date: Date, time: Angle, lati: Angle, longi: Angle) -> (Angle, Angle) {
        let (ra, de) = self.equatorial();
        let ha = time.gst(date) + longi - ra;
        let alt = Angle::asin(de.sin() * lati.sin() + de.cos() * lati.cos() * ha.cos());
        let azip = Angle::acos((de.sin() - lati.sin() * alt.sin()) / (lati.cos() * alt.cos()));
        let azi = match ha.sin() < 0.0 {
            true => azip,
            false => Angle::from_degrees(360.0 - azip.degrees()),
        };
        (azi, alt)
    }
    /// Azimuth and Altitude, dependent on longitude, Latitude and time
    ///
    /// From Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn from_horizon(
        azi: Angle,
        alt: Angle,
        date: Date,
        time: Angle,
        lati: Angle,
        longi: Angle,
    ) -> Self {
        let de = Angle::asin(alt.sin() * lati.sin() + alt.cos() * lati.cos() * azi.cos());
        let hap = Angle::acos((alt.sin() - lati.sin() * de.sin()) / (lati.cos() * de.cos()));
        let ha = match azi.sin() < 0.0 {
            true => hap,
            false => Angle::from_degrees(360.0 - hap.degrees()),
        };
        Coord::from_equatorial(time.gst(date) + longi - ha, de)
    }

    /// Used in solar calculations, based on the plane of the orbit of the earth
    ///
    /// From Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn ecliptic(self, d: Date) -> (Angle, Angle) {
        let (ra, de) = self.equatorial();
        let e = mean_obliquity_ecl(d);
        let beta = Angle::asin(de.sin() * e.cos() - de.cos() * e.sin() * ra.sin());
        let y = ra.sin() * e.cos() + de.tan() * e.sin();
        let x = ra.cos();
        let lambda = Angle::atan2(y, x);
        (lambda, beta)
    }
    /// Used in solar calculations, based on the plane of the orbit of the earth
    ///
    /// From Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn from_ecliptic(lambda: Angle, beta: Angle, d: Date) -> Self {
        let e = mean_obliquity_ecl(d);
        let de = Angle::asin(beta.sin() * e.cos() + beta.cos() * e.sin() * lambda.sin());
        let ra = Angle::atan2(lambda.sin() * e.cos() - beta.tan() * e.sin(), lambda.cos());
        Coord::from_equatorial(ra, de)
    }

    /// Convert Rectangular Coordinates to RA/Dec
    ///
    /// Note how this has no pair function that converts to rectangular coords
    pub fn from_cartesian(x: f64, y: f64, z: f64) -> Self {
        let (tx, ty, tz) = (x, y, z);
        let r = (tx * tx + ty * ty + tz * tz).sqrt();
        let l = Angle::atan2(ty, tx);
        let t2 = Angle::from_radians(0.5 * std::f64::consts::PI - (tz / r).acos());

        Coord::from_equatorial(l, t2)
    }

    /// Returns the angle between two objects
    pub fn dist(self, from: Self) -> Angle {
        let ((a1, d1), (a2, d2)) = (self.equatorial(), from.equatorial());
        Angle::acos(d1.sin() * d2.sin() + d1.cos() * d2.cos() * (a1 - a2).cos())
    }
    /// Returns (Rise, Set) UT, This function will fail for locations in the sky that never appear over the horizon
    ///
    /// From Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn riseset(self, date: Date, lati: Angle, longi: Angle) -> Option<(Angle, Angle)> {
        let (ra, de) = self.equatorial();
        let ar = Angle::acos(de.sin() / lati.cos());
        let h = Angle::acos(-lati.tan() * de.tan());
        if h.radians().is_nan() || ar.radians().is_nan() {
            return None;
        }
        let lsts = (ra - h - longi).ungst(date);
        let lstr = (ra + h - longi).ungst(date);
        Some((lsts, lstr))
    }

    /// (Roughly) Accounts for precession in coordinates.
    pub fn precess(self, epoch: Date, d: Date) -> Self {
        let (ra, de) = self.equatorial();
        let diff = (d.julian() - epoch.julian()) / 365.25;
        let m =
            Angle::from_clock(0, 0, 3.07234) + Angle::from_clock(0, 0, 0.00186) * epoch.centuries();
        let n = Angle::from_degminsec(0, 0, 20.0468)
            + Angle::from_degminsec(0, 0, 0.0085) * epoch.centuries();
        let deltara = m.degrees() + n.degrees() * ra.sin() * de.tan();
        let deltade = n.degrees() * ra.cos();
        Coord::from_equatorial(
            ra + Angle::from_degrees(deltara * diff),
            de + Angle::from_degrees(deltade * diff),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Many of these tests do not conform with data you can pull out of stellarium/other tools, they are correct nonetheless.
    // * How do you know?: By personal confirmation of the result data with other resources
    // * Then how can I use these functions?: In conjunction with the functions for correction (procession, nutation, abberation, refraction)
    // Note that even without correction, these tests are almost always within 16' (half the moons diameter).
    #[test]
    fn test_horiz() {
        let arcturus = Coord::from_equatorial(
            Angle::from_clock(14, 16, 50.0),
            Angle::from_degminsec(19, 02, 50.1),
        );
        let sirius = Coord::from_equatorial(
            Angle::from_clock(6, 46, 13.1),
            Angle::from_degminsec(-16, 45, 06.8),
        );
        assert_eq!(
            arcturus.horizon(
                Date::from_calendar(2025, 3, 10, Angle::default()),
                Angle::from_clock(19, 52, 25.0),
                Angle::from_degrees(55.47885),
                Angle::from_degrees(133.94531)
            ),
            (
                Angle::from_degminsec(219, 35, 16.2),
                Angle::from_degminsec(48, 24, 46.1)
            )
        );
        assert_eq!(
            sirius.horizon(
                Date::from_calendar(2025, 3, 7, Angle::default()),
                Angle::from_clock(23, 36, 52.0),
                Angle::from_degrees(5.0),
                Angle::from_degrees(-1.0)
            ),
            (
                Angle::from_degminsec(249, 20, 18.2),
                Angle::from_degminsec(29, 28, 54.8)
            )
        );
        assert_eq!(
            sirius.horizon(
                Date::from_calendar(2025, 3, 11, Angle::default()),
                Angle::from_clock(2, 0, 0.0),
                Angle::from_degrees(44.8714),
                Angle::from_degrees(-93.20801)
            ),
            (
                Angle::from_degminsec(184, 42, 2.3),
                Angle::from_degminsec(29, 45, 27.2)
            )
        );
        assert_eq!(
            Coord::from_horizon(
                Angle::from_degminsec(184, 42, 2.3),
                Angle::from_degminsec(29, 45, 27.2),
                Date::from_calendar(2025, 3, 11, Angle::default()),
                Angle::from_clock(2, 0, 0.0),
                Angle::from_degrees(44.8714),
                Angle::from_degrees(-93.20801)
            ),
            sirius
        );
        assert_eq!(sirius.dist(arcturus), Angle::from_degminsec(115, 55, 5.17));
    }

    #[test]
    fn test_riseset() {
        let c = Coord::from_equatorial(
            Angle::from_clock(23, 39, 20.0),
            Angle::from_degminsec(21, 42, 00.0),
        );
        assert_eq!(
            c.riseset(
                Date::from_calendar(1980, 8, 24, Angle::default()),
                Angle::from_degrees(30.0),
                Angle::from_degrees(64.0)
            ),
            Some((Angle::from_clock(14, 18, 9.0), Angle::from_clock(4, 6, 5.0)))
        );
        assert_eq!(
            c.riseset(
                Date::from_calendar(1980, 8, 24, Angle::default()),
                Angle::from_degrees(-85.0),
                Angle::from_degrees(0.0),
            ),
            None
        );
    }

    #[test]
    fn test_ecliptic() {
        let star1 = Coord::from_equatorial(
            Angle::from_clock(9, 34, 53.6),
            Angle::from_degminsec(19, 32, 14.2),
        );
        assert_eq!(
            star1.ecliptic(Date::from_calendar(1950, 0, 1, Angle::default())),
            (
                Angle::from_degminsec(139, 41, 10.0),
                Angle::from_degminsec(4, 52, 31.0)
            )
        );
        assert_eq!(
            Coord::from_ecliptic(
                Angle::from_degminsec(139, 41, 10.0),
                Angle::from_degminsec(4, 52, 31.0),
                Date::from_calendar(1950, 0, 1, Angle::default())
            ),
            star1
        );
    }

    #[test]
    fn test_precession() {
        let star1 = Coord::from_equatorial(
            Angle::from_clock(3, 8, 10.6),
            Angle::from_degminsec(40, 57, 20.2),
        );
        /*assert_eq!(
            star1.precess(
                Date::from_calendar(2000, 1, 1, Angle::default()),
                Date::from_calendar(2024, 04, 04, Angle::default())
            ),
            star1
        );*/
    }
}
