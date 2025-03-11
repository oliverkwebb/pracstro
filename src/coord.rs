use crate::time::*;

/// Pair of period values, Representing "How far up" and "How far round"
///
/// Base Value is right ascension, and declination
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Coord(Period, Period);
/// Interfaces:
/// - Celestial (Right Ascension, Declination)
/// - Equatorial (Hour Angle, Declination)
///	- Horizon (Azimuth, Altitude)
/// - Ecliptic (Ecliptic Latitude, Ecliptic Longitude)
/// The book also specifies Galactic Coordinates, but never uses them
impl Coord {
    /// Right Ascension and Declination
    pub fn celestial(self) -> (Period, Period) {
        (self.0, self.1)
    }
    pub fn from_celestial(x: Period, y: Period) -> Self {
        Coord(x, y)
    }

    /// Hour angle and Declination, dependent on longitude and time
    pub fn equatorial(self, date: Date, time: Period, longi: Period) -> (Period, Period) {
        let (ra, de) = self.celestial();
        (time.gst(date).add(longi).sub(ra), de)
    }
    pub fn from_equatorial(
        ha: Period,
        de: Period,
        date: Date,
        time: Period,
        longi: Period,
    ) -> Self {
        Coord::from_celestial(time.gst(date).add(longi).sub(ha), de)
    }

    pub fn horizon(
        self,
        date: Date,
        time: Period,
        lati: Period,
        longi: Period,
    ) -> (Period, Period) {
        let (ha, de) = self.equatorial(date, time, longi);
        let alt = Period::asin(de.sin() * lati.sin() + de.cos() * lati.cos() * ha.cos());
        let azip = Period::acos((de.sin() - lati.sin() * alt.sin()) / (lati.cos() * alt.cos()));
        let azi = match ha.sin() < 0.0 {
            true => azip,
            false => Period::from_degrees(360.0 - azip.degrees()),
        };
        (azi, alt)
    }
    pub fn from_horizon(
        azi: Period,
        alt: Period,
        date: Date,
        time: Period,
        lati: Period,
        longi: Period,
    ) -> Self {
        let de = Period::asin(alt.sin() * lati.sin() + alt.cos() * lati.cos() * azi.cos());
        let hap = Period::acos((alt.sin() - lati.sin() * de.sin()) / (lati.cos() * de.cos()));
        let ha = match azi.sin() < 0.0 {
            true => hap,
            false => Period::from_degrees(360.0 - hap.degrees()),
        };
        Coord::from_equatorial(ha, de, date, time, longi)
    }

    /// Returns the angle between two objects
    pub fn dist(self, from: Self) -> Period {
        let ((a1, d1), (a2, d2)) = (self.celestial(), from.celestial());
        Period::acos(d1.sin() * d2.sin() + d1.cos() * d2.cos() * a1.sub(a2).cos())
    }
    /// Returns (Rise, Set) UT, This function will fail for locations in the sky that never appear over the horizon
    pub fn riseset(self, date: Date, lati: Period, longi: Period) -> Option<(Period, Period)> {
        let (ra, de) = self.celestial();
        let ar = Period::acos(de.sin() / lati.cos());
        let h = Period::acos(-lati.tan() * de.tan());
        if h.radians().is_nan() || ar.radians().is_nan() {
            return None;
        }
        let lsts = ra.sub(h).sub(longi).ungst(date);
        let lstr = ra.add(h).sub(longi).ungst(date);
        Some((lsts, lstr))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Many of these tests do not conform with data you can pull out of stellarium/other tools, they are correct nonetheless.
    // * How do you know?: By personal confirmation of the result data with other resources
    // * Then how can I use these functions?: In conjunction with the functions for correction (procession, nutation, abberation, refraction)
    // Note that even without correction, these tests are almost always within 15' (half the moons diameter).
    #[test]
    fn test_ra_ha() {
        let star1 =
            Coord::from_celestial(Period::from_clock(18, 32, 21.0), Period::from_degrees(23.4));
        assert_eq!(
            star1
                .equatorial(
                    Date::from_calendar(1980, 4, 22.0),
                    Period::from_clock(14, 36, 51.67),
                    Period::from_degrees(-64.0)
                )
                .0,
            Period::from_radians(1.5325395556005414) // Period::from_clock(5, 51, 44.0)
        );
        assert_eq!(
            Coord::from_equatorial(
                Period::from_clock(5, 51, 44.0),
                star1.celestial().1,
                Date::from_calendar(1980, 4, 22.0),
                Period::from_clock(14, 36, 51.67),
                Period::from_degrees(-64.0)
            )
            .0,
            Period::from_radians(4.853000580733089) // Period::from_clock(18, 32, 21.0)
        );
    }

    #[test]
    fn test_horiz() {
        let arcturus = Coord::from_celestial(
            Period::from_clock(14, 16, 50.0),
            Period::from_degminsec(19, 02, 50.1),
        );
        let sirius = Coord::from_celestial(
            Period::from_clock(6, 46, 13.1),
            Period::from_degminsec(-16, 45, 06.8),
        );
        assert_eq!(
            arcturus.horizon(
                Date::from_calendar(2025, 3, 10.0),
                Period::from_clock(19, 52, 25.0),
                Period::from_degrees(55.47885),
                Period::from_degrees(133.94531)
            ),
            (
                Period::from_degminsec(219, 42, 17.3),
                Period::from_degminsec(48, 21, 43.1)
            )
        );
        assert_eq!(
            sirius.horizon(
                Date::from_calendar(2025, 3, 7.0),
                Period::from_clock(23, 36, 52.0),
                Period::from_degrees(5.0),
                Period::from_degrees(-1.0)
            ),
            (
                Period::from_degminsec(249, 17, 18.2),
                Period::from_degminsec(29, 37, 54.8)
            )
        );
        assert_eq!(
            sirius.horizon(
                Date::from_calendar(2025, 3, 11.0),
                Period::from_clock(2, 0, 0.0),
                Period::from_degrees(44.8714),
                Period::from_degrees(-93.20801)
            ),
            (
                Period::from_degminsec(184, 45, 12.1),
                Period::from_degminsec(29, 45, 9.25)
            )
        );
        /* assert_eq!(
            Coord::from_horizon(
                Period::from_degminsec(184, 45, 12.1),
                Period::from_degminsec(29, 45, 9.25),
                Date::from_calendar(2025, 3, 11.0),
                Period::from_clock(2, 0, 0.0),
                Period::from_degrees(44.8714),
                Period::from_degrees(-93.20801)
            ),
            sirius
        );*/
        assert_eq!(sirius.dist(arcturus), Period::from_degminsec(115, 46, 0.5));
    }

    #[test]
    fn test_riseset() {
        let c = Coord::from_celestial(
            Period::from_clock(23, 39, 20.0),
            Period::from_degminsec(21, 42, 00.0),
        );
        assert_eq!(
            c.riseset(
                Date::from_calendar(1980, 8, 24.0),
                Period::from_degrees(30.0),
                Period::from_degrees(64.0)
            )
            .unwrap(),
            (
                Period::from_degminsec(214, 27, 17.0),
                Period::from_degminsec(61, 26, 22.0)
            )
        );
    }
}
