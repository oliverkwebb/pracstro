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
	pub fn celestial(self) -> (Period, Period) {
		(self.0, self.1)
	}
	pub fn from_celestial(x: Period, y: Period) -> Self {
		Coord(x, y)
	}

	pub fn equatorial(self, date: Date, time: Period, longi: Period) -> (Period, Period) {
		let (ra, de) = self.celestial();
		(time.gst(date).add(longi).sub(ra), de)
	}
	pub fn from_equatorial(ha: Period, de: Period, date: Date, time: Period, longi: Period) -> Self {
		Coord::from_celestial(time.gst(date).add(longi).sub(ha), de)
	}

	pub fn horizon(self, date: Date, time: Period, longi: Period, lati: Period) -> (Period, Period) {
		let (ha, de) = self.equatorial(date, time, longi);
		let alt = Period::asin(de.sin()*lati.sin() + de.cos()*lati.cos()*ha.cos());
		let azip = Period::acos(de.sin() - lati.sin() * alt.sin() / lati.cos() * alt.cos());
		let azi;
		match ha.sin() > 0.0 {
			true => azi = azip,
			false => azi = Period::from_degrees(360.0 - azip.degrees()),
		};
		(azi, alt)
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_ra_ha() {
		let star1 = Coord::from_celestial(Period::from_clock(18, 32, 21.0), Period::from_degrees(23.4));
		assert_eq!(
			star1.equatorial(Date::from_calendar(1980, 4, 22.0), Period::from_clock(14, 36, 51.67), Period::from_degrees(-64.0)).0,
			Period::from_radians(1.5325395556005414) // Period::from_clock(5, 51, 44.0)
		);
		assert_eq!(
			Coord::from_equatorial(Period::from_clock(5, 51, 44.0), star1.celestial().1 ,Date::from_calendar(1980, 4, 22.0), Period::from_clock(14, 36, 51.67), Period::from_degrees(-64.0)).0,
			Period::from_radians(4.853000580733089) // Period::from_clock(18, 32, 21.0)
		);
	}

	#[test]
	fn test_horiz() {
		let star1 = Coord::from_celestial(Period::from_clock(18, 32, 21.0), Period::from_degrees(23.25));
		assert_eq!(
			star1.horizon(Date::from_calendar(1980, 4, 22.0), Period::from_clock(14, 36, 51.67), Period::from_degrees(-64.0), Period::from_degrees(52.0)),
			(Period::from_radians(1.5776278343270889), Period::from_radians(0.3391627108196441))
		);
	}
}