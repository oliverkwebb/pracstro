//! # Time, Date, and Angle Handling
//!
//! This module contains functions for the handling and conversion of Times, Dates, and Angles.
//!
//! This data can be represented in two types:
//! - The [`Period`] type, which represents anything modulo arithmetic should be used to handle
//! - The [`Date`] type, which represents an instant in continuous time
//!
//! ```rust
//! use pracstro::*;
//! time::Date::from_calendar(2024, 06, 30.0).julian(); // Gets the julian date at 2024-06-30
//! ```

use std::fmt;
use std::ops::{Add, Sub};

/// True Modulus Operation, Least Positive Residue
fn lpr(x: f64, y: f64) -> f64 {
    let z = x % y;
    match z < 0.0 {
        true => z + y,
        false => z,
    }
}

/**
Angles and Time are the most prominent use for this type

| Property          | To Method               | From Method                  |
|-------------------|-------------------------|------------------------------|
| Degrees (Decimal) | [`Period::degrees()`]   | [`Period::from_degrees()`]   |
| Radians           | [`Period::radians()`]   | [`Period::from_radians()`]   |
| Turns (\[0,1\])   | [`Period::turns()`]     | [`Period::from_turns()`]     |
| Hours (Decimal)   | [`Period::decimal()`]   | [`Period::from_decimal()`]   |
| Clock Time        | [`Period::clock()`]     | [`Period::from_clock()`]     |
| Degrees (DMS)     | [`Period::degminsec()`] | [`Period::from_degminsec()`] |
| Sine              | [`Period::sin()`]       | [`Period::asin()`]           |
| Cosine            | [`Period::cos()`]       | [`Period::acos()`]           |
| Tangent           | [`Period::tan()`]       | [`Period::atan2()`]          |

Additional Methods:
* Latitude displaying: [`Period::to_latitude()`]
* GST Correction: [`Period::gst()`] and [`Period::ungst()`]
* Inverse of angle: [`Period::inverse()`]
*/
#[derive(Clone, Copy, Default)]
pub struct Period(f64);
impl Period {
    /// Returns the angle as radians.
    ///
    /// This is the only function that should directly access the fields of the type.
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_degrees(180.0).radians(); // Pi
    /// ```
    pub fn radians(self) -> f64 {
        self.0
    }
    /// Constructs a angle from radians, reducing it to the range of \[0, 2*PI\].
    ///
    /// This is the one of the two only functions that directly access [`Period`].
    /// Reduction to the desired range is done with Least Positive Residue instead of the remainder operator.
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_radians(std::f64::consts::PI).degrees(); // 180.0
    /// ```
    pub fn from_radians(x: f64) -> Self {
        Period(lpr(x, std::f64::consts::TAU))
    }
    /// Converts angles internally so that formatting them as latitudes makes sense
    ///
    /// This should be used directly in conjunction with a formatting function.
    /// The main use case for this is in coordinate handling where angles of latitude can be negative.
    /// This is the only other function that should be able to directly access [`Period`] than [`Period::from_radians()`].
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_degrees(-25.0).degrees(); // 335.0
    /// Period::from_degrees(-25.0).to_latitude().degrees(); // -25.0
    /// ```
    pub fn to_latitude(self) -> Self {
        match self.radians() > std::f64::consts::PI {
            true => Period(self.radians() - std::f64::consts::TAU),
            false => self,
        }
    }

    /// Returns the angle as fractional degrees.
    ///
    /// A wrapper around [`f64::to_degrees()`]
    ///
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_degrees(-25.0).degrees(); // 335.0
    /// ```
    pub fn degrees(self) -> f64 {
        self.radians().to_degrees()
    }
    /// Constructs an angle from fractional degrees.
    ///
    /// A wrapper around [`f64::to_radians()`]
    ///
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_degrees(-25.0).degrees(); // 335.0
    /// ```
    pub fn from_degrees(x: f64) -> Self {
        Period::from_radians(x.to_radians())
    }

    /// Returns the angle in the range between 0 and 1 (i.e. turns)
    ///
    /// Used in the calculation for the illuminated fraction of planets
    pub fn turns(self) -> f64 {
        self.degrees() / 360.0
    }
    /// Constructs an angle from a value between 0 and 1
    ///
    /// Used in the calculation for the illuminated fraction of planets
    pub fn from_turns(x: f64) -> Self {
        Self::from_degrees(x * 360.0)
    }

    /// Returns the angle in fractional number of hours
    ///
    /// A wrapper around [`Period::degrees()`] and one division by 15. Since one hour is 15 degrees of rotation
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_degrees(120.0).decimal(); // 8.0
    /// ```
    pub fn decimal(self) -> f64 {
        self.degrees() / 15.0 // 1 hour <-> 15 degrees, 360/24 = 15
    }
    /// Constructs an angle from a fractional number of hours
    ///
    /// A wrapper around [`Period::from_degrees()`] and one multiplication by 15. Since one hour is 15 degrees of rotation
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_decimal(8.0).degrees(); // 120.00
    /// ```
    pub fn from_decimal(x: f64) -> Self {
        Period::from_degrees(x * 15.0)
    }

    /// Returns (Hour, Minute, Second) of a time/angle
    ///
    /// Used in hour-angle displays for some coordinate systems, and in times.
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_decimal(8.0).clock(); // (8, 0, 0.0)
    /// ```
    pub fn clock(self) -> (u8, u8, f64) {
        let y = self.decimal();
        (
            y.trunc() as u8,
            (y.fract() * 60.0).trunc() as u8,
            (y.fract() * 60.0).fract() * 60.0,
        )
    }
    /// Constructs an angle out of an hour, minute, and second
    ///
    /// Used in hour-angle displays for some coordinate systems, and in times.
    /// ```
    /// use pracstro::time::Period;
    /// Period::from_clock(8, 0, 0.0).decimal(); // 8.0
    /// ```
    pub fn from_clock(h: u8, m: u8, s: f64) -> Self {
        Period::from_decimal((h as f64) + (((m as f64) + (s / 60.0)) / 60.0))
    }

    /// Converts an angle to a degree with arcminutes and arcseconds
    pub fn degminsec(self) -> (i16, u8, f64) {
        let y = self.degrees();
        (
            y.trunc() as i16,
            (y.fract() * 60.0).trunc() as u8,
            (y.fract() * 60.0).fract() * 60.0,
        )
    }
    /// Identical to from_clock in math
    pub fn from_degminsec(d: i16, m: u8, s: f64) -> Self {
        Period::from_degrees((d as f64) + (((m as f64) + (s / 60.0)) / 60.0))
    }

    /// Handles the discontinuity created by the orbit of the earth as compared to its rotation.
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn gst(self, date: Date) -> Self {
        let t = date.centuries();
        let t0 = lpr(
            6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
            24.0,
        );
        Period::from_decimal(lpr(t0 + (self.decimal() * 1.002737909), 24.0))
    }
    /// Handles the discontinuity created by the orbit of the earth as compared to its rotation.
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn ungst(self, date: Date) -> Self {
        let t = date.centuries();
        let t0 = lpr(
            6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
            24.0,
        );
        Period::from_decimal(lpr(self.decimal() - t0, 24.0) * 0.9972695663)
    }

    /// Gets the hour angle from the angle as if it were a right ascension, and vice versa.
    ///
    /// Used in horizontal coordinates.
    pub fn hourangle_rightas(self, date: Date, time: Period, longi: Period) -> Self {
        time.gst(date) + longi - self
    }

    /// Sine of Period
    pub fn sin(self) -> f64 {
        self.radians().sin()
    }
    /// Cosine of Period
    pub fn cos(self) -> f64 {
        self.radians().cos()
    }
    /// Tangent of Period
    pub fn tan(self) -> f64 {
        self.radians().tan()
    }
    /// Period from Arcsine
    pub fn asin(x: f64) -> Self {
        Period::from_radians(x.asin())
    }
    /// Period from Arccosine
    pub fn acos(x: f64) -> Self {
        Period::from_radians(x.acos())
    }
    /// Period from 2-argument Arctangent
    pub fn atan2(x: f64, y: f64) -> Self {
        Period::from_radians(x.atan2(y))
    }
    /// Reverses the angle
    pub fn inverse(self) -> Self {
        Period::from_degrees(360.0 - self.degrees())
    }
}
/// Used in testing
impl fmt::Debug for Period {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (d, m, s) = self.degminsec();
        // Hubble has a resolution of 0.1", this is more than sufficent
        write!(f, "{}°{}'{:.2}\"", d, m, s)
    }
}
/// Does not check if arcseconds are equal
impl PartialEq for Period {
    fn eq(&self, other: &Self) -> bool {
        let (d, m, _) = self.degminsec();
        let (d2, m2, _) = other.degminsec();
        d == d2 && m == m2
    }
}
impl Add<Period> for Period {
    type Output = Period;
    /// Addition, For timezones and LST
    fn add(self, x: Self) -> Self {
        Period::from_radians(self.radians() + x.radians())
    }
}
impl Sub<Period> for Period {
    type Output = Period;
    /// Subtraction, For timezones and LST
    fn sub(self, x: Self) -> Self {
        Period::from_radians(self.radians() - x.radians())
    }
}

/**
Continuous Instant in Time

| Property          | To Method             | From Method                |
|-------------------|-----------------------|----------------------------|
| Julian Day        | [`Date::julian()`]    | [`Date::from_julian()`]    |
| Calendar          | [`Date::calendar()`]  | [`Date::from_calendar()`]  |
| Unix Time         | [`Date::unix()`]      | [`Date::from_unix()`]      |

Additional Methods
* Get the current time: [`Date::now()`]
*/
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Date(f64);
impl Date {
    /// Direct interface to type
    ///
    /// This is the only function that should directly read the fields of the type
    pub fn julian(self) -> f64 {
        self.0
    }
    /// Direct interface to type
    ///
    /// This is the only function that should directly write the fields of the type
    pub fn from_julian(x: f64) -> Self {
        Date(x)
    }

    /// Returns Julian Centuries since 1900.
    ///
    /// Used heavily in astronomical estimation of things that change slowly
    pub fn centuries(self) -> f64 {
        (self.julian() - 2451545.0) / 36525.0
    }

    /// Returns Year, Month, Day (time is Period::from_decimal(day.fract()))
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn calendar(self) -> (u32, u8, f64) {
        let j = self.julian() + 0.5;
        let (i, f) = (j.trunc(), j.fract());

        let b = if i > 2_299_160.0 {
            let a = ((i - 1867216.25) / 36524.25).trunc();
            i + 1.0 + a - (a / 4.0).trunc()
        } else {
            1.0
        };

        let c = b + 1524.0;
        let d = ((c - 122.2) / 365.25).trunc();
        let e = (365.25 * d).trunc();
        let g = ((c - e) / 30.6001).trunc();

        let m = if g < 13.5 { g - 1.0 } else { g - 13.0 };
        let y = if m > 2.5 { d - 4716.0 } else { d - 4715.0 };
        let d = c - e + f - (30.6001 * g).trunc();

        (y as u32, m as u8, d)
    }
    /// Takes Year, Month, and Day
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn from_calendar(y: u64, m: u8, d: f64) -> Self {
        let (mut year, mut month, day): (f64, f64, f64) = (y as f64, m as f64, d);
        if month < 3.0 {
            year -= 1.0;
            month += 12.0;
        }
        let (yp, mp) = (year, month);

        let b = if year >= 1582.0 {
            let a = (yp / 100.0).trunc();
            2.0 - a + (a / 4.0).trunc()
        } else {
            0.0
        };

        let c = (if yp < 0.0 {
            (365.25 * yp) - 0.75
        } else {
            365.25 * yp
        })
        .trunc();

        Date::from_julian(b + c + (30.6001 * (mp + 1.0)).trunc() + day + 1_720_994.5)
    }
    /// Gets the time of day in a current calendar date
    pub fn time(self) -> Period {
    	Period::from_decimal(self.calendar().2.fract()*24.0)
    }

    /// Interface for unix time, Does not correct for the 1582 Julain/Gregorian split
    pub fn unix(self) -> f64 {
        (self.julian() - 2440587.5) * 86400.0
    }
    /// Interface for unix time, Does not correct for the 1582 Julain/Gregorian split
    pub fn from_unix(t: f64) -> Self {
        Date::from_julian((t / 86400.0) + 2440587.5)
    }

    /// Gets the current date
    ///
    /// Since this function relies on SystemTime and duration_since, it does not work for dates before 1970.
    pub fn now() -> Self {
        use std::time::{SystemTime, UNIX_EPOCH};
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Expected pre-1970-01-01 date")
            .as_secs() as f64;
        Date::from_unix(now)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lpr() {
        assert_eq!(lpr(-1.0, 5.0), 4.0);
        assert_eq!(lpr(7.0, 5.0), 2.0);
    }

    #[test]
    fn test_julian() {
        assert_eq!(
            Date::from_julian(2_446_113.75),
            Date::from_calendar(1985, 2, 17.25)
        );
        assert_eq!(Date::from_julian(2_446_113.75).calendar(), (1985, 2, 17.25));
        assert_eq!(Date::from_calendar(1967, 04, 12.6).time().decimal(), 14.400000002235174);
    }

    #[test]
    fn test_decimalhrs() {
        assert_eq!(
            Period::from_clock(18, 31, 27.0),
            Period::from_decimal(18.52417)
        );
        assert_eq!(Period::from_decimal(11.75), Period::from_clock(11, 45, 0.0));
    }

    #[test]
    fn test_gst() {
        assert_eq!(
            Period::from_decimal(14.614_353).gst(Date::from_julian(2_444_351.5)),
            Period::from_decimal(4.668119549708194)
        );
        assert_eq!(
            Period::from_decimal(4.668119549708194).ungst(Date::from_julian(2_444_351.5)),
            Period::from_decimal(14.614352994461141)
        );
    }

    #[test]
    fn test_lst() {
        assert_eq!(
            Period::from_decimal(4.668_119).add(Period::from_degrees(-64.0)),
            Period(0.10509997509720659)
        );
        assert_eq!(
            Period::from_decimal(0.401_453).sub(Period::from_degrees(-64.0)),
            Period(1.2221108709065023)
        );
    }
    #[test]
    fn test_lati() {
        assert_eq!(
            Period::from_degrees(-25.0).to_latitude().degrees(),
            -25.00000000000002 // God I hate floating point
        );
    }

    #[test]
    fn test_turn() {
        assert_eq!(Period::from_turns(0.5), Period::from_degrees(180.0));
        assert_eq!(Period::from_degrees(180.0), Period::from_turns(0.5));
    }
}
