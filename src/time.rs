//! # Time, Date, and Angle Handling
//!
//! This module contains functions for the handling and conversion of Times, Dates, and Angles.
//!
//! This data can be represented in two types:
//! - The [`Angle`] type, which represents anything modulo arithmetic should be used to handle
//! - The [`Date`] type, which represents an instant in continuous time
//!
//! ```rust
//! # use pracstro::*;
//! time::Date::from_calendar(2024, 06, 30, time::Angle::from_clock(16, 30, 0.0)).julian(); // Gets the julian date at 2024-06-30T16:30:00Z
//! ```

use std::f64::consts::{PI, TAU};
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

/**
Angles and Time are the most prominent use for this type

| Property          | To Method               | From Method                  |
|-------------------|-------------------------|------------------------------|
| Degrees (Decimal) | [`Angle::degrees()`]   | [`Angle::from_degrees()`]   |
| Radians           | [`Angle::radians()`]   | [`Angle::from_radians()`]   |
| Turns (\[0,1\])   | [`Angle::turns()`]     | [`Angle::from_turns()`]     |
| Hours (Decimal)   | [`Angle::decimal()`]   | [`Angle::from_decimal()`]   |
| Clock Time        | [`Angle::clock()`]     | [`Angle::from_clock()`]     |
| Degrees (DMS)     | [`Angle::degminsec()`] | [`Angle::from_degminsec()`] |
| Sine              | [`Angle::sin()`]       | [`Angle::asin()`]           |
| Cosine            | [`Angle::cos()`]       | [`Angle::acos()`]           |
| Tangent           | [`Angle::tan()`]       | [`Angle::atan2()`]          |

Additional Methods:
* Latitude displaying: [`Angle::to_latitude()`]
* Inverse of angle: [`Angle::inverse()`]
* GST Correction: [`Angle::gst()`] and [`Angle::ungst()`]
* Approx. Atmosphereic Refraction: [`Angle::refract()`] and [`Angle::refractdelta()`]
*/
#[derive(Clone, Copy, Default, PartialOrd)]
pub struct Angle(f64);
impl Angle {
    /// Returns the angle as radians.
    ///
    /// This is the only function that should directly access the fields of the type.
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_degrees(180.0).radians(); // Pi
    /// ```
    pub const fn radians(self) -> f64 {
        self.0
    }
    /// Constructs a angle from radians, reducing it to the range of \[0, 2*PI\].
    ///
    /// This is the one of the two only functions that directly access [`Angle`].
    /// Reduction to the desired range is done with Least Positive Residue instead of the remainder operator.
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_radians(std::f64::consts::PI).degrees(); // 180.0
    /// ```
    pub const fn from_radians(x: f64) -> Self {
        const fn lpr(x: f64, y: f64) -> f64 {
            let z = x % y;
            z + if z < 0.0 { y } else { 0.0 }
        }
        Angle(lpr(x, TAU))
    }
    /// Converts angles internally so that formatting them as latitudes makes sense
    ///
    /// This should be used directly in conjunction with a formatting function.
    /// The main use case for this is in coordinate handling where angles of latitude can be negative.
    /// This is the only other function that should be able to directly access [`Angle`] than [`Angle::from_radians()`].
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_degrees(-25.0).degrees(); // 335.0
    /// Angle::from_degrees(-25.0).to_latitude().degrees(); // -25.0
    /// ```
    pub const fn to_latitude(self) -> Self {
        match self.radians() > PI {
            true => Angle(self.radians() - TAU),
            false => self,
        }
    }

    /// Returns the angle as fractional degrees.
    ///
    /// A wrapper around [`f64::to_degrees()`]
    ///
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_degrees(-25.0).degrees(); // 335.0
    /// ```
    pub const fn degrees(self) -> f64 {
        self.radians().to_degrees()
    }
    /// Constructs an angle from fractional degrees.
    ///
    /// A wrapper around [`f64::to_radians()`]
    ///
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_degrees(-25.0).degrees(); // 335.0
    /// ```
    pub const fn from_degrees(x: f64) -> Self {
        Angle::from_radians(x.to_radians())
    }

    /// Returns the angle in the range between 0 and 1 (i.e. turns)
    ///
    /// Used in the calculation for the illuminated fraction of planets
    pub const fn turns(self) -> f64 {
        self.radians() / TAU
    }
    /// Constructs an angle from a value between 0 and 1
    ///
    /// Used in the calculation for the illuminated fraction of planets
    pub const fn from_turns(x: f64) -> Self {
        Self::from_radians(x * TAU)
    }

    /// Returns the angle in fractional number of hours
    ///
    /// A wrapper around [`Angle::degrees()`] and one division by 15. Since one hour is 15 degrees of rotation
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_degrees(120.0).decimal(); // 8.0
    /// ```
    pub const fn decimal(self) -> f64 {
        self.degrees() / 15.0 // 1 hour <-> 15 degrees, 360/24 = 15
    }
    /// Constructs an angle from a fractional number of hours
    ///
    /// A wrapper around [`Angle::from_degrees()`] and one multiplication by 15. Since one hour is 15 degrees of rotation
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_decimal(8.0).degrees(); // 120.00
    /// ```
    pub const fn from_decimal(x: f64) -> Self {
        Angle::from_degrees(x * 15.0)
    }

    /// Returns (Hour, Minute, Second) of a time/angle
    ///
    /// Used in hour-angle displays for some coordinate systems, and in times.
    /// ```
    /// # use pracstro::time::Angle;
    /// Angle::from_decimal(8.0).clock(); // (8, 0, 0.0)
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
    /// use pracstro::time::Angle;
    /// Angle::from_clock(8, 0, 0.0).decimal(); // 8.0
    /// ```
    pub const fn from_clock(h: u8, m: u8, s: f64) -> Self {
        Angle::from_decimal((h as f64) + (((m as f64) + (s / 60.0)) / 60.0))
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
    pub const fn from_degminsec(d: i16, m: u8, s: f64) -> Self {
        Angle::from_degrees((d as f64) + (((m as f64) + (s / 60.0)) / 60.0))
    }

    /// Handles the discontinuity created by the orbit of the earth as compared to its rotation.
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn gst(self, date: Date) -> Self {
        let t = date.centuries();
        Angle::from_decimal(
            6.697374558
                + (2400.051336 * t)
                + (0.000025862 * (t * t))
                + (self.decimal() * 1.002737909),
        )
    }
    /// Handles the discontinuity created by the orbit of the earth as compared to its rotation.
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn ungst(self, date: Date) -> Self {
        let t = date.centuries();
        let t0 = Angle::from_decimal(6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)));
        (self - t0) * 0.9972695663
    }

    /// Sine of Angle
    pub fn sin(self) -> f64 {
        self.radians().sin()
    }
    /// Cosine of Angle
    pub fn cos(self) -> f64 {
        self.radians().cos()
    }
    /// Tangent of Angle
    pub fn tan(self) -> f64 {
        self.radians().tan()
    }
    /// Angle from Arcsine
    pub fn asin(x: f64) -> Self {
        Angle::from_radians(x.asin())
    }
    /// Angle from Arccosine
    pub fn acos(x: f64) -> Self {
        Angle::from_radians(x.acos())
    }
    /// Angle from 2-argument Arctangent
    pub fn atan2(x: f64, y: f64) -> Self {
        Angle::from_radians(x.atan2(y))
    }
    /// Reverses the angle
    pub fn inverse(self) -> Self {
        Angle::from_radians(TAU - self.radians())
    }

    /// Calculates the approximate atmospheric refraction
    ///
    /// In reality, this is an complex calculation dependent on factors such as temperature and pressure. But it can be
    /// Approximated to a reasonable extent assuming some factors.
    ///
    /// From <https://www.celestialprogramming.com/snippets/atmosphericrefraction.html>
    pub fn refractdelta(self) -> Self {
        Angle::from_degminsec(
            0,
            0,
            (1.02
                / (self.degrees() + (10.3 / (self.degrees() + 5.11)))
                    .to_radians()
                    .tan())
                * 60.0,
        )
    }
    /// Accounts for atmospheric refraction
    ///
    /// This does not calculate refraction on altitudes under the horizon
    ///
    /// This should only be used on the altitude value of a horizontal coordinate.
    pub fn refract(self) -> Self {
        if self.to_latitude().degrees() > 0.0 {
            self + self.refractdelta()
        } else {
            self
        }
    }
}
/// Used in testing
impl fmt::Debug for Angle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (d, m, s) = self.degminsec();
        // Hubble has a resolution of 0.1", this is more than sufficent
        write!(f, "{}°{}'{:.2}\"", d, m, s)
    }
}
/// Does not check if arcseconds are equal
impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        let (d, m, _) = self.degminsec();
        let (d2, m2, _) = other.degminsec();
        d == d2 && m == m2
    }
}
impl Add<Angle> for Angle {
    type Output = Angle;
    /// Addition, For timezones and LST
    fn add(self, x: Self) -> Self {
        Angle::from_radians(self.radians() + x.radians())
    }
}
impl Sub<Angle> for Angle {
    type Output = Angle;
    /// Subtraction, For timezones and LST
    fn sub(self, x: Self) -> Self {
        Angle::from_radians(self.radians() - x.radians())
    }
}
impl Mul<f64> for Angle {
    type Output = Angle;
    /// Multiplication
    fn mul(self, x: f64) -> Self {
        Angle::from_radians(self.radians() * x)
    }
}
impl Div<f64> for Angle {
    type Output = Angle;
    /// Multiplication
    fn div(self, x: f64) -> Self {
        Angle::from_radians(self.radians() / x)
    }
}

/**
Continuous Instant in Time

| Property          | To Method             | From Method                |
|-------------------|-----------------------|----------------------------|
| Julian Day        | [`Date::julian()`]    | [`Date::from_julian()`]    |
| Calendar          | [`Date::calendar()`]  | [`Date::from_calendar()`]  |
| Unix Time         | [`Date::unix()`]      | [`Date::from_unix()`]      |
| Date/Time         | [`Date::time()`]      | [`Date::from_time()`]      |

Additional Methods
* Get the current time: [`Date::now()`]
* Julian Centuries since J2000: [`Date::centuries()`]
*/
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Date(f64);
impl Date {
    /// Direct interface to type
    ///
    /// This is the only function that should directly read the fields of the type
    pub const fn julian(self) -> f64 {
        self.0
    }
    /// Direct interface to type
    ///
    /// This is the only function that should directly write the fields of the type
    pub const fn from_julian(x: f64) -> Self {
        Date(x)
    }

    /// Returns Julian Centuries since 1900.
    ///
    /// Used heavily in astronomical estimation of things that change slowly
    pub const fn centuries(self) -> f64 {
        (self.julian() - 2451545.0) / 36525.0
    }

    /// Returns Year, Month, Day (time is Angle::from_decimal(day.fract()))
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, Although similar algorithms exist in other sources
    pub fn calendar(self) -> (i64, u8, u8, Angle) {
        let j = self.julian() + 0.5;
        let (i, f) = (j.trunc(), j.fract());

        let b = if i > 2_299_160.0 {
            let a = ((i - 1867216.25) / 36524.25).trunc();
            i + 1.0 + a - (a / 4.0).trunc()
        } else {
            i
        } + 1524.0;

        let d = ((b - 122.2) / 365.25).trunc();
        let e = (365.25 * d).trunc();
        let g = ((b - e) / 30.6001).trunc();

        let m = if g < 13.5 { g - 1.0 } else { g - 13.0 };
        let y = if m > 2.5 { d - 4716.0 } else { d - 4715.0 };
        let d = b - e + f - (30.6001 * g).trunc();

        (y as i64, m as u8, d as u8, Angle::from_turns(d.fract()))
    }
    /// Takes Year, Month, and Day
    ///
    /// Algorithm from Practical Astronomy with Your Calculator, although similar algorithms exist in other sources
    pub fn from_calendar(y: i64, m: u8, day: u8, t: Angle) -> Self {
        let (year, month) = if m < 3 { (y - 1, m + 12) } else { (y, m) };

        Date::from_julian(
            if year >= 1582 {
                2 - (year / 100) + (year / 400)
            } else {
                0
            } as f64
                + (365.25 * year as f64 - if year < 0 { 0.75 } else { 0.0 }).trunc()
                + (30.6001 * (month + 1) as f64).trunc()
                + day as f64
                + t.turns()
                + 1_720_994.5,
        )
    }
    /// Gets the time of day in a current calendar date
    pub fn time(self) -> Angle {
        self.calendar().3
    }
    /// Constructs a date out of a date and a time.
    pub fn from_time(d: Self, t: Angle) -> Self {
        let (y, m, d, _) = d.calendar();
        Self::from_calendar(y, m, d, t)
    }

    /// Interface for unix time, Does not correct for the 1582 Julain/Gregorian split
    pub const fn unix(self) -> f64 {
        (self.julian() - 2440587.5) * 86400.0
    }
    /// Interface for unix time, Does not correct for the 1582 Julain/Gregorian split
    pub const fn from_unix(t: f64) -> Self {
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

/// Calculate the date of Easter
pub fn easter(year: i32) -> (i32, i32) {
    let a = year % 19;
    let (b, c) = (year / 100, year % 100);
    let d = b / 4;
    let g = (b - (b + 8) / 25 + 1) / 3;
    let h = ((19 * a) + b - d - g + 15) % 30;
    let l = (32 + (2 * b % 4) + (2 * c / 4) - h - c % 4) % 7;
    let m = (a + (11 * h) + (22 * l)) / 451;
    let n = (h + l - (7 * m) + 114) / 31;
    let p = (h + l - (7 * m) + 114) % 31;
    (n, p + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_julian() {
        assert_eq!(
            Date::from_julian(2_446_113.75),
            Date::from_calendar(1985, 2, 17, Angle::from_decimal(6.0))
        );
        assert_eq!(
            Date::from_julian(2_446_113.75).calendar(),
            (1985, 2, 17, Angle::from_decimal(6.0))
        );
        assert_eq!(
            Date::from_calendar(1967, 04, 12, Angle::from_turns(0.6))
                .time()
                .decimal(),
            14.400000002235174
        );
    }

    #[test]
    fn test_calendar() {
        assert_eq!(
            Date::from_calendar(1000, 10, 20, Angle::default()).calendar(),
            (1000, 10, 20, Angle::default())
        );
        assert_eq!(
            Date::from_calendar(1000, 14, 20, Angle::default()).calendar(),
            (1001, 2, 20, Angle::default())
        );
    }

    #[test]
    fn test_decimalhrs() {
        assert_eq!(
            Angle::from_clock(18, 31, 27.0),
            Angle::from_decimal(18.52417)
        );
        assert_eq!(Angle::from_decimal(11.75), Angle::from_clock(11, 45, 0.0));
    }

    #[test]
    fn test_gst() {
        assert_eq!(
            Angle::from_clock(14, 36, 51.6).gst(Date::from_julian(2_444_351.5)),
            Angle::from_clock(4, 40, 5.23)
        );
        assert_eq!(
            Angle::from_clock(4, 40, 5.23).ungst(Date::from_julian(2_444_351.5)),
            Angle::from_clock(14, 36, 51.6)
        );
    }

    #[test]
    fn test_lati() {
        assert_eq!(
            Angle::from_degrees(-25.0).to_latitude().degrees(),
            -25.00000000000002 // God I hate floating point
        );
    }

    #[test]
    fn test_turn() {
        assert_eq!(Angle::from_turns(0.5), Angle::from_degrees(180.0));
        assert_eq!(Angle::from_degrees(180.0), Angle::from_turns(0.5));
    }

    #[test]
    fn test_easter() {
        assert_eq!(easter(2000), (4, 23));
        assert_eq!(easter(2024), (3, 31));
    }

    #[test]
    fn test_refract() {
        assert_eq!(
            Angle::from_degrees(25.0).refractdelta(),
            Angle::from_degminsec(0, 2, 9.2)
        );
        assert_eq!(
            Angle::from_degrees(-25.0).refract(),
            Angle::from_degminsec(-25, 0, 0.0)
        );
    }
}
