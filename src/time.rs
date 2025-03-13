//!
//! This module contains functions for the handling and conversion of Times, Dates, and Angles.
//! This data can be represented in two types
//! - The `Date` type, which represents an instance in continuous time
//! - The `Period` type, which represents anything which modulo arithmetic

use std::fmt;

/// True Modulus Operation
fn lpr(x: f64, y: f64) -> f64 {
    let z = x % y;
    match z < 0.0 {
        true => z + y,
        false => z,
    }
}

/// Continious Instant in time
/// Julian Day, Epoch is Jan 0 4713 BC
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Date(f64);
/// Properties of concern:
///     Calendar (Year, Month, Day, etc)
///     Julian Date
impl Date {
    pub fn julian(self) -> f64 {
        self.0
    }
    pub fn from_julian(x: f64) -> Self {
        Date(x)
    }

    /// Year, Month, Day (time is Period::from_decimal(day.fract()))
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

    pub fn unix(self) -> f64 {
        (self.julian() - 2440587.5) * 86400.0
    }
    pub fn from_unix(t: f64) -> Self {
        Date::from_julian((t / 86400.0) + 2440587.5)
    }

    /// Sunday is 0
    pub fn dow(self) -> u8 {
        ((self.julian() + 1.5) / 7.0).fract().round() as u8
    }
}

/// Angles and Time being the most prominent use for this type
/// Radians, [0, Tau (2pi)]
#[derive(Clone, Copy)]
pub struct Period(f64);
/// Properties of concern:
///     Radians
///     Degrees
///     Decimal Hours (i.e. 11.5)
///     Hour, Minute, (Second[.Subsecond]) (i.e. 11:30)
impl Period {
    /// This is the only function that directly reads to Period
    pub fn radians(self) -> f64 {
        self.0
    }
    /// This is the only function that directly writes to Period
    /// ALWAYS DO THIS UNDER A (TRUE) MODULO GUARD WITH LPR
    pub fn from_radians(x: f64) -> Self {
        Period(lpr(x, std::f64::consts::TAU))
    }

    pub fn degrees(self) -> f64 {
        self.radians().to_degrees()
    }
    pub fn from_degrees(x: f64) -> Self {
        Period::from_radians(x.to_radians())
    }

    pub fn decimal(self) -> f64 {
        self.degrees() / 15.0 // 1 hour <-> 15 degrees, 360/24 = 15
    }
    pub fn from_decimal(x: f64) -> Self {
        Period::from_degrees(x * 15.0)
    }

    /// Returns (Hour, Minute, Second)
    pub fn clock(self) -> (u8, u8, f64) {
        let y = self.decimal();
        (
            y.trunc() as u8,
            (y.fract() * 60.0).trunc() as u8,
            (y.fract() * 60.0).fract() * 60.0,
        )
    }
    // Converts from clocktime
    pub fn from_clock(h: u8, m: u8, s: f64) -> Self {
        Period::from_decimal((h as f64) + (((m as f64) + (s / 60.0)) / 60.0))
    }

    pub fn degminsec(self) -> (i16, u8, f64) {
        let y = self.degrees();
        (
            y.trunc() as i16,
            (y.fract() * 60.0).trunc() as u8,
            (y.fract() * 60.0).fract() * 60.0,
        )
    }
    pub fn from_degminsec(d: i16, m: u8, s: f64) -> Self {
        Period::from_degrees((d as f64) + (((m as f64) + (s / 60.0)) / 60.0))
    }

    // Converts to siderial time
    pub fn gst(self, date: Date) -> Self {
        let jday = date.julian();
        let s = jday - 2451545.0;
        let t = s / 36525.0;
        let t0 = lpr(
            6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
            24.0,
        );
        Period::from_decimal(lpr(t0 + (self.decimal() * 1.002737909), 24.0))
    }
    // Converts from siderial time
    pub fn ungst(self, date: Date) -> Self {
        let jday = date.julian();
        let s = jday - 2451545.0;
        let t = s / 36525.0;
        let t0 = lpr(
            6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
            24.0,
        );
        Period::from_decimal(lpr(self.decimal() - t0, 24.0) * 0.9972695663)
    }

    pub fn fixquad(self) -> Self {
        let a = self.degrees();
        Period::from_degrees(a - 360.0 * (a / 360.0).floor())
    }

    /// Addition, For timezones and LST
    pub fn add(self, x: Self) -> Self {
        Period::from_radians(self.radians() + x.radians())
    }
    /// Subtraction, For timezones and LST
    pub fn sub(self, x: Self) -> Self {
        Period::from_radians(self.radians() - x.radians())
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
    /// Arcsine of Period
    pub fn asin(x: f64) -> Self {
        Period::from_radians(x.asin())
    }
    /// Arccosine of Period
    pub fn acos(x: f64) -> Self {
        Period::from_radians(x.acos())
    }
    /// Arctangent of Period
    pub fn atan(x: f64) -> Self {
        Period::from_radians(x.atan())
    }
    pub fn atan2(x: f64, y: f64) -> Self {
        Period::from_radians(x.atan2(y))
    }
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
impl PartialEq for Period {
    fn eq(&self, other: &Self) -> bool {
        let (d, m, _) = self.degminsec();
        let (d2, m2, _) = other.degminsec();
        d == d2 && m == m2
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
}
