/// Calculate the date of easter
pub fn easter(year: i32) -> (i32, i32) {
    let a = year % 19;
    let b = year / 100;
    let c = year % 100;
    let d = b / 4;
    let e = b % 4;
    let f = (b + 8) / 25;
    let g = (b - f + 1) / 3;
    let h = ((19 * a) + b - d - g + 15) % 30;
    let i = c / 4;
    let k = c % 4;
    let l = (32 + (2 * e) + (2 * i) - h - k) % 7;
    let m = (a + (11 * h) + (22 * l)) / 451;
    let n = (h + l - (7 * m) + 114) / 31;
    let p = (h + l - (7 * m) + 114) % 31;
    (n, p + 1)
}

#[derive(Debug, PartialEq)]
pub struct YearMonthDate {
    year: f64,
    month: f64,
    day: f64,
}

/// Days since January 0th 4713 BC
#[derive(Debug, PartialEq)]
pub struct JulDate(f64);

#[derive(Debug, PartialEq)]
pub struct ClockTime {
    hour: u8,
    minute: u8,
    second: f64,
}

#[derive(Debug, PartialEq)]
pub struct DecimalHrs(f64);

/// True Modulus Operation
pub fn lpr(x: f64, y: f64) -> f64 {
    let z = x % y;
    match z < 0.0 {
        true => z + y,
        false => z,
    }
}

/// Julian Day -> Day/Month/Year
impl From<JulDate> for YearMonthDate {
    fn from(jday: JulDate) -> Self {
        let j = jday.0 + 0.5;
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

        YearMonthDate {
            day: c - e + f - (30.6001 * g).trunc(),
            month: m,
            year: y,
        }
    }
}

/// Day, Month, Year -> Julian Day
impl From<YearMonthDate> for JulDate {
    fn from(t: YearMonthDate) -> Self {
        let (day, mut month, mut year) = (t.day, t.month, t.year);
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

        JulDate(b + c + (30.6001 * (mp + 1.0)).trunc() + day + 1_720_994.5)
    }
}

impl From<ClockTime> for DecimalHrs {
    fn from(t: ClockTime) -> Self {
        DecimalHrs((t.hour as f64) + (((t.minute as f64) + (t.second / 3600.0)) / 60.0))
    }
}

impl From<DecimalHrs> for ClockTime {
    fn from(t: DecimalHrs) -> Self {
        ClockTime {
            hour: t.0.trunc() as u8,
            minute: ((t.0.fract()) * 60.0).trunc() as u8,
            second: ((t.0.fract()) * 60.0).fract() * 60.0,
        }
    }
}

impl DecimalHrs {
    pub fn correctz(self, tzoff: i32) -> Self {
        DecimalHrs(lpr(self.0 + tzoff as f64, 24.0))
    }
    pub fn ucorrectz(self, tzoff: i32) -> Self {
        DecimalHrs(lpr(self.0 - tzoff as f64, 24.0))
    }
}

pub fn gst_to_dec(decimal: DecimalHrs, jd: JulDate) -> DecimalHrs {
    let jday = jd.0;
    let s = jday - 2451545.0;
    let t = s / 36525.0;
    let t0 = lpr(
        6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
        24.0,
    );
    DecimalHrs(lpr(decimal.0 - t0, 24.0) * 0.9972695663)
}

pub fn jdec_to_gst(decimal: DecimalHrs, jd: JulDate) -> DecimalHrs {
    let jday = jd.0;
    let s = jday - 2451545.0;
    let t = s / 36525.0;
    let t0 = lpr(
        6.697374558 + (2400.051336 * t) + (0.000025862 * (t * t)),
        24.0,
    );
    DecimalHrs(lpr(t0 + (decimal.0 * 1.002737909), 24.0))
}

pub fn gst_to_lst(gst: DecimalHrs, longit: f64) -> DecimalHrs {
    DecimalHrs(lpr(gst.0 + (longit / 15.0), 24.0))
}

pub fn lst_to_gst(lst: DecimalHrs, longit: f64) -> DecimalHrs {
    DecimalHrs(lpr(lst.0 - (longit / 15.0), 24.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_easter() {
        assert_eq!(easter(2000), (4, 23));
        assert_eq!(easter(2024), (3, 31));
    }

    #[test]
    fn test_lpr() {
        assert_eq!(lpr(-1.0, 5.0), 4.0);
    }

    #[test]
    fn test_julian() {
        assert_eq!(
            YearMonthDate::from(JulDate(2_446_113.75)),
            YearMonthDate {
                day: 17.25,
                month: 2.0,
                year: 1985.0
            }
        );
        assert_eq!(
            JulDate::from(YearMonthDate {
                day: 17.25,
                month: 2.0,
                year: 1985.0
            }),
            JulDate(2_446_113.75)
        );
    }

    #[test]
    fn test_decimalhrs() {
        assert_eq!(
            DecimalHrs::from(ClockTime {
                hour: 18,
                minute: 31,
                second: 27.0
            }),
            DecimalHrs(18.516791666666666)
        );
        assert_eq!(
            ClockTime::from(DecimalHrs(11.75)),
            ClockTime {
                hour: 11,
                minute: 45,
                second: 0.0
            }
        );
    }

    #[test]
    fn test_gst() {
        assert_eq!(
            jdec_to_gst(DecimalHrs(14.614_353), JulDate(2_444_351.5)),
            DecimalHrs(4.668119549708194)
        );
        assert_eq!(
            gst_to_dec(DecimalHrs(4.668119549708194), JulDate(2_444_351.5)),
            DecimalHrs(14.614352994461141)
        );
    }

    #[test]
    fn test_lst() {
        assert_eq!(
            gst_to_lst(DecimalHrs(4.668_119), -64.0),
            DecimalHrs(0.4014523333333333)
        );
        assert_eq!(
            lst_to_gst(DecimalHrs(0.401_453), -64.0),
            DecimalHrs(4.668119666666667)
        );
    }
}
