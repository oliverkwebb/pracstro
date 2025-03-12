/// Miscellaneous non-dependent algorithms

/// Solve the equation of Kepler
pub fn kepler(mut m: f64, ecc: f64) -> f64 {
    m = m.to_radians();
    let mut e = m;
    let mut delta: f64;

    delta = e - ecc * e.sin() - m;
    e -= delta / (1.0 - ecc * e.cos());
    while delta.abs() > 1E-6 {
        delta = e - ecc * e.sin() - m;
        e -= delta / (1.0 - ecc * e.cos());
    }
    return e;
}

/// Calculate the date of Easter
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_easter() {
        assert_eq!(easter(2000), (4, 23));
        assert_eq!(easter(2024), (3, 31));
    }
}
