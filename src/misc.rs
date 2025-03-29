//! Miscellaneous non-dependent algorithms

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
    fn test_easter() {
        assert_eq!(easter(2000), (4, 23));
        assert_eq!(easter(2024), (3, 31));
    }
}
