use crate::{coord, time};

/// Calculate the coordinates of the sun at a given time
///
/// From https://www.celestialprogramming.com/snippets/sunPositionVsop.html
pub fn whereis_sun(d: time::Date) -> coord::Coord {
    let t = (d.julian() - 2451545.0) / 365250.0;
    let mut y = 0.00010466965 * (0.09641690558 + 18849.22754997420 * t).cos();
    y = y + 0.00835292314 * (0.13952878991 + 12566.15169998280 * t).cos() - 0.02442699036;
    y = y + 0.99989211030 * (0.18265890456 + 6283.07584999140 * t).cos();
    y = y + (0.00093046324 + 0.00051506609 * (4.43180499286 + 12566.15169998280 * t).cos()) * t;

    let mut x = 0.00561144206 + 0.00010466628 * (1.66722645223 + 18849.22754997420 * t).cos();
    x = x + 0.00835257300 * (1.71034539450 + 12566.15169998280 * t).cos();
    x = x + 0.99982928844 * (1.75348568475 + 6283.07584999140 * t).cos();
    x = x + (0.00123403056 + 0.00051500156 * (6.00266267204 + 12566.15169998280 * t).cos()) * t;

    let z = 0.00227822442 * (3.41372504278 + 6283.07584999140 * t).cos() * t;

    let tx = -(x + y * 0.000000440360 + z * -0.000000190919);
    let ty = -(x * -0.000000479966 + y * 0.917482137087 + z * -0.397776982902);
    let tz = -(y * 0.397776982902 + z * 0.917482137087);

    let r = (tx * tx + ty * ty + tz * tz).sqrt();
    let l = time::Period::from_radians(ty.atan2(tx));
    let t2 = time::Period::from_radians(0.5 * std::f64::consts::PI - (tz / r).acos());

    coord::Coord::from_celestial(l, t2)
}

/// Calculates the sun rise and set time
///
/// Calculating the second-level time for any altitude and weather is essentially impossible.
/// Your altitude and temperature will effect the result far more than angular size and parallax.
pub fn sun_rise_set(d: time::Date, lati: time::Period, longi: time::Period) -> (time::Period, time::Period) {
	whereis_sun(d).riseset(d, lati, longi).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sunpos() {
        assert_eq!(
            whereis_sun(time::Date::from_julian(2268932.541667)),
            coord::Coord::from_celestial(
                time::Period::from_degrees(298.49306),
                time::Period::from_degrees(-20.91664)
            )
        );
    }

    // "Is this a reliable way of getting the ecliptic longitude of the sun?"
    #[test]
    fn test_lambdasun() {
    	assert_eq!(
            whereis_sun(time::Date::from_calendar(1980, 7, 27.0)).ecliptic(time::Date::from_calendar(1980, 7, 27.0)).0,
            time::Period::from_degminsec(124, 23, 40.8)
    	)
    }

}
