use pracstro::{sol, time};

fn main() {
    let now = time::Date::now();
    for p in sol::PLANETS {
        let (ra, de) = p.location(now).celestial();
        let ((rah, ram, _), (ded, dem, _)) = (ra.clock(), de.to_latitude().degminsec());
        println!(
            "{:<10} {:>2}h{:02} RA {:>3}°{:02}' De {:.2} AU",
            p.name,
            rah,
            ram,
            ded,
            dem,
            p.distance(now)
        );
    }
}
