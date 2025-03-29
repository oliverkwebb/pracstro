use pracstro::{moon, sol, time};
use std::time::Instant;

fn run_test(name: &str, n: u32, run: fn() -> ()) {
    let start = Instant::now();
    (1..n).for_each(|_| run());
    let time = start.elapsed();
    println!("{},{:?},{:?}", name, time / n, time);
}

fn main() {
    let n = 40000;
    println!("Test,Mean (n={}),Full", n);

    run_test("Control (NOP)", n, || ());
    run_test("Current time", n, || {
        time::Date::now();
    });
    fn ephem() {
        let now = time::Date::now();
        for p in sol::PLANETS {
            let (ra, de) = p.location(now).equatorial();
            (ra.clock(), de.to_latitude().degminsec());
        }
    }
    run_test("Full ephemeris", n, ephem);
    run_test("Moon Phase", n, || {
        moon::MOON.illumfrac(time::Date::now());
    });
}
