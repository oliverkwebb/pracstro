use pracstro::{moon,time};

const HALFMON: f64 = 14.76529434;
const EMOJIS: [&str; 8] = ["🌑", "🌒", "🌓", "🌔", "🌕", "🌖", "🌗", "🌘"];
const SEMOJI: [&str; 8] = ["🌑", "🌘", "🌗", "🌖", "🌕", "🌔", "🌓", "🌒"];
const PNAMES: [&str; 8] = [
    "New",
    "Waxing Crescent",
    "First Quarter",
    "Waxing Gibbous",
    "Full",
    "Waning Gibbous",
    "Last Quarter",
    "Waning Crescent",
];

fn phaseidx(ilumfrac: f64, mage: f64) -> usize {
    let half: bool = mage > HALFMON;
    match (ilumfrac, half) {
        (0.00..0.04, _) => 0,
        (0.96..1.00, _) => 4,
        (0.46..0.54, true) => 6,
        (0.46..0.54, false) => 2,
        (0.54..0.96, true) => 5,
        (0.54..0.96, false) => 3,
        (_, true) => 7,
        (_, false) => 1,
    }
}

fn main() {
	let p = moon::MOON.phase(time::Date::now());
    println!("{} {} ({:.2}%)", PNAMES[phaseidx(p.0, p.1)], EMOJIS[phaseidx(p.0, p.1)], p.0 * 100.0);
}
