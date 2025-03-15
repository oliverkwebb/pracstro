#![warn(missing_docs)]
/*!
# pracstro - Compact Astronomy Library and Ephemeris Generator

pracstro is an astronomy library made from a collection of other algorithms that is
compact, principled, and easy to understand. It's made for calculating properties
of celestial objects, such as the moon, sun, planets, and stars.

```
use pracstro::*;

let now_date = time::Date::from_calendar(2025, 4, 16.0);
let now_time = time::Period::from_clock(19, 41, 11.0);
let my_latitude = time::Period::from_degrees(30.5);
let my_longitude = time::Period::from_degrees(-110.0);

sol::VENUS.location(now_date).horizon(now_date, now_time, my_latitude, my_longitude); // Get the horizontal coordinates of Venus
moon::MOON.phase(now_date).0; // The illuminated fraction of the moons surface
time::Period::from_degrees(120.0).clock(); // 16h00m00s
```

# Benchmarks

Although speed is not a direct goal of this library, the simplicity of the algorithms used often makes
the library much faster than other astronomy libraries. The tradeoff for this simplicity is
some accuracy, though the library remains accurate enough for most real-world use.

Average of many tests, n = 40,000:

| Test           | `pracstro`           | [`astro`](https://crates.io/crates/astro) |
|----------------|----------------------|---------------------|
| Moon Phase     | 558ns                | 2,979ns (3µs)       |
| Jupiter Coords | 601ns                | 70,860ns (70µs)     |
| Full ephemeris | 3,406ns (3.4µs)      | 1,208,833ns (1.2ms) |

# Structure
This library contains 4 primary modules, which build upon the ones before them:
1. [`time`] for the conversion and representation of times, dates, and angles.
2. [`coord`] for the conversion and representation of coordinates.
3. [`sol`] for the calculation of properties of planets and the sun.
4. [`moon`] for the calculation of properties of the moon.

Each of these have one or two types that represent a certain kind of data:
- [`Date`](time::Date) - An instant in continuous time.
- [`Period`](time::Period) - An angle automatically corrected to be between \[0°, 360°\]. Which can also represent a time of day.
- [`Coord`](coord::Coord) - A pair of angles, representing latitude/longitude on a sphere.
- [`Planet`](sol::Planet) - A planets orbital properties, along with data required for orbital correction.
- [`Moon`](moon::Moon) - The moons orbital properties.

These types have methods to get the properties of this data. Primarily in pairs of methods that convert to/from a certain
representation of that data. Although lone methods that get certain data for a type do exist.
*/

pub mod time;

pub mod coord;

pub mod sol;

/// Lunar Dynamics
pub mod moon;

/// Utility Functions
pub mod misc;
