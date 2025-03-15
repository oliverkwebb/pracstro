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

# Goals
* **Simplicity** - The algorithms in this library should be parsimonious enough to be transcribed into several different languages.
* **0-dependency** - People using this library should not have to pull in 400 crates, it should do its job on its own accord.
* **Understandable Design** - The API for this library should be modular and understandable.
* **Correctness** - Although calculating the exact position of a planet or even star is not feasible in a library designed to be simple, algorithms
provided by this library should be good enough most real-world use.
* **Testing and documentation** - The bulk of the source code of this library is tests and documentation. As a general rule, bugs are a bad thing.

# Non-Goals
* **Complete Accuracy** - There are an extremely large amount of factors which can throw values off over a long enough time. It is fortunate that most of these factors are negligible over time-spans of hundreds, or even in some cases thousands of years. Additionally, some of these effects (altitude, temperature, air-pressure) are not predictable.
* **Cataloging** - This library contains information about the Moon, Sun, and Planets (including Pluto). it does not however catalog minor planets or stars

# Precision

This library aims to be accurate enough across a large enough range of time for most user-end applications.
Its not perfect in its answers down to the arcsecond (nor could it be on a reasonable scale),
but is does provide enough accuracy for widgets, telescope pointing, scheduling ("Will I be able to see Venus this weekend?"), etc.

High precision astronomy libraries that work over vast spans of time contain impressive levels of code in

```
use pracstro::*;
time::Period::from_radians(time::Period::from_decimal(16.0).radians()).clock();
```

# Resources

Any substantial algorithm implemented from another source will be annotated with its source in this documentation.
This library started as an attempt to implement the algorithms from *Practical astronomy with Your Calculator*,
but branched out into different techniques. Most test cases in this library were either pulled from *Practical Astronomy with Your Calculator*
or created from data from Stellariun.

The algorithms in this library are from several sources, mainly:
 * <https://www.celestialprogramming.com/> by Greg Miller
 * *Practical astronomy with Your Calculator* by Peter Duffett-Smith
 * *Astronomical Formulae for Calculators* by Jean Meeus

Additional tools similar to this in job:
 * JPL's SPICE toolkit: https://naif.jpl.nasa.gov/naif/toolkit.html
 * JPL's Online Horizons Ephemeris Generator: https://ssd.jpl.nasa.gov/horizons/
 * Stellarium - Full GUI Planetarium: https://stellarium.org/
