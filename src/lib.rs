#![warn(missing_docs)]
//! # Introduction
//! pracstro is an astronomy library that is practical, principled, and easy to understand.
//! It can be used for calculating properties of celestial objects, like the moon, sun,
//! planets, and stars. It's accurate enough across a large enough range of time for most user-end
//! applications. Its not perfect in its answers down to the arcsecond (nor could it be on a reasonable scale),
//! but is does provide enough accuracy for widgets, telescope pointing, scheduling ("Will I be able to see Venus this weekend?"), etc.
//!
//! The algorithms in this library are from several sources, mainly:
//!  * *Practical astronomy with Your Calculator* by Peter Duffett-Smith
//!  * *Astronomical Formulae for Calculators* by Jean Meeus
//!  * <https://www.celestialprogramming.com/> by Greg Miller
//!
//! Any substantial algorithm will be annotated with its source in this documentation.
//! This library started as an attempt to implement the algorithms from *Practical astronomy with Your Calculator*,
//! but branched out into different techniques.
//!
//! ```rust
//! use pracstro::*;
//! let now_date = time::Date::from_calendar(2025, 4, 16.0);
//! let now_time = time::Period::from_clock(19, 41, 11.0);
//! let my_latitude = time::Period::from_degrees(30.5);
//! let my_longitude = time::Period::from_degrees(-110.0);
//! sol::VENUS.location(now_date).horizon(now_date, now_time, my_latitude, my_longitude); // Get the horizonal coordinates of Venus
//! moon::MOON.phase(now_date).0; // The illuminated fraction of the moons surface
//! time::Period::from_degrees(120.0).clock(); // 16h00m00s
//! ```
//!
//! # Structure
//! This library consists of 4 main modules, which build upon the ones before them:
//! 1. [`time`], which contains functions for conversion and handling of times and dates to and from different formats
//! 2. [`coord`], which contains functions for handling coordinates (pairs of angles) and transforming them into several desired representations.
//! 3. [`sol`], which contains algorithms that can calculate the position of planets to a accurate precision.
//! 4. [`moon`], which contains algorithms for getting the position and phase of the moon.
//!
//! Each of these modules have one or two core types that represent a certain kind of data::
//! - [`Date`](time::Date) - An instant in continuous time.
//! - [`Period`](time::Period) - An angle automatically corrected to be between \[0°, 360°\]. Which can also represent a time of day.
//! - [`Coord`](coord::Coord) - A pair of angles, representing latitude/longitude on a sphere.
//! - [`Planet`](sol::Planet) - A planet and its orbital properties, along with properties required for orbital correction.
//! - [`Moon`](moon::Moon) - The moons orbital properties at a epoch.
//!
//! These types have methods to get the properties of this data. Primarily in pairs of methods that convert to/from a certain
//! representation of that data. Although lone methods that get certain data for a type do exist.
//! ```rust
//! use pracstro::*;
//! time::Period::from_radians(time::Period::from_decimal(16.0).radians()).clock();
//! ```

pub mod time;

pub mod coord;

pub mod sol;

/// Lunar Dynamics
pub mod moon;

/// Utility Functions
pub mod misc;
