//! pracstro is an astronomy library made to be practical, principled, and easy to understand.
//! Its main use case is calculating the properties of celestial objects. Such as the Moon,
//! Sun, Planets, Stars, and other objects.
//!
//! It provides enough percision in its answers over a large enough range of time to be suitable
//! for most user-end applications, it is not perfect in its answers down to the smallest fractions
//! of arcseconds (nor could it be on a reasonable scale), but it does provide enough accuracy for
//! most use.
//!
//! This library is inspired off several old books, namely:
//!  * *Practical astronomy with your calculator* by Peter Duffett-Smith
//!  * *Astronomical Formulae for Calculators* by Jean Meeus

/// Time, Date, and Angle handling
pub mod time;

/// Coordinate handling
pub mod coord;

/// Solar Dynamics
pub mod sol;

/// Lunar Dynamics
pub mod moon;

/// Utility Functions
pub mod misc;
