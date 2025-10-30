//! Celestial object trait for generics

use crate::coord::Coord;
use crate::time;

/// A celestial object in pracstro is defined by the ability to query its cartesian coordinates from time
pub trait CelObj {
    /// The cartesian coordinates of the object
    fn locationcart(&self, d: time::Date) -> (f64, f64, f64);

    /// The 2D Polar Coordinates of the object
    fn location(&self, d: time::Date) -> Coord {
        let (x, y, z) = self.locationcart(d);
        Coord::from_cartesian(x, y, z)
    }

    /// The distance from the reference frame to the object, in AU
    fn distance(&self, d: time::Date) -> f64 {
        let (x, y, z) = self.locationcart(d);
        (x * x + y * y + z * z).sqrt()
    }
}
