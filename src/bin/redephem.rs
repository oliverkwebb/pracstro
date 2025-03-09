use std::io::{BufRead};
use std::fmt;
use pracstro::time;

#[derive(Debug, Copy, Clone)]
enum Property {}

#[derive(Debug, Copy, Clone)]
enum Value {
    Number(f64),
    Period(time::Period),
    ParamRef(Property),

    // Explicit Units
    Radians(f64),
    Degrees(f64),
    Clock(u8, u8, f64),
    Julian(f64),
    Calendar(u8, u8, f64),
}
impl fmt::Display for Value {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Value::*;
        match self {
            Number(x) => write!(f, "{}", x),

            Radians(x) => write!(f, "{} rad", x),
            Degrees(x) => write!(f, "{}°", x),
            Clock(x, y, z) => write!(f, "{}:{}:{}", x, y, z),
            _ => write!(f, "mrrow"),
        }
    }
}

fn reduce(a: Value, b: Value) -> Value {
    use Value::*;
    match a {
        Period(x) => match b {
            _ => a,
        },
        _ => a,
    }
}

fn resolve(query: String) -> Value {
    Value::Period(time::Period::from_radians(3.0))
}

fn main() {
    let stdin = std::io::stdin().lock().lines();
    stdin.for_each(|x| println!("{}", resolve(x.unwrap())));
}
