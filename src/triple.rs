use std::rc::Rc;

use crate::numbers::quad_field::*;
use crate::coord::*;
use crate::util::*;

pub enum Pos {
    A,
    B,
    C
}

#[derive(Debug)]
pub struct MarkoffTriple {
    a: Rc<Coord>,
    b: Rc<Coord>,
    c: Rc<Coord>
}

impl MarkoffTriple {
    pub fn make(a: u128, b: u128, c: u128, f: &QuadField) -> Self {
        MarkoffTriple {
            a: Rc::new(Coord::new(a, f)),
            b: Rc::new(Coord::new(b, f)),
            c: Rc::new(Coord::new(c, f))
        }
    }

    pub fn new(a: &Rc<Coord>, b: &Rc<Coord>, c: &Rc<Coord>) -> Self {
        MarkoffTriple {
            a: Rc::clone(a),
            b: Rc::clone(b),
            c: Rc::clone(c)
        }
    }

    pub fn ord(&self, pos: Pos, f: &QuadField) -> u128 {
        match pos {
            Pos::A => &self.a,
            Pos::B => &self.b,
            Pos::C => &self.c
        }.get_ord(f).value()
    }

    pub fn get_from_ab(x: &Rc<Coord>, y: &Rc<Coord>, f: &QuadField) -> Vec<Self> {
        let a = x.v();
        let b = y.v();
        let m = f.p();
        let a2 = intpow(a, 2, m);
        let b2 = intpow(b, 2, m);
        let mut disc = f.pminusone().from_int(long_multiply(long_multiply(9, a2, m), b2, m));
        disc.value += m;
        disc.value -= long_multiply(4, a2 + b2, m);
        disc.value %= m;

        let mut res = long_multiply(3, long_multiply(a, b, m), m);
        if res % 2 == 1 { res += m };
        res /= 2;
        match disc.int_sqrt(f.pminusone()).map(|z| z.value) {
            None => Vec::new(),
            Some(0) => vec![res],
            Some(mut d) => {
                if d % 2 == 1 { d += m; } 
                d /= 2;
                vec![(res + d) % m, (res + m - d) % m]
            }
        }.iter().map(|c| 
            MarkoffTriple::new(x, y, &Rc::new(Coord::new(*c, f)))
        ).collect()
    }
}

