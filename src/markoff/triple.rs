use std::rc::Rc;

use crate::markoff::Coord;
use crate::numbers::*;
use crate::util::*;

pub enum Pos {
    A,
    B,
    C,
}

#[derive(Debug)]
pub struct MarkoffTriple<const P: u128> {
    pub a: Rc<Coord<P>>,
    pub b: Rc<Coord<P>>,
    pub c: Rc<Coord<P>>,
}

impl<const P: u128> MarkoffTriple<P> {
    pub fn make(a: u128, b: u128, c: u128, f: &QuadField<P>) -> Self {
        MarkoffTriple {
            a: Rc::new(Coord::new(a, f)),
            b: Rc::new(Coord::new(b, f)),
            c: Rc::new(Coord::new(c, f)),
        }
    }

    pub fn new(a: &Rc<Coord<P>>, b: &Rc<Coord<P>>, c: &Rc<Coord<P>>) -> Self {
        MarkoffTriple {
            a: Rc::clone(a),
            b: Rc::clone(b),
            c: Rc::clone(c),
        }
    }

    pub fn ord(
        &self,
        pos: Pos,
        f: &QuadField<P>,
        minusonesize: &Factorization,
        plusonesize: &Factorization,
    ) -> u128 {
        match pos {
            Pos::A => &self.a,
            Pos::B => &self.b,
            Pos::C => &self.c,
        }
        .get_ord(f, minusonesize, plusonesize)
        .value()
    }

    pub fn rot<'a>(
        &'a self,
        pos: Pos,
        f: &'a QuadField<P>,
    ) -> impl Iterator<Item = MarkoffTriple<P>> + 'a {
        match pos {
            Pos::A => self.a.rot(&self.b, &self.c, f),
            Pos::B => self.b.rot(&self.a, &self.c, f),
            Pos::C => self.c.rot(&self.a, &self.b, f),
        }
        .map(move |(x, y)| match pos {
            Pos::A => MarkoffTriple {
                a: Rc::clone(&self.a),
                b: x,
                c: y,
            },
            Pos::B => MarkoffTriple {
                a: x,
                b: Rc::clone(&self.b),
                c: y,
            },
            Pos::C => MarkoffTriple {
                a: x,
                b: y,
                c: Rc::clone(&self.c),
            },
        })
    }

    pub fn get_from_ab(x: &Rc<Coord<P>>, y: &Rc<Coord<P>>, f: &QuadField<P>) -> Vec<Self> {
        let a = x.v();
        let b = y.v();
        let a2 = intpow(a, 2, P);
        let b2 = intpow(b, 2, P);
        let mut disc = FpNum::<P>::from(long_multiply(long_multiply(9, a2, P), b2, P));
        disc += P;
        disc -= long_multiply(4, a2 + b2, P);

        let mut res = long_multiply(3, long_multiply(a, b, P), P);
        if res % 2 == 1 {
            res += P
        };
        res /= 2;
        match disc.int_sqrt().map(|z| z.into()) {
            None => Vec::new(),
            Some(0) => vec![res],
            Some(mut d) => {
                if d % 2 == 1 {
                    d += P;
                }
                d /= 2;
                vec![(res + d) % P, (res + P - d) % P]
            }
        }
        .iter()
        .map(|c| MarkoffTriple::new(x, y, &Rc::new(Coord::new(*c, f))))
        .collect()
    }
}
