use std::rc::Rc;

use crate::numbers::quad_field::*;
use crate::coord::*;

pub enum Pos {
    A,
    B,
    C
}

#[derive(Debug)]
pub struct MarkoffTriple<'a> {
    f: &'a QuadFieldExt<'a>,
    a: Rc<Coord<'a>>,
    b: Rc<Coord<'a>>,
    c: Rc<Coord<'a>>
}

impl<'a> MarkoffTriple<'a> {
    pub fn make(f: &'a Rc<QuadFieldExt<'a>>, a: u128, b: u128, c: u128) -> Self {
        MarkoffTriple {
            f,
            a: Rc::new(Coord::new(a, f)),
            b: Rc::new(Coord::new(b, f)),
            c: Rc::new(Coord::new(c, f))
        }
    }

    pub fn new(f: &'a Rc<QuadFieldExt<'a>>, a: &'a Rc<Coord>, b: &'a Rc<Coord>, c: &'a Rc<Coord>) -> Self {
        MarkoffTriple {
            f,
            a: Rc::clone(a),
            b: Rc::clone(b),
            c: Rc::clone(c)
        }
    }

    pub fn ord(&self, pos: Pos) -> u128 {
        match pos {
            Pos::A => &self.a,
            Pos::B => &self.b,
            Pos::C => &self.c
        }.get_ord(&self.f).value()
    }
}

