use std::rc::Rc;

use crate::quad_field::*;
use crate::coord::*;

enum Pos {
    A,
    B,
    C
}

struct MarkoffTriple {
    fp: Rc<QuadFieldExt>,
    a: Coord,
    b: Coord,
    c: Coord
}

impl MarkoffTriple {
    pub fn ord(&mut self, pos: Pos) -> u128 {
        match pos {
            Pos::A => &mut self.a,
            Pos::B => &mut self.b,
            Pos::C => &mut self.c
        }.get_ord(&self.fp).value()
    }
}

