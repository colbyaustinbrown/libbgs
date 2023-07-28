use std::rc::Rc;

use crate::quad_field::*;
use crate::coord::*;

enum Pos {
    A,
    B,
    C
}

struct MarkoffTriple<'a> {
    fp: Rc<QuadFieldExt<'a>>,
    a: Coord<'a>,
    b: Coord<'a>,
    c: Coord<'a>
}

impl<'a> MarkoffTriple<'a> {
    pub fn ord(&'a mut self, pos: Pos) -> u128 {
        match pos {
            Pos::A => &mut self.a,
            Pos::B => &mut self.b,
            Pos::C => &mut self.c
        }.get_ord(&self.fp).value()
    }
}

