use either::*;
use std::rc::Rc;

use crate::quad_field::*;
use crate::fp::*;
use crate::factorization::*;

pub struct Coord {
    a: u128,
    chi: Option<Either<QuadNumber, FpNumber>>
}

impl Coord {
    pub fn new(a: u128) -> Coord {
        Coord {
            a,
            chi: None
        }
    }

    pub fn get_chi(&mut self, fp: &Rc<QuadFieldExt>) -> &Either<QuadNumber, FpNumber> {
        self.chi.get_or_insert_with(|| fp.int_sqrt_either(self.a))
    }

    pub fn get_ord(&mut self, fp: &Rc<QuadFieldExt>) -> Factorization {
        self.get_chi(fp).as_ref().either(|l| l.order(), |r| r.order())
    }

    pub fn rot<'a>(&'a self, b: u128, c: u128) -> impl Iterator<Item = (Coord, Coord)> + 'a {
        std::iter::successors(Some((b,c)), move |(y,z)| {
            let next = (*z, 3 * self.a * z - y);
            if next == (b, c) { None } else { Some(next) }
        })
        .map(|(y,z)| {
            (Coord::new(y), Coord::new(z))
        })
    }
}

fn either_multiply(a: &mut Either<QuadNumber, FpNumber>, b: &Either<QuadNumber, FpNumber>) {
    match (a, b) {
        (Left(x), Left(y)) => {
            x.multiply(y);
        },
        (Left(x), Right(y)) => {
            x.multiply(&x.group().from_ints(y.value(), 0)); 
        },
        (r @ Right(_), Left(y)) => {
            let v = r.clone().unwrap_right().value();
            let mut res = y.group().from_ints(v, 0);
            res.multiply(y);
            *r = Left(res);
        },
        (Right(x), Right(y)) => {
            x.multiply(y);
        }
    }
}
