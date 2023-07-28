use either::*;

use crate::quad_field::*;
use crate::fp::*;
use crate::factorization::*;

pub struct Coord<'a> {
    a: u128,
    chi: Option<Either<QuadNumExt<'a>, FpNum>>
}

impl<'a> Coord<'a> {
    pub fn new(a: u128) -> Coord<'a> {
        Coord {
            a,
            chi: None
        }
    }

    pub fn get_chi(&mut self, fp: &'a QuadFieldExt) -> &Either<QuadNumExt, FpNum> {
        self.chi.get_or_insert_with(|| fp.int_sqrt_either(self.a))
    }

    pub fn get_ord(&mut self, fp: &'a QuadFieldExt) -> Factorization {
        self.get_chi(fp).as_ref().either(|l| l.order(fp), |r| r.order(fp.pminusone()))
    }

    pub fn rot(&self, b: u128, c: u128) -> impl Iterator<Item = (Coord, Coord)> {
        std::iter::successors(Some((b,c)), move |(y,z)| {
            let next = (*z, 3 * self.a * z - y);
            if next == (b, c) { None } else { Some(next) }
        })
        .map(|(y,z)| {
            (Coord::new(y), Coord::new(z))
        })
    }
}

fn either_multiply<F>(a: &mut Either<QuadNum<F>, FpNum>, b: &Either<QuadNum<F>, FpNum>, fp2: &F, fp: &FpStar) 
where F: QuadField<Elem = QuadNum<F>> {
    match (a, b) {
        (Left(x), Left(y)) => {
            x.multiply(y, fp2);
        },
        (Left(x), Right(y)) => {
            x.multiply(&fp2.from_ints(y.value(), 0), fp2); 
        },
        (r @ Right(_), Left(y)) => {
            let v = r.clone().unwrap_right().value();
            let mut res = fp2.from_ints(v, 0);
            res.multiply(y, fp2);
            *r = Left(res);
        },
        (Right(x), Right(y)) => {
            x.multiply(y, fp);
        }
    }
}

