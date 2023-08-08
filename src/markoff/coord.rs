use either::*;
use std::rc::Rc;

use crate::numbers::quad_field::*;
use crate::numbers::fp::*;
use crate::numbers::factorization::*;

#[derive(Debug)]
pub struct Coord {
    v: u128,
    chi: Either<QuadNum, FpNum>
}

impl Coord {
    pub fn new(v: u128, fp: &QuadField) -> Coord {
        let v3 = long_multiply(v, 3, fp.p());
        let disc = intpow(v3, 2, fp.p());
        let disc = (disc + fp.p() - 4) % fp.p();
        let chi = match fp.int_sqrt_either(disc) {
            Left(mut x) => {
                x.a0 += v3;
                if x.a0 % 2 == 1 { x.a0 += fp.p(); }
                if x.a1 % 2 == 1 { x.a1 += fp.p(); }
                x.a0 /= 2;
                x.a1 /= 2;
                Left(x)
            },
            Right(mut x) => {
                x.value += v3;
                if x.value % 2 == 1 { x.value += fp.p(); } 
                x.value /= 2;
                Right(x)
            }
        };
        Coord {
            v: v % fp.p(),
            chi
        }
    }

    pub fn from_chi(chis: Either<(QuadNum, QuadNum), (FpNum, FpNum)>, fp: &QuadField) -> Coord {
        Coord {
            v: chis.clone().either(|mut l| {
                l.0.add(l.1, fp);
                l.0.a0
            }, |mut r| {
                r.0.value = (r.0.value + r.1.value) % fp.p();
                r.0.value
            }),
            chi: chis
                .map_left(|x| x.0)
                .map_right(|x| x.0)
        }
    }

    pub fn v(&self) -> u128 {
        self.v 
    }

    pub fn chi(&self) -> &Either<QuadNum, FpNum> {
        &self.chi
    }

    pub fn get_ord(&self, fp: &QuadField) -> Factorization {
        self.chi().as_ref().either(|l| l.order(fp), |r| r.order(fp.pminusone()))
    }

    pub fn rot<'a>(self: &'a Rc<Self>, b: &'a Rc<Coord>, c: &'a Rc<Coord>, fp: &'a QuadField) -> impl Iterator<Item = (Rc<Coord>, Rc<Coord>)> + 'a {
        std::iter::successors(Some((Rc::clone(b),Rc::clone(c))), move |(y,z)| {
            let tmp = long_multiply(3, self.v(), fp.p());
            let tmp = long_multiply(tmp, z.v(), fp.p());
            let (b_,c_) = (Rc::clone(z), Rc::new(Coord::new(tmp + fp.p() - y.v(), fp)));
            if &b_ == b && &c_ == c { None } else { Some((b_,c_)) }
        })
    }
}

impl PartialEq for Coord {
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}
impl Eq for Coord {}

fn either_multiply(a: &mut Either<QuadNum, FpNum>, b: &Either<QuadNum, FpNum>, fp2: &QuadField, fp: &FpStar) {
    match (a, b) {
        (Left(x), Left(y)) => {
            x.multiply(y, fp2);
        },
        (Left(x), Right(y)) => {
            x.multiply(&QuadNum::from_ints(y.value(), 0), fp2); 
        },
        (r @ Right(_), Left(y)) => {
            let v = r.clone().unwrap_right().value();
            let mut res = QuadNum::from_ints(v, 0);
            res.multiply(y, fp2);
            *r = Left(res);
        },
        (Right(x), Right(y)) => {
            x.multiply(y, fp);
        }
    }
}

