use either::*;
use std::rc::Rc;

use crate::numbers::quad_field::*;
use crate::numbers::fp::*;
use crate::numbers::factorization::*;

#[derive(Debug)]
pub struct Coord<const P: u128> {
    v: u128,
    chi: Either<QuadNum<P>, FpNum<P>>
}

impl<const P: u128> Coord<P> {
    pub fn new(v: u128, fp: &QuadField<P>) -> Coord<P> {
        let v3 = long_multiply(v, 3, P);
        let disc = intpow(v3, 2, P);
        let disc = (disc + P - 4) % P;
        let chi = match fp.int_sqrt_either(disc) {
            Left(mut x) => {
                x.a0 += v3;
                if x.a0 % 2 == 1 { x.a0 += P; }
                if x.a1 % 2 == 1 { x.a1 += P; }
                x.a0 /= 2;
                x.a1 /= 2;
                Left(x)
            },
            Right(mut x) => {
                x.value += v3;
                if x.value % 2 == 1 { x.value += P; } 
                x.value /= 2;
                Right(x)
            }
        };
        Coord {
            v: v % P,
            chi
        }
    }

    pub fn from_chi(chis: Either<(QuadNum<P>, QuadNum<P>), (FpNum<P>, FpNum<P>)>) -> Coord<P> {
        Coord {
            v: chis.clone().either(|mut l| {
                l.0.add(l.1);
                l.0.a0
            }, |mut r| {
                r.0.value = (r.0.value + r.1.value) % P;
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

    pub fn chi(&self) -> &Either<QuadNum<P>, FpNum<P>> {
        &self.chi
    }

    /*
    pub fn get_ord(&self, fp: &QuadField) -> Factorization {
        self.chi().as_ref().either(|l| l.order(fp), |r| r.order(fp.pminusone()))
    }
    */

    pub fn rot<'a>(self: &'a Rc<Self>, b: &'a Rc<Coord<P>>, c: &'a Rc<Coord<P>>, fp: &'a QuadField<P>) -> impl Iterator<Item = (Rc<Coord<P>>, Rc<Coord<P>>)> + 'a {
        std::iter::successors(Some((Rc::clone(b),Rc::clone(c))), move |(y,z)| {
            let tmp = long_multiply(3, self.v(), P);
            let tmp = long_multiply(tmp, z.v(), P);
            let (b_,c_) = (Rc::clone(z), Rc::new(Coord::new(tmp + P - y.v(), fp)));
            if &b_ == b && &c_ == c { None } else { Some((b_,c_)) }
        })
    }
}

impl<const P: u128> PartialEq for Coord<P> {
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}
impl<const P: u128> Eq for Coord<P> {}

fn either_multiply<const P: u128>(a: &mut Either<QuadNum<P>, FpNum<P>>, b: &Either<QuadNum<P>, FpNum<P>>, fp2: &QuadField<P>, fp: &FpStar<P>) {
    match (a, b) {
        (Left(x), Left(y)) => {
            x.multiply(y, fp2);
        },
        (Left(x), Right(y)) => {
            x.multiply(&QuadNum::from_ints(y.value, 0), fp2); 
        },
        (r @ Right(_), Left(y)) => {
            let v = r.clone().unwrap_right().value;
            let mut res = QuadNum::from_ints(v, 0);
            res.multiply(y, fp2);
            *r = Left(res);
        },
        (Right(x), Right(y)) => {
            x.multiply(y, fp);
        }
    }
}
