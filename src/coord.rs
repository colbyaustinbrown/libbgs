use either::*;
use std::rc::Rc;

use crate::numbers::quad_field::*;
use crate::numbers::fp::*;
use crate::numbers::factorization::*;

#[derive(Debug)]
pub struct Coord<'a> {
    v: u128,
    chi: Either<QuadNumExt<'a>, FpNum>
}

impl<'a> Coord<'a> {
    pub fn new(v: u128, fp: &'a QuadFieldExt) -> Coord<'a> {
        let disc = intpow(v, 2, fp.p());
        let disc = (disc * disc + fp.p() - 4) % fp.p();
        let chi = match fp.int_sqrt_either(disc) {
            Left(mut x) => {
                x.a0 += v;
                if x.a0 % 2 == 1 { x.a0 += fp.p(); }
                if x.a1 % 2 == 1 { x.a1 += fp.p(); }
                x.a0 /= 2;
                x.a1 /= 2;
                Left(x)
            },
            Right(mut x) => {
                x.value += v;
                if x.value % 2 == 1 { x.value += fp.p(); } 
                x.value /= 2;
                Right(x)
            }
        };
        Coord {
            v,
            chi
        }
    }

    pub fn v(&self) -> u128 {
        self.v 
    }

    pub fn chi(&self) -> &Either<QuadNumExt, FpNum> {
        &self.chi
    }

    pub fn get_ord(&self, fp: &QuadFieldExt) -> Factorization {
        self.chi().as_ref().either(|l| l.order(fp), |r| r.order(fp.pminusone()))
    }

    pub fn rot(self: &'a Rc<Self>, b: &'a Rc<Coord<'a>>, c: &'a Rc<Coord<'a>>, fp: &'a QuadFieldExt) -> impl Iterator<Item = (Rc<Coord<'a>>, Rc<Coord<'a>>)> + 'a {
        std::iter::successors(Some((Rc::clone(b),Rc::clone(c))), move |(y,z)| {
            let (b_,c_) = (Rc::clone(z), Rc::new(Coord::new(3 * self.v() * z.v() - y.v(), fp)));
            if Rc::as_ptr(&b_) == Rc::as_ptr(b) { None } else { Some((b_,c_)) }
        })
    }
}

impl<'a> PartialEq for Coord<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}
impl<'a> Eq for Coord<'a> {}

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

