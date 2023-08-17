use either::*;
use std::rc::Rc;

use crate::numbers::factorization::*;
use crate::numbers::fp::*;
use crate::numbers::quad_field::*;
use crate::numbers::sylow::*;

#[derive(Debug)]
pub struct Coord<const P: u128> {
    v: u128,
    chi: Either<QuadNum<P>, FpNum<P>>,
}

impl<const P: u128> Coord<P> {
    pub fn new(v: u128, fp: &QuadField<P>) -> Coord<P> {
        let v3 = long_multiply(v, 3, P);
        let disc = intpow(v3, 2, P);
        let disc = (disc + P - 4) % P;
        let chi = fp.int_sqrt_either(disc).map_either(
            |mut x| {
                x.0 += v3;
                if x.0 % 2 == 1 {
                    x.0 += P;
                }
                if x.1 % 2 == 1 {
                    x.1 += P;
                }
                x.0 /= 2;
                x.1 /= 2;
                x
            },
            |mut x| {
                x.0 += v3;
                if x.0 % 2 == 1 {
                    x.0 += P;
                }
                x.0 /= 2;
                x
            },
        );
        Coord { v: v % P, chi }
    }

    pub fn from_chi_fp<'a, 'b: 'a>(
        chi: &SylowElem<'a, FpStar<P>>,
        decomp: &SylowDecomp<'b, FpStar<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.invert(decomp).to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord {
            v: (chi + chi_inv).into(),
            chi: Right(chi),
        }
    }

    pub fn from_chi_quad<'a, 'b: 'a>(
        chi: &SylowElem<'a, QuadField<P>>,
        decomp: &SylowDecomp<'b, QuadField<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.invert(decomp).to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord {
            v: (chi + chi_inv).0,
            chi: Left(chi),
        }
    }

    pub fn v(&self) -> u128 {
        self.v
    }

    pub fn chi(&self) -> &Either<QuadNum<P>, FpNum<P>> {
        &self.chi
    }

    pub fn rot<'a>(
        self: &'a Rc<Self>,
        b: &'a Rc<Coord<P>>,
        c: &'a Rc<Coord<P>>,
        fp: &'a QuadField<P>,
    ) -> impl Iterator<Item = (Rc<Coord<P>>, Rc<Coord<P>>)> + 'a {
        std::iter::successors(Some((Rc::clone(b), Rc::clone(c))), move |(y, z)| {
            let tmp = long_multiply(3, self.v(), P);
            let tmp = long_multiply(tmp, z.v(), P);
            let (b_, c_) = (Rc::clone(z), Rc::new(Coord::new(tmp + P - y.v(), fp)));
            if &b_ == b && &c_ == c {
                None
            } else {
                Some((b_, c_))
            }
        })
    }

    /*
    pub fn get_ord(&self, fp: &QuadField) -> Factorization {
        self.chi().as_ref().either(|l| l.order(fp), |r| r.order(fp.pminusone()))
    }
    */
}

impl<const P: u128> PartialEq for Coord<P> {
    fn eq(&self, other: &Self) -> bool {
        self.v == other.v
    }
}
impl<const P: u128> Eq for Coord<P> {}
