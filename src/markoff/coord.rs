use either::*;

use crate::numbers::*;
use crate::util::*;

#[derive(PartialEq, Clone, Copy, Eq, Debug)]
pub struct Coord<const P: u128>(FpNum<P>);

impl<const P: u128> Coord<P> {
    pub fn to_chi(&self, fp: &QuadField<P>) -> Either<QuadNum<P>, FpNum<P>> {
        let v3 = long_multiply(self.0.0, 3, P);
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
        chi
    }

    pub fn from_chi_fp<'a, 'b: 'a, S, const L: usize>(
        chi: &SylowElem<'a, S, L, FpStar<P>>,
        decomp: &SylowDecomp<'b, S, L, FpStar<P>>,
    ) -> Coord<P> 
    where
        S: Eq,
        FpStar<P>: Factored<S, L>
    {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);

        // We use the non-normalized equation:
        // x^2 + y^2 + z^2 - xyz = 0
        Coord(chi + chi_inv)
    }

    pub fn from_chi_quad<'a, 'b: 'a, S, const L: usize>(
        chi: &SylowElem<'a, S, L, QuadField<P>>,
        decomp: &SylowDecomp<'b, S, L, QuadField<P>>,
    ) -> Coord<P> 
    where
        S: Eq,
        QuadField<P>: Factored<S, L>
    {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);

        // We use the non-normalized equation:
        // x^2 + y^2 + z^2 - xyz = 0
        Coord(FpNum::from((chi + chi_inv).0))
    }

    pub fn rot<'a>(
        self,
        b: Coord<P>,
        c: Coord<P>,
    ) -> impl Iterator<Item = (Coord<P>, Coord<P>)> {
        std::iter::successors(Some((b, c)), move |(y, z)| {
            let (b_, c_) = (*z, Coord::from(self.0 * z.0 + P.into() - y.0));
            if b_ == b && c_ == c {
                None
            } else {
                Some((b_, c_))
            }
        })
    }

    pub fn get_ord<const L: usize>(
        &self,
        fp: &QuadField<P>,
        minusonesize: &Factorization<L>,
        plusonesize: &Factorization<L>,
    ) -> Factorization<L> {
        self.to_chi(fp).as_ref().either(
            |l| l.order(plusonesize),
            |r| r.order(minusonesize),
        )
    }
}

impl<const P: u128> From<u128> for Coord<P> {
    fn from(value: u128) -> Coord<P> {
        Coord(FpNum::from(value))
    }
}

impl<const P: u128> From<FpNum<P>> for Coord<P> {
    fn from(value: FpNum<P>) -> Coord<P> {
        Coord(value)
    }
}

impl<const P: u128> From<Coord<P>> for u128 {
    fn from(value: Coord<P>) -> u128 {
        u128::from(value.0)
    }
}

impl<const P: u128> std::fmt::Display for Coord<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        self.0.fmt(f)
    }
}
