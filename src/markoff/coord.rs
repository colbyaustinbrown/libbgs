use either::*;

use crate::numbers::*;

/// A coordinate for a Markoff triple.
/// May represent any of $a$, $b$, or $c$ in a Markoff triple $(a, b, c)$.
/// This is a single field struct containing only an `FpNum<P>` for prime `P`.
// pub struct Coord<const P: u128>(FpNum<P>);
pub type Coord<const P: u128> = FpNum<P>;

impl<const P: u128> Coord<P> {
    /// Returns an element $\chi$ such that, for a coordinate $a$, $a = \chi + \chi^{-1}$.
    /// If $a$ is a quadratic residue modulo `P`, then $\chi \in \mathbb{F}\_p$, and the result
    /// will be a `Right<FpNum<P>>`. Otherwise, $\chi \in \mathbb{F}\_{p^2}$, and the result will
    /// be a `Left<QuadNum<P>>`.
    pub fn to_chi(&self) -> Either<QuadNum<P>, FpNum<P>> {
        let v3 = self.multiply(&FpNum::from(3));
        let disc = v3.pow(2);
        let disc = disc + FpNum::from(P) - FpNum::from(4);
        let two_inv = FpNum::from(2).inverse();
        QuadNum::int_sqrt_either(disc).map_either(
            |mut x| {
                x.0 += v3;
                x.0 = x.0.multiply(&two_inv);
                x.1 = x.1.multiply(&two_inv);
                x
            },
            |mut x| {
                x += v3;
                x.multiply(&two_inv)
            },
        )
    }

    /// Returns the coordinate $a = \chi + \chi^{-1}$, where $\chi \in \mathbb{F}_p$.
    pub fn from_chi_fp<S, const L: usize>(
        chi: &SylowElem<S, L, FpNum<P>>,
        decomp: &SylowDecomp<S, L, FpNum<P>>,
    ) -> Coord<P>
    where
        FpNum<P>: Factor<S, L>,
    {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);

        // We use the non-normalized equation:
        // x^2 + y^2 + z^2 - xyz = 0
        chi + chi_inv
    }

    /// Returns the coordinate $a = \chi + \chi^{-1}$, where $\chi \in \mathbb{F}_{p^2}$.
    pub fn from_chi_quad<S, const L: usize>(
        chi: &SylowElem<S, L, QuadNum<P>>,
        decomp: &SylowDecomp<S, L, QuadNum<P>>,
    ) -> Coord<P>
    where
        QuadNum<P>: Factor<S, L>,
    {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);

        // We use the non-normalized equation:
        // x^2 + y^2 + z^2 - xyz = 0
        FpNum::from((chi + chi_inv).0)
    }

    /// Returns an iterator yielding the coordinates $(b, c)$ contained in the orbit with fixed coordinate
    /// $a$ (the coordinate on which `rot` is called), beginning with $(a, b, c)$.
    pub fn rot(self, b: Coord<P>, c: Coord<P>) -> impl Iterator<Item = (Coord<P>, Coord<P>)> {
        std::iter::successors(Some((b, c)), move |(y, z)| {
            let (b_, c_) = (*z, self * z + FpNum::from(P) - *y);
            if b_ == b && c_ == c {
                None
            } else {
                Some((b_, c_))
            }
        })
    }

    /// Returns the order of the map $\text{rot}\_a$, that is, $\lvert \langle \text{rot}\_a \rangle \rvert$.
    pub fn get_ord<
        S1, 
        S2, 
        const L1: usize, 
        const L2: usize
    >(&self) -> u128
    where
        FpNum<P>: Factor<S1, L1>,
        QuadNum<P>: Factor<S2, L2>,
    {
        self.to_chi()
            .as_ref()
            .either(|l| l.order().value(), |r| r.order().value())
    }
}
