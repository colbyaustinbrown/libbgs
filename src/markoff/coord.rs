use either::*;

use crate::numbers::*;

/// A coordinate for a Markoff triple.
/// May represent any of $a$, $b$, or $c$ in a Markoff triple $(a, b, c)$.
/// This is a single field struct containing only an `FpNum<P>` for prime `P`.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Coord<const P: u128>(pub FpNum<P>);

#[derive(PartialEq, Eq)]
/// The three types of conic sections an orbit could lie inside.
pub enum Conic {
    /// An orbit of order dividing $p - 1$ (and not equal to 2).
    Hyperbola,
    /// An orbit of order dividing $p + 1$ (and not equal to 2).
    Ellipse,
    /// An orbit of order exactly 1 or 2.
    Parabola,
}

impl<const P: u128> Coord<P> {
    /// Returns an element $\chi$ such that, for a coordinate $a$, $a = \chi + \chi^{-1}$.
    /// If $a$ is a quadratic residue modulo `P`, then $\chi \in \mathbb{F}\_p$, and the result
    /// will be a `Right<FpNum<P>>`. Otherwise, $\chi \in \mathbb{F}\_{p^2}$, and the result will
    /// be a `Left<QuadNum<P>>`.
    pub fn to_chi(&self) -> Either<QuadNum<P>, FpNum<P>> {
        let disc = self.0.pow(2) - FpNum::from(4);
        let two_inv = FpNum::from(2).inverse();
        QuadNum::int_sqrt_either(disc).map_either(
            |x| (QuadNum::<P>::from(self.0) + x) * QuadNum::<P>::from(two_inv),
            |x| (self.0 + x) * two_inv,
        )
    }

    /// Returns an iterator yielding the coordinates $(b, c)$ contained in the orbit with fixed coordinate
    /// $a$ (the coordinate on which `rot` is called), beginning with $(a, b, c)$.
    pub fn rot(self, b: Coord<P>, c: Coord<P>) -> impl Iterator<Item = (Coord<P>, Coord<P>)> {
        std::iter::successors(Some((b, c)), move |(y, z)| {
            let (b_, c_) = (*z, self.0 * z.0 + FpNum::from(P) - y.0);
            if b_ == b && c_ == c.0 {
                None
            } else {
                Some((b_, Coord(c_)))
            }
        })
    }

    /// Returns an iterator yielding the values $b$ such that $(a, b, c)$ is a Markoff triple for
    /// some value $c$.
    pub fn part(self, b: Coord<P>) -> Option<Box<dyn Iterator<Item = Coord<P>>>> {
        let a = self.0;
        let Some(disc) = (a * a * b.0 * b.0 - 4 * (a * a + b.0 * b.0)).int_sqrt() else {
            return None;
        };
        let c = (a * b.0 + disc) * FpNum::from(2).inverse();
        Some(Box::new(self.rot(b, Coord(c)).map(|x| x.1)))
    }

    /// Returns the order of the map $\text{rot}\_a$, that is, $\lvert \langle \text{rot}\_a \rangle \rvert$, along with the type of [`Conic`] that it is.
    pub fn rot_order<S1, S2>(&self) -> (u128, Conic)
    where
        FpNum<P>: Factor<S1>,
        QuadNum<P>: Factor<S2>,
    {
        match self
            .to_chi()
            .as_ref()
            .map_either(|l| l.order(), |r| r.order())
        {
            Left(1) | Right(1) => (1, Conic::Parabola),
            Left(2) | Right(2) => (2, Conic::Parabola),
            Left(d) => (d, Conic::Hyperbola),
            Right(d) => (d, Conic::Ellipse),
        }
    }

    /// Returns an upper bound on the endgame breakpoint.
    /// That is, every triple with order larger than the value returned by this method is
    /// guarenteed to lie in the endgame.
    pub fn endgame<S>() -> (u128, u128)
    where
        FpNum<P>: Factor<S>,
        QuadNum<P>: Factor<S>,
    {
        let tmp = 8.0 * (P as f64).sqrt();
        let hyper = tmp * ((P - 1) * FpNum::FACTORS.tau()) as f64;
        let hyper = hyper / (FpNum::FACTORS.phi() as f64);
        let ellip = tmp * ((P + 1) * QuadNum::FACTORS.tau()) as f64;
        let ellip = ellip / (QuadNum::FACTORS.phi() as f64);
        let res = (hyper.ceil() as u128, ellip.ceil() as u128);
        (std::cmp::min(res.0, P - 2), std::cmp::min(res.1, P))
    }
}

impl<const P: u128> From<u128> for Coord<P> {
    fn from(src: u128) -> Coord<P> {
        Coord(FpNum::from(src))
    }
}

impl<const P: u128> From<Coord<P>> for u128 {
    fn from(src: Coord<P>) -> u128 {
        u128::from(src.0)
    }
}

/// Common trait for the `from_chi` and `from_chi_conj` methods to be defined on both `FpNum` and
/// `QuadNum`.
pub trait FromChi<S, const P: u128, C> 
where
    C: SylowDecomposable<S>
{
    /// Returns the coordinate $a = \chi + \chi^{-1}$.
    fn from_chi<const L: usize>(
        chi: &SylowElem<S, L, C>,
        decomp: &SylowDecomp<S, L, C>,
    ) -> Self;

    /// Returns the coordinate $a = \chi - \chi^{-1}$.
    fn from_chi_conj<const L: usize>(
        chi: &SylowElem<S, L, C>,
        decomp: &SylowDecomp<S, L, C>,
    ) -> Self;
}

impl<S, const P: u128> FromChi<S, P, FpNum<P>> for Coord<P>
where
    FpNum<P>: Factor<S>,
{
    fn from_chi<const L: usize>(
        chi: &SylowElem<S, L, FpNum<P>>,
        decomp: &SylowDecomp<S, L, FpNum<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord(chi + chi_inv)
    }

    fn from_chi_conj<const L: usize>(
        chi: &SylowElem<S, L, FpNum<P>>,
        decomp: &SylowDecomp<S, L, FpNum<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord(chi - chi_inv)
    }
}

impl<S, const P: u128> FromChi<S, P, QuadNum<P>> for Coord<P>
where
    QuadNum<P>: Factor<S>,
{
    fn from_chi<const L: usize>(
        chi: &SylowElem<S, L, QuadNum<P>>,
        decomp: &SylowDecomp<S, L, QuadNum<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord((chi + chi_inv).0)
    }

    fn from_chi_conj<const L: usize>(
        chi: &SylowElem<S, L, QuadNum<P>>,
        decomp: &SylowDecomp<S, L, QuadNum<P>>,
    ) -> Coord<P> {
        let chi_inv = chi.inverse().to_product(decomp);
        let chi = chi.to_product(decomp);
        Coord((chi - chi_inv).0)
    }
}
