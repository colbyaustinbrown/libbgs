use std::fmt;
use std::marker::PhantomData;

use crate::numbers::*;
use crate::util::*;

/// A decomposition of a finite cyclic group into the direct sum of its Sylow subgroups.
/// In particular, this group represents the right hand side of the isomorphism
/// $$G \cong \bigoplus_{i = 1}^n \mathbb{Z} / p_i^{t_i} \mathbb{Z}$$
/// where
/// $$|G| = \prod_{i = 1}^n p_i^{t_i}$$
/// and $G$ is a finite cyclic group.
#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<S, const L: usize, C: SylowDecomposable<S, L>> {
    /// Deprecated. Do not access this field directly.
    #[deprecated]
    generators: [C; L],
    _phantom: PhantomData<S>,
}

#[derive(Eq, PartialEq)]
/// An element of the decomposition of a finite cyclic group into the direct sum of its Sylow
/// subgroups.
pub struct SylowElem<S: Eq, const L: usize, C: SylowDecomposable<S, L>> {
    /// The powers on the generators of the Sylow subgroups.
    /// In particular, if an element of a group $G$ with generators $g_1,\ldots,g_n$ is
    /// $$g = \prod_{i = 1}^n g_i^{r_i},$$
    /// then the coordinates of that element are $r_1,\ldots,r_n$.
    pub coords: [u128; L],
    _phantom: PhantomData<(C, S)>,
}

/// Groups that can be decomposed into a direct sum of cyclic Sylow subgroups.
/// In particular, these groups must be finite and cyclic.
pub trait SylowDecomposable<S, const L: usize>: Factored<S, L> + GroupElem + Eq {
    /// Finds a Sylow generator for the Sylow subgroup of prime power index `i`.
    fn find_sylow_generator(i: usize) -> Self;

    /// True if the given element is a generator of the Sylow subgroup of the prime power
    /// represented by `d`.
    fn is_sylow_generator(candidate: &Self, d: (u128, u128)) -> Option<Self> {
        let pow = Self::size() / intpow(d.0, d.1, 0);
        let res = candidate.pow(pow);
        if res.pow(intpow(d.0, d.1 - 1, 0)).is_one() {
            None
        } else {
            Some(res)
        }
    }
}

impl<const L: usize, S: Eq, C: SylowDecomposable<S, L>> SylowDecomp<S, L, C> {
    /// Returns a decomposition for the group.
    /// This method may be expensive because it calls `find_sylow_generator` for each Sylow
    /// subgroup.
    pub fn new() -> SylowDecomp<S, L, C> {
        let generators = (0..L)
            .map(|i| C::find_sylow_generator(i))
            .collect::<Vec<C>>()
            .try_into()
            .unwrap();
        SylowDecomp {
            generators,
            _phantom: PhantomData,
        }
    }

    pub fn generators(&self) -> &[C] {
        &self.generators
    }
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> Factored<S, L> for SylowElem<S, L, C> {
    const FACTORS: Factorization<L> = <C as Factored<S, L>>::FACTORS;
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> SylowDecomposable<S, L>
    for SylowElem<S, L, C>
{
    fn find_sylow_generator(i: usize) -> Self {
        let mut coords = [0; L];
        coords[i] = 1;
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> SylowElem<S, L, C> {
    /// Returns an element of the Sylow decomposition with the given coordinates.
    pub fn new(coords: [u128; L]) -> SylowElem<S, L, C> {
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    /// Returns the element of the original group with the given coordinates.
    pub fn to_product(&self, g: &SylowDecomp<S, L, C>) -> C {
        (0..L)
            .filter(|i| self.coords[*i] > 0)
            .fold(C::one(), |x, i| {
                let y = g.generators[i].pow(self.coords[i]);
                x.multiply(&y)
            })
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub fn order(&self) -> Factorization<L> {
        let mut prime_powers = [(0, 0); L];
        for i in 0..L {
            let mut x = self.clone();
            for j in 0..L {
                if j == i {
                    continue;
                }
                x = x.pow(<C as Factored<S, L>>::FACTORS.factor(j));
            }

            let mut r = 0;
            while !x.is_one() {
                x = x.pow(<C as Factored<S, L>>::FACTORS[i].0);
                r += 1;
            }
            prime_powers[i] = (<C as Factored<S, L>>::FACTORS[i].0, r);
        }
        Factorization::new(prime_powers)
    }
}

impl<S, const L: usize, C: Eq> GroupElem for SylowElem<S, L, C>
where
    S: Eq,
    C: SylowDecomposable<S, L>,
{
    fn is_one(&self) -> bool {
        self.coords.iter().all(|x| *x == 0)
    }

    fn multiply(&self, other: &SylowElem<S, L, C>) -> SylowElem<S, L, C> {
        let mut coords = self.coords;
        for i in 0..L {
            coords[i] = (coords[i] + other.coords[i]) % <C as Factored<S, L>>::FACTORS.factor(i);
        }
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    fn square(&self) -> SylowElem<S, L, C> {
        let mut coords = self.coords;
        for i in 0..L {
            coords[i] = (coords[i] * 2) % <C as Factored<S, L>>::FACTORS.factor(i);
        }
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    fn inverse(&self) -> SylowElem<S, L, C> {
        let mut coords = self.coords;
        for i in 0..L {
            coords[i] = <C as Factored<S, L>>::FACTORS.factor(i) - coords[i];
        }
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    fn one() -> SylowElem<S, L, C> {
        SylowElem {
            coords: [0; L],
            _phantom: PhantomData,
        }
    }

    fn size() -> u128 {
        C::size()
    }
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> Clone for SylowElem<S, L, C> {
    fn clone(&self) -> Self {
        SylowElem {
            coords: self.coords,
            _phantom: PhantomData,
        }
    }
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> fmt::Debug for SylowElem<S, L, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

/// Utility methods for use in other tests.
/// These methods should probably not be used outside of this crate.
pub mod tests {
    use super::*;

    /// True if `x` is of order `d`, False otherwise.
    /// Expensive and intended for use only with small values of `d`.
    pub fn test_is_generator_small<S, const L: usize, C: SylowDecomposable<S, L>>(
        x: &C,
        d: u128,
    ) -> bool {
        let mut y = x.clone();
        for _ in 1..d {
            if y.is_one() {
                return false;
            }
            y = y.multiply(x);
        }
        y.is_one()
    }

    /// True if `x` is not of order prime power dividing `d`, but is a prime power of `d`.
    /// Much cheaper than `test_is_generator_small`, but may return a false positive.
    pub fn test_is_generator_big<S, const L: usize, C: SylowDecomposable<S, L>>(
        x: &C,
        d: (u128, u128),
    ) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(!y.is_one());
            y = y.pow(d.0);
        }
        y = y.pow(d.0);
        assert!(y.is_one());
    }
}
