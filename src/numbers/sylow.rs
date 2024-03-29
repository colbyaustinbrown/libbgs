use std::marker::PhantomData;

use crate::numbers::*;
use libbgs_util::*;

/// A decomposition of a finite cyclic group into the direct sum of its Sylow subgroups.
/// In particular, this group represents the right hand side of the isomorphism
/// $$G \cong \bigoplus_{i = 1}^n \mathbb{Z} / p_i^{t_i} \mathbb{Z}$$
/// where
/// $$|G| = \prod_{i = 1}^n p_i^{t_i}$$
/// and $G$ is a finite cyclic group.
pub struct SylowDecomp<S, const L: usize, C: SylowDecomposable<S>> {
    precomputed: [[C; 256]; L],
    generators_powered: [C; L],
    _phantom: PhantomData<S>,
}

/// An element of the decomposition of a finite cyclic group into the direct sum of its Sylow
/// subgroups.
pub struct SylowElem<S, const L: usize, C: SylowDecomposable<S>> {
    /// The powers on the generators of the Sylow subgroups.
    /// In particular, if an element of a group $G$ with generators $g_1,\ldots,g_n$ is
    /// $$g = \prod_{i = 1}^n g_i^{r_i},$$
    /// then the coordinates of that element are $r_1,\ldots,r_n$.
    pub coords: [u128; L],
    _phantom: PhantomData<(C, S)>,
}

/// Groups that can be decomposed into a direct sum of cyclic Sylow subgroups.
/// In particular, these groups must be finite and cyclic.
pub trait SylowDecomposable<S>: Factor<S> + GroupElem + Eq {
    /// Finds a Sylow generator for the Sylow subgroup of prime power index `i`.
    fn find_sylow_generator(i: usize) -> Self;

    /// True if the given element is a generator of the Sylow subgroup of the prime power
    /// represented by `d`.
    fn is_sylow_generator(candidate: &Self, d: (u128, usize)) -> Option<Self> {
        let pow = Self::SIZE / intpow::<0>(d.0, d.1 as u128);
        let res = candidate.pow(pow);
        if res.pow(intpow::<0>(d.0, (d.1 - 1) as u128)) == Self::ONE {
            None
        } else {
            Some(res)
        }
    }

    /// Returns the number of elements of a particular order.
    /// The argument is the powers of the prime factors of the group's order.
    fn count_elements_of_order(ds: &[usize]) -> u128 {
        let mut total = 1;
        for (d, (p, t)) in ds.iter().zip(Self::FACTORS.factors()) {
            if *d > *t {
                return 0;
            } else if *d > 0 {
                let tmp = intpow::<0>(*p, (*d - 1) as u128);
                total *= tmp * *p - tmp;
            }
        }
        total
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> SylowDecomp<S, L, C> {
    /// Returns a decomposition for the group.
    /// This method may be expensive because it calls `find_sylow_generator` for each Sylow
    /// subgroup.
    pub fn new() -> SylowDecomp<S, L, C> {
        let mut generators_powered = [C::ONE; L];
        #[allow(deprecated)]
        let mut precomputed = [C::ONE_256; L];
        let mut i = 0;
        while i < L {
            let x = C::find_sylow_generator(i);
            let mut g = C::ONE;
            let mut j = 0;
            while j < 256 {
                precomputed[i][j] = g.clone();
                g = g.multiply(&x);
                j += 1;
            }
            generators_powered[i] = g;
            i += 1;
        }
        SylowDecomp {
            precomputed,
            generators_powered,
            _phantom: PhantomData,
        }
    }

    /// Get the generators for decomposition.
    /// The index of each generator corresponds to the index of the prime power in the
    /// factorization. That is, if the prime power at index `i` of the factorization is $(p, t)$,
    /// then the generator at index `i` of the array returned by the `generators` method is a
    /// generator of the Sylow subgroup of order $p^t$.
    pub fn generator(&self, i: usize) -> &C {
        &self.precomputed[i][1]
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> Factor<S> for SylowElem<S, L, C> {
    const FACTORS: Factorization = <C as Factor<S>>::FACTORS;
}

impl<S, const L: usize, C: SylowDecomposable<S>> SylowDecomposable<S> for SylowElem<S, L, C> {
    fn find_sylow_generator(i: usize) -> Self {
        let mut coords = [0; L];
        coords[i] = 1;
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> SylowElem<S, L, C> {
    /// Returns an element of the Sylow decomposition with the given coordinates.
    pub const fn new(coords: [u128; L]) -> SylowElem<S, L, C> {
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    /// Returns the element of the original group with the given coordinates.
    pub fn to_product(&self, g: &SylowDecomp<S, L, C>) -> C {
        (0..L).filter(|i| self.coords[*i] > 0).fold(C::ONE, |x, i| {
            let mut y = g.precomputed[i][(self.coords[i] & 0xFF) as usize].clone();
            if self.coords[i] > 0xFF {
                y = y.multiply(&g.generators_powered[i].pow(self.coords[i] >> 8));
            }
            x.multiply(&y)
        })
    }

    /// Returns the positive integer represented by this `Factorization`.
    pub fn order(&self) -> u128 {
        let mut res = 1;
        for i in 0..L {
            let mut x = *self;
            for j in 0..L {
                if j == i {
                    continue;
                }
                x = x.pow(C::FACTORS.factor(j));
            }

            let mut r = 0;
            while x != Self::ONE {
                x = x.pow(C::FACTORS[i].0);
                r += 1;
            }
            res *= intpow::<0>(C::FACTORS[i].0, r);
        }
        res
    }
}

impl<S, const L: usize, C: Eq> GroupElem for SylowElem<S, L, C>
where
    C: SylowDecomposable<S>,
{
    const ONE: Self = SylowElem {
        coords: [0; L],
        _phantom: PhantomData,
    };

    const SIZE: u128 = C::SIZE;

    fn multiply(&self, other: &SylowElem<S, L, C>) -> SylowElem<S, L, C> {
        let mut coords = self.coords;
        for i in 0..L {
            coords[i] = (coords[i] + other.coords[i]) % C::FACTORS.factor(i);
        }
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    fn inverse(&self) -> SylowElem<S, L, C> {
        let mut coords = self.coords;
        for i in 0..L {
            coords[i] = C::FACTORS.factor(i) - coords[i];
        }
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }
}

impl<S, const L: usize, C: SylowDecomposable<S>> PartialEq for SylowElem<S, L, C> {
    fn eq(&self, other: &Self) -> bool {
        self.coords == other.coords
    }
}
impl<S, const L: usize, C: SylowDecomposable<S>> Eq for SylowElem<S, L, C> {}

impl<S, const L: usize, C: SylowDecomposable<S>> Clone for SylowElem<S, L, C> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<S, const L: usize, C: SylowDecomposable<S>> Copy for SylowElem<S, L, C> {}

impl<S, const L: usize, C: SylowDecomposable<S>> std::fmt::Debug for SylowElem<S, L, C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.coords.fmt(f)
    }
}

/// Utility methods for use in other tests.
/// These methods should probably not be used outside of this crate.
pub mod tests {
    use super::*;

    /// True if `x` is of order `d`, False otherwise.
    /// Expensive and intended for use only with small values of `d`.
    pub fn test_is_generator_small<S, const L: usize, C: SylowDecomposable<S>>(
        x: &C,
        d: usize,
    ) -> bool {
        let mut y = x.clone();
        for _ in 1..d {
            if y == C::ONE {
                return false;
            }
            y = y.multiply(x);
        }
        y == C::ONE
    }

    /// True if `x` is not of order prime power dividing `d`, but is a prime power of `d`.
    /// Much cheaper than `test_is_generator_small`, but may return a false positive.
    pub fn test_is_generator_big<S, const L: usize, C: SylowDecomposable<S>>(
        x: &C,
        d: (u128, usize),
    ) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(y != C::ONE);
            y = y.pow(d.0);
        }
        y = y.pow(d.0);
        assert!(y == C::ONE);
    }
}
