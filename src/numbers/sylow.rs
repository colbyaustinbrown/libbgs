use std::fmt;
use std::marker::PhantomData;

pub use crate::numbers::group::*;

use crate::numbers::factorization::*;

/// A decomposition of a finite cyclic group into the direct sum of its Sylow subgroups.
/// In particular, this group represents the right hand side of the isomorphism
/// $$G \cong \bigoplus_{i = 1}^n \mathbb{Z} / p_i^{t_i} \mathbb{Z}$$
/// where
/// $$|G| = \prod_{i = 1}^n p_i^{t_i}$$
/// and $G$ is a finite cyclic group.
#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<'a, C: SylowDecomposable> {
    /// Deprecated. Do not access this field directly.
    #[deprecated]
    pub parent: &'a C,
    /// Deprecated. Do not access this field directly.
    #[deprecated]
    pub fact: Factorization,
    /// Deprecated. Do not access this field directly.
    #[deprecated]
    pub generators: Vec<C::Elem>,
}

/// An element of the decomposition of a finite cyclic group into the direct sum of its Sylow
/// subgroups.
#[derive(Eq, PartialEq)]
pub struct SylowElem<'a, C: SylowDecomposable> {
    _group: PhantomData<fn(&'a C) -> &'a C>,
    /// The powers on the generators of the Sylow subgroups.
    /// In particular, if an element of a group $G$ with generators $g_1,\ldots,g_n$ is
    /// $$g = \prod_{i = 1}^n g_i^{r_i},$$
    /// then the coordinates of that element are $r_1,\ldots,r_n$.
    pub coords: Vec<u128>,
}

/// Groups that can be decomposed into a direct sum of cyclic Sylow subgroups.
/// In particular, these groups must be finite and cyclic.
pub trait SylowDecomposable: Group {
    /// Finds a Sylow generator for the Sylow subgroup of prime power index `i`.
    fn find_sylow_generator(&self, i: usize, fact: &Factorization) -> Self::Elem;

    /// True if the given element is a generator of the Sylow subgroup of the prime power
    /// represented by `d`.
    fn is_sylow_generator(&self, candidate: &Self::Elem, d: (u128, u128)) -> Option<Self::Elem> {
        let pow = self.size() / intpow(d.0, d.1, 0);
        let res = candidate.pow(pow, self);
        if res.pow(intpow(d.0, d.1 - 1, 0), self).is_one(self) {
            None
        } else {
            Some(res)
        }
    }
}

impl<'a, C: SylowDecomposable> SylowDecomp<'a, C> {
    /// Returns a decomposition for the group.
    /// This method may be expensive because it calls `find_sylow_generator` for each Sylow
    /// subgroup.
    pub fn new(parent: &C, fact: Factorization) -> SylowDecomp<C> {
        let length = fact.len();
        let generators = (0..length)
            .map(|i| parent.find_sylow_generator(i, &fact))
            .collect();
        SylowDecomp {
            parent,
            fact,
            generators,
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }

    fn len(&self) -> usize {
        self.generators.len()
    }
}

impl<'a, C: SylowDecomposable> Group for SylowDecomp<'a, C> {
    type Elem = SylowElem<'a, C>;

    fn one(&self) -> SylowElem<'a, C> {
        SylowElem {
            _group: PhantomData,
            coords: vec![0; self.factors().len()],
        }
    }

    fn size(&self) -> u128 {
        self.parent.size()
    }
}

impl<'a, C: SylowDecomposable> Factored for SylowDecomp<'a, C> {
    fn factors(&self) -> &Factorization {
        &self.fact
    }
}

impl<'a, C: SylowDecomposable> SylowDecomposable for SylowDecomp<'a, C> {
    fn find_sylow_generator(&self, i: usize, _: &Factorization) -> Self::Elem {
        let mut coords = vec![0; self.factors().len()];
        coords[i] = 1;
        SylowElem {
            _group: PhantomData,
            coords,
        }
    }
}

impl<'a, C: SylowDecomposable> SylowElem<'a, C> {
    /// Returns an element of the Sylow decomposition with the given coordinates.
    pub fn new(coords: Vec<u128>) -> SylowElem<'a, C> {
        SylowElem {
            _group: PhantomData,
            coords,
        }
    }

    /// Returns the element of the original group with the given coordinates.
    pub fn to_product(&self, g: &SylowDecomp<C>) -> C::Elem {
        (0..g.len())
            .filter(|i| self.coords[*i] > 0)
            .fold(g.parent.one(), |x, i| {
                let y = g.generators[i].pow(self.coords[i], g.parent);
                x.multiply(&y, g.parent)
            })
    }

    /// Returns the order of this element, that is, the smallest positive power to which this
    /// element raised is equal to one.
    pub fn order(&self, g: &'a SylowDecomp<C>) -> Factorization {
        let prime_powers: Vec<(u128, u128)> = (0..g.len())
            .map(|i| {
                let mut x = self.clone();
                for j in 0..g.len() {
                    if j == i {
                        continue;
                    }
                    x = x.pow(g.fact.factor(j), g);
                }

                let mut r = 0;
                while !x.is_one(g) {
                    x = x.pow(g.fact[i].0, g);
                    r += 1;
                }
                (g.fact[i].0, r)
            })
            .collect();
        Factorization::new(prime_powers)
    }
}

impl<'a, C> GroupElem for SylowElem<'a, C>
where
    C: SylowDecomposable + 'a,
{
    type Group = SylowDecomp<'a, C>;

    fn is_one(&self, _: &SylowDecomp<C>) -> bool {
        self.coords.iter().all(|x| *x == 0)
    }

    fn multiply(&self, other: &SylowElem<C>, g: &SylowDecomp<C>) -> SylowElem<'a, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .enumerate()
                .map(|(i, (x, y))| (x + y) % g.fact.factor(i))
                .collect(),
            _group: PhantomData,
        }
    }

    fn square(&self, g: &SylowDecomp<C>) -> SylowElem<'a, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .enumerate()
                .map(|(i, x)| x * 2 % g.fact.factor(i))
                .collect(),
            _group: PhantomData,
        }
    }

    fn invert(&self, g: &Self::Group) -> SylowElem<'a, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .enumerate()
                .map(|(i, x)| g.fact.factor(i) - x)
                .collect(),
            _group: PhantomData,
        }
    }
}

impl<'a, C: SylowDecomposable> Clone for SylowElem<'a, C> {
    fn clone(&self) -> Self {
        SylowElem {
            _group: PhantomData,
            coords: self.coords.clone(),
        }
    }
}

impl<'a, C: SylowDecomposable> fmt::Debug for SylowElem<'a, C> {
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
    pub fn test_is_generator_small<C: SylowDecomposable>(x: &C::Elem, d: u128, g: &C) -> bool {
        let mut y = x.clone();
        for _ in 1..d {
            if y.is_one(g) {
                return false;
            }
            y = y.multiply(x, g);
        }
        y.is_one(g)
    }

    /// True if `x` is not of order prime power dividing `d`, but is a prime power of `d`.
    /// Much cheaper than `test_is_generator_small`, but may return a false positive.
    pub fn test_is_generator_big<C: SylowDecomposable>(x: &C::Elem, d: (u128, u128), g: &C) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(!y.is_one(g));
            y = y.pow(d.0, g);
        }
        y = y.pow(d.0, g);
        assert!(y.is_one(g));
    }
}
