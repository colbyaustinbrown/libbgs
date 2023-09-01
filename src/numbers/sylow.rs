use std::fmt;
use std::marker::PhantomData;

use crate::numbers::*;
use crate::util::*;

#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<S, const L: usize, C: SylowDecomposable<S, L>> {
    generators: Vec<C::Elem>,
    _phantom: PhantomData<S>,
}

#[derive(Eq, PartialEq)]
pub struct SylowElem<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> {
    _group: PhantomData<fn(&'a C) -> &'a C>,
    pub coords: Vec<u128>,
    _phantom: PhantomData<S>,
}

pub trait SylowDecomposable<S, const L: usize>: Factored<S, L> + Eq {
    type Elem: GroupElem;

    fn find_sylow_generator(i: usize) -> Self::Elem;

    fn is_sylow_generator(candidate: &Self::Elem, d: (u128, u128)) -> Option<Self::Elem> {
        let pow = Self::Elem::size() / intpow(d.0, d.1, 0);
        let res = candidate.pow(pow);
        if res.pow(intpow(d.0, d.1 - 1, 0)).is_one() {
            None
        } else {
            Some(res)
        }
    }
}

impl<const L: usize, S: Eq, C: SylowDecomposable<S, L>> SylowDecomp<S, L, C> {
    pub fn new() -> SylowDecomp<S, L, C> {
        let generators = (0..L)
            .map(|i| C::find_sylow_generator(i))
            .collect();
        SylowDecomp {
            generators,
            _phantom: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.generators.len()
    }

    pub fn generators(&self) -> &[C::Elem] {
        &self.generators
    }
}

impl<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> Factored<S, L> for SylowElem<'a, S, L, C> {
    const FACTORS: Factorization<L> = <C as Factored<S, L>>::FACTORS;
}

impl<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> SylowDecomposable<S, L> for SylowElem<'a, S, L, C> {
    type Elem = SylowElem<'a, S, L, C>;

    fn find_sylow_generator(i: usize) -> Self::Elem {
        let mut coords = vec![0; L];
        coords[i] = 1;
        SylowElem {
            _group: PhantomData,
            coords,
            _phantom: PhantomData,
        }
    }
}

impl<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> SylowElem<'a, S, L, C> {
    pub fn new(coords: Vec<u128>) -> SylowElem<'a, S, L, C> {
        SylowElem {
            _group: PhantomData,
            coords,
            _phantom: PhantomData,
        }
    }

    pub fn to_product(&self, g: &SylowDecomp<S, L, C>) -> C::Elem {
        (0..g.len())
            .filter(|i| self.coords[*i] > 0)
            .fold(C::Elem::one(), |x, i| {
                let y = g.generators[i].pow(self.coords[i]);
                x.multiply(&y)
            })
    }

    pub fn order(&self, g: &'a SylowDecomp<S, L, C>) -> Factorization<L> {
        let mut prime_powers = [(0,0); L];
        for i in 0..L {
            let mut x = self.clone();
            for j in 0..g.len() {
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

impl<'a, S, const L: usize, C: Eq> GroupElem for SylowElem<'a, S, L, C>
where
    S: Eq,
    C: SylowDecomposable<S, L> + 'a,
{
    fn is_one(&self) -> bool {
        self.coords.iter().all(|x| *x == 0)
    }

    fn multiply(&self, other: &SylowElem<S, L, C>) -> SylowElem<'a, S, L, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .zip(&other.coords)
                .enumerate()
                .map(|(i, (x, y))| (x + y) % <C as Factored<S, L>>::FACTORS.factor(i))
                .collect(),
            _group: PhantomData,
            _phantom: PhantomData,
        }
    }

    fn square(&self) -> SylowElem<'a, S, L, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .enumerate()
                .map(|(i, x)| x * 2 % <C as Factored<S, L>>::FACTORS.factor(i))
                .collect(),
            _group: PhantomData,
            _phantom: PhantomData,
        }
    }

    fn inverse(&self) -> SylowElem<'a, S, L, C> {
        SylowElem {
            coords: self
                .coords
                .iter()
                .enumerate()
                .map(|(i, x)| <C as Factored<S, L>>::FACTORS.factor(i) - x)
                .collect(),
            _group: PhantomData,
            _phantom: PhantomData,
        }
    }

    fn one() -> SylowElem<'a, S, L, C> {
        SylowElem {
            _group: PhantomData,
            coords: vec![0; L],
            _phantom: PhantomData,
        }
    }

    fn size() -> u128 {
        C::Elem::size()
    }
}

impl<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> Clone for SylowElem<'a, S, L, C> {
    fn clone(&self) -> Self {
        SylowElem {
            _group: PhantomData,
            coords: self.coords.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<'a, S: Eq, const L: usize, C: SylowDecomposable<S, L>> fmt::Debug for SylowElem<'a, S, L, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

pub mod tests {
    use super::*;

    // utility method for external tests
    pub fn test_is_generator_small<S, const L: usize, C: SylowDecomposable<S, L>>(x: &C::Elem, d: u128) -> bool {
        let mut y = x.clone();
        for _ in 1..d {
            if y.is_one() {
                return false;
            }
            y = y.multiply(x);
        }
        y.is_one()
    }

    // utility method for external tests
    pub fn test_is_generator_big<S, const L: usize, C: SylowDecomposable<S, L>>(x: &C::Elem, d: (u128, u128)) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(!y.is_one());
            y = y.pow(d.0);
        }
        y = y.pow(d.0);
        assert!(y.is_one());
    }
}
