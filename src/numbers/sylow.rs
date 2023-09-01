use std::fmt;
use std::marker::PhantomData;

use crate::numbers::*;
use crate::util::*;

#[derive(PartialEq, Eq, Debug)]
pub struct SylowDecomp<S, const L: usize, C: SylowDecomposable<S, L>> {
    generators: [C; L],
    _phantom: PhantomData<S>,
}

#[derive(Eq, PartialEq)]
pub struct SylowElem<S: Eq, const L: usize, C: SylowDecomposable<S, L>> {
    pub coords: [u128; L],
    _phantom: PhantomData<(C, S)>,
}

pub trait SylowDecomposable<S, const L: usize>: Factored<S, L> + GroupElem + Eq {
    fn find_sylow_generator(i: usize) -> Self;

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

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> SylowDecomposable<S, L> for SylowElem<S, L, C> {
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
    pub fn new(coords: [u128; L]) -> SylowElem<S, L, C> {
        SylowElem {
            coords,
            _phantom: PhantomData,
        }
    }

    pub fn to_product(&self, g: &SylowDecomp<S, L, C>) -> C {
        (0..L)
            .filter(|i| self.coords[*i] > 0)
            .fold(C::one(), |x, i| {
                let y = g.generators[i].pow(self.coords[i]);
                x.multiply(&y)
            })
    }

    pub fn order(&self) -> Factorization<L> {
        let mut prime_powers = [(0,0); L];
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
            coords: self.coords.clone(),
            _phantom: PhantomData,
        }
    }
}

impl<S: Eq, const L: usize, C: SylowDecomposable<S, L>> fmt::Debug for SylowElem<S, L, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.coords.fmt(f)
    }
}

pub mod tests {
    use super::*;

    // utility method for external tests
    pub fn test_is_generator_small<S, const L: usize, C: SylowDecomposable<S, L>>(x: &C, d: u128) -> bool {
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
    pub fn test_is_generator_big<S, const L: usize, C: SylowDecomposable<S, L>>(x: &C, d: (u128, u128)) {
        let mut y = x.clone();
        for _ in 0..d.1 {
            assert!(!y.is_one());
            y = y.pow(d.0);
        }
        y = y.pow(d.0);
        assert!(y.is_one());
    }
}
