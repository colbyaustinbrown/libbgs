use std::rc::Rc;
// use either::*;

use crate::factorization::*;
use crate::util::*;
use crate::semigroup::*;
use crate::sylow::*;
use crate::fp::*;

#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldExt {
    pminusone: Rc<FpStar>,
    pplusone: Factorization,
    r: u128
}

type QuadSubgroup = QuadFieldExt;

impl Semigroup for QuadSubgroup {
    type Elem = QuadNumber;
    fn size(&self) -> &Factorization {
        &self.pplusone
    }
    fn one(self: &Rc<Self>) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0: 1,
            a1: 0
        }
    }

}

impl QuadFieldExt {
    pub fn new(pminusone: &Rc<FpStar>, pplusone: Factorization, r: u128) -> QuadFieldExt {
        QuadFieldExt {
            pminusone: Rc::clone(pminusone),
            pplusone,
            r
        }
    }
    pub fn make(pminusone: &Rc<FpStar>, pplusone: Factorization) -> QuadFieldExt {
        let p = pplusone.value() - 1;
        let r = if p % 4 == 3 {
            p - 1
        } else if p % 8 == 3 || p % 8 == 5 {
            2
        } else {
            let mut res = 0;
            for i in 0..p {
                let a = standard_affine_shift(p, i);
                if intpow(a, (p - 1) / 2, p) == p - 1 {
                    res = a;
                    break;
                }
            }
            res 
        };
        QuadFieldExt::new(pminusone, pplusone, r)
    }

    pub fn change_r(&self, r: u128) -> QuadFieldExt {
        QuadFieldExt {
            r: r,
            pminusone: Rc::clone(&self.pminusone),
            pplusone: self.pplusone.clone()
        }
    }

    fn steinitz(self: &Rc<Self>, i: u128) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0: i % self.p(),
            a1: i / self.p()
        }
    }

    fn from_ints(self: &Rc<Self>, a0: u128, a1: u128) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(self),
            a0,
            a1
        }
    }

    pub fn p(&self) -> u128 {
        self.pminusone.value() + 1
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct QuadNumber {
    subgroup: Rc<QuadSubgroup>,
    a0: u128,
    a1: u128
}

impl QuadNumber {
    fn is_zero(&self) -> bool {
        self.a0 == 0 && self.a1 == 0
    }

    fn from_int(subgroup: Rc<QuadSubgroup>, a0: u128) -> QuadNumber {
        QuadNumber {
            subgroup,
            a0,
            a1: 0
        }
    }
}

impl SemigroupElem for QuadNumber {
    type Group = QuadSubgroup;
    fn is_one(&self) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn group(&self) -> &Rc<QuadSubgroup> {
        &self.subgroup
    }

    fn multiply(&mut self, other: &QuadNumber) {
        let f = &self.subgroup;
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r, p), p)) % p;
        self.a1 = (long_multiply(self.a1, other.a0, p) + long_multiply(a0_old, other.a1, p)) % p;
    }

    fn square(&mut self) {
        let f = &self.subgroup;
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r, p), p)) % p;
        self.a1 = (long_multiply(self.a1, a0_old, p) + long_multiply(a0_old, self.a1, p)) % p;
    }
}

impl SylowDecomposable for QuadSubgroup {
    fn find_sylow_generator(self: &Rc<Self>, i: usize) -> QuadNumber {
        let pow = self.pminusone.value();
        // should be self.p * self.p, but maybe this works?
        (1..self.p() * 2)
            .map(|i| {
                let j = standard_affine_shift(self.p() * 2, i);
                let mut p = self.steinitz(j);
                p.pow(pow);
                p
            })
            .find_map(|c| self.is_sylow_generator(&c, i))
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sylow::tests::*;

    fn p7() -> QuadFieldExt {
        QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2,1), (3,1)]
            )),
            Factorization ::new(
                vec![(2,3)]
            )
        )
    }

    #[test]
    fn one_is_one() {
        let one = Rc::new(p7()).one();
        assert!(one.is_one());
    }

    #[test]
    fn calculates_r_as_nonresidue() {
        let f = p7();
        for i in 2..f.p() {
            assert_ne!((i * i) % f.p(), f.r);
        }
    }

    #[test]
    fn powers_up() {
        let mut x = QuadNumber {
            subgroup: Rc::new(p7()),
            a0: 3,
            a1: 4
        };
        x.pow(48);
        assert!(x.is_one());
    }

    #[test]
    fn powers_up_big() {
        let fp = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            )),
            Factorization::new(
                vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
            )
        ));
        let mut x = QuadNumber {
            subgroup: Rc::clone(&fp),
            a0: 3,
            a1: 5
        };
        x.pow(fp.pminusone.value());
        x.pow(fp.pplusone.value());
        assert!(x.is_one());
    }

    #[test]
    fn sylow_finds_generators() {
        let pplusone = Factorization::new(
            vec![(2, 1), (3, 2)]
        );
        let f = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 4)]
            )),
            pplusone.clone()
        ));
        let g = Rc::new(SylowDecomp::new(&f));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().factor(i);
            test_is_generator_small::<QuadFieldExt>(gen, d);
        }
    }

    #[test]
    fn sylow_finds_generators_big() {
        let pplusone = Factorization::new(
            vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
        );
        let fp = Rc::new(QuadFieldExt::make(
            &Rc::new(Factorization::new(
                vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            )),
            pplusone.clone()
        ));
        let g = Rc::new(SylowDecomp::new(&fp));
        for i in 0..g.generators.len() {
            let gen = &g.generators[i];
            let d = g.size().prime_powers()[i];
            test_is_generator_big::<QuadFieldExt>(gen, d);
        }
    }
}

