use std::rc::Rc;

use crate::factorization::*;
use crate::util::*;
use crate::semigroup::*;

#[derive(Debug)]
pub struct QuadFieldExt {
    p: u128,
    pub pplusone: Factorization,
    pub pminusone: Factorization,
    r: u128
}

type QuadSubgroup = QuadFieldExt;

impl Semigroup for QuadSubgroup {
    type Elem = QuadNumber;
    fn order(&self) -> u128 {
        self.pplusone.value
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
    pub fn new(pminusone: Factorization, pplusone: Factorization) -> QuadFieldExt {
        let p = pplusone.value - 1;
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
        QuadFieldExt {
            p: p,
            pplusone,
            pminusone,
            r: r
        }
    }
}

#[derive(Clone, Debug)]
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

    fn steinitz(subgroup: &Rc<QuadSubgroup>, i: u128) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(subgroup),
            a0: i % subgroup.p,
            a1: i / subgroup.p
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
        self.a0 = (long_multiply(self.a0, other.a0, f.p) + long_multiply(self.a1, long_multiply(other.a1, f.r, f.p), f.p)) % f.p;
        self.a1 = (long_multiply(self.a1, other.a0, f.p) + long_multiply(a0_old, other.a1, f.p)) % f.p;
    }

    fn square(&mut self) {
            let f = &self.subgroup;
            let a0_old = self.a0;
            self.a0 = (long_multiply(self.a0, self.a0, f.p) + long_multiply(self.a1, long_multiply(self.a1, f.r, f.p), f.p)) % f.p;
            self.a1 = (long_multiply(self.a1, a0_old, f.p) + long_multiply(a0_old, self.a1, f.p)) % f.p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p7() -> QuadFieldExt {
        QuadFieldExt::new(
            Factorization {
                value: 6, 
                factors: vec![2, 3], 
                prime_powers: vec![(2,1), (3,1)]
            },
            Factorization { 
                value: 8,
                factors: vec![8],
                prime_powers: vec![(2,3)]
            }
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
        for i in 2..f.p {
            assert_ne!((i * i) % f.p, f.r);
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
        let fp = Rc::new(QuadFieldExt::new(
            Factorization {
                value: 1_000_000_000_000_000_124_398,
                factors: vec![2, 7, 13, 841, 43, 705737, 215288719],
                prime_powers: vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
            },
            Factorization {
                value: 1_000_000_000_000_000_124_400,
                factors: vec![16, 3, 25, 121, 17, 19, 23, 97, 757, 1453, 8689],
                prime_powers: vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
            }
        ));
        let mut x = QuadNumber {
            subgroup: Rc::clone(&fp),
            a0: 3,
            a1: 5
        };
        x.pow(fp.pminusone.value);
        x.pow(fp.pplusone.value);
        assert!(x.is_one());
    }
}

