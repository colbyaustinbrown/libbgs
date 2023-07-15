use std::rc::Rc;

use crate::factorization::*;
use crate::util::*;
use crate::semigroup::*;

pub struct QuadFieldExt {
    p: u128,
    pub pplusone: Factorization,
    pub pminusone: Factorization,
    r: u128
}

type QuadSubgroup = QuadFieldExt;

impl Semigroup for QuadSubgroup {
    fn order(&self) -> u128 {
        self.pplusone.value
    }
}

impl QuadFieldExt {
    fn new(pplusone: Factorization, pminusone: Factorization) {
        let p = pplusone.value * pminusone.value;
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
        };
    }
}

#[derive(Clone)]
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
    fn one(param: &Rc<QuadSubgroup>) -> QuadNumber {
        QuadNumber {
            subgroup: Rc::clone(param),
            a0: 1,
            a1: 0
        }
    }

    fn is_one(&self) -> bool {
        self.a0 == 1 && self.a0 == 0
    }

    fn param(&self) -> &Rc<QuadSubgroup> {
        &self.subgroup
    }

    fn multiply(&mut self, other: &QuadNumber) {
        let f = &self.subgroup;
        self.a0 = (self.a0 * other.a0 + self.a1 * other.a1 * f.r) % f.p;
        self.a1 = (self.a1 * other.a0 + self.a0 * other.a1) % f.p;
    }

    fn square(&mut self) {
            let f = &self.subgroup;
            let a0_old = self.a0;
            self.a0 = (self.a0 * self.a0 + self.a1 * self.a1 * f.r) % f.p;
            self.a1 = (self.a1 * a0_old + a0_old * self.a1) % f.p;
    }
}

type Fp = u128;
impl Semigroup for Fp {
    fn order(&self) -> u128 {
        *self
    }
}
type FpNumber<'a> = (u128, Rc<u128>);

impl<'a> SemigroupElem for FpNumber<'a> {
    type Group = Fp;
    fn one(param: &Rc<Fp>) -> FpNumber<'a> {
        (1, Rc::clone(param))
    }
    fn is_one(&self) -> bool {
        self.0 == 1
    }

    fn param(&self) -> &Rc<Fp> {
        &self.1
    }

    fn multiply(&mut self, other: &FpNumber) {
        self.0 = self.0 * other.0 % (*self.1);
    }

    fn square(&mut self) {
        self.0 = self.0 * self.0 % (*self.1);
    }
}

