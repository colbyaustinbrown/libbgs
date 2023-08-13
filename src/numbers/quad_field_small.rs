use crate::util::*;
use crate::numbers::semigroup::*;

#[derive(PartialEq, Eq, Debug)]
pub struct QuadFieldSml {
    p: u128,
    r: u128
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct QuadNumSml {
    pub a0: u128,
    pub a1: u128
}

impl QuadFieldSml {
    pub fn new(p: u128, r: u128) -> QuadFieldSml {
        QuadFieldSml {
            p,
            r
        }
    }

    fn r(&self) -> u128 {
        self.r
    }

    fn change_r(&self, r: u128) -> QuadFieldSml {
        QuadFieldSml {
            p: self.p,
            r
        }
    }

    fn p(&self) -> u128 {
        self.p
    }
}

impl QuadNumSml {
    pub fn from_ints(a0: u128, a1: u128) -> Self {
        QuadNumSml {
            a0,
            a1
        }
    }
}

impl Semigroup for QuadFieldSml {
    type Elem = QuadNumSml;

    fn size(&self) -> u128 {
        self.p - 1
    }
    fn one(&self) -> QuadNumSml {
        QuadNumSml {
            a0: 1,
            a1: 0
        }
    }
}

impl SemigroupElem for QuadNumSml {
    type Group = QuadFieldSml;

    fn is_one(&self) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&mut self, other: &QuadNumSml, f: &QuadFieldSml) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, other.a0, p) + long_multiply(a0_old, other.a1, p)) % p;
    }

    fn square(&mut self, f: &QuadFieldSml) {
        let a0_old = self.a0;
        let p = f.p();
        self.a0 = (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r(), p), p)) % p;
        self.a1 = (long_multiply(self.a1, a0_old, p) + long_multiply(a0_old, self.a1, p)) % p;
    }
}

