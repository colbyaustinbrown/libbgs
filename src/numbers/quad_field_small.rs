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

    fn is_one(&self, _f: &QuadFieldSml) -> bool {
        self.a0 == 1 && self.a1 == 0
    }

    fn multiply(&self, other: &QuadNumSml, f: &QuadFieldSml) -> QuadNumSml {
        let p = f.p();
        QuadNumSml {
            a0: (long_multiply(self.a0, other.a0, p) + long_multiply(self.a1, long_multiply(other.a1, f.r(), p), p)) % p,
            a1: (long_multiply(self.a1, other.a0, p) + long_multiply(self.a0, other.a1, p)) % p
        }
    }

    fn square(&self, f: &QuadFieldSml) -> QuadNumSml {
        let p = f.p();
        QuadNumSml {
            a0: (long_multiply(self.a0, self.a0, p) + long_multiply(self.a1, long_multiply(self.a1, f.r(), p), p)) % p,
            a1: (long_multiply(self.a1, self.a0, p) + long_multiply(self.a0, self.a1, p)) % p
        }
    }
}

