use crate::factorization::*;
use crate::util::*;
use crate::semigroup::*;

pub struct QuadField {
    p: u128,
    pub pplusone: Factorization,
    pub pminusone: Factorization,
    r: u128
}

impl QuadField {
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
        QuadField {
            p: p,
            pplusone,
            pminusone,
            r: r
        };
    }
}

pub struct QuadNumber<'a> {
    field: &'a QuadField,
    a0: u128,
    a1: u128
}

impl<'a> QuadNumber<'a> {
    fn from_int(field: &QuadField, a0: u128) -> QuadNumber {
        QuadNumber {
            field,
            a0,
            a1: 0
        }
    }
}

impl<'a> ParameterizedSemigroup<&'a QuadField> for QuadNumber<'a> {
    fn one(param: &'a QuadField) -> QuadNumber<'a> {
        QuadNumber {
            field: param,
            a0: 1,
            a1: 0
        }
    }

    fn param(&self) -> &'a QuadField {
        self.field
    }

    fn multiply(&mut self, other: &QuadNumber) {
        let f = self.field;
        self.a0 = (self.a0 * other.a0 + self.a1 * other.a1 * f.r) % f.p;
        self.a1 = (self.a1 * other.a0 + self.a0 * other.a1) % f.p;
    }

    fn square(&mut self) {
            let f = self.field;
            let a0_old = self.a0;
            self.a0 = (self.a0 * self.a0 + self.a1 * self.a1 * f.r) % f.p;
            self.a1 = (self.a1 * a0_old + a0_old * self.a1) % f.p;
    }
}

type FpNumber<'a> = (u128, &'a u128);

impl<'a> ParameterizedSemigroup<&'a u128> for FpNumber<'a> {
    fn one(param: &'a u128) -> FpNumber<'a> {
        (1, param)
    }

    fn param(&self) -> &'a u128 {
        self.1
    }

    fn multiply(&mut self, other: &FpNumber) {
        self.0 = self.0 * other.0 % self.1;
    }

    fn square(&mut self) {
        self.0 = self.0 * self.0 % self.1;
    }
}

