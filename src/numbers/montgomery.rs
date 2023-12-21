use std::ops::*;

use libbgs_util::*;

/// A number in the Montgomery modular arithmetic system, modulo `N`.
/// For more informtion, see: Montgomery, Peter (April 1985). "Modular Multiplication Without Trial
/// Division". Mathematics of Computation. 44 (170): 519-521.
///
/// # Example
/// ```
/// use libbgs::numbers::Montgomery;
/// let x = Montgomery::<7>::from(6);
/// let y = Montgomery::<7>::from(8);
/// let z = x * y;
/// assert_eq!(u128::from(z), 48 % 7);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Montgomery<const N: u128>(u128);

impl<const N: u128> Montgomery<N> {
    /// Returns the Montgomery representation of this number.
    #[inline(always)]
    pub fn raw(&self) -> u128 {
        self.0
    }

    // Elements are represented as aR + N in this system.
    // Chosen to be a power of two to turn multiplication / division into shift operations.
    const R: u128 = {
        let mut n = N.next_power_of_two();
        while n.ilog2() % 4 != 0 {
            n <<= 1;
        }
        n
    };

    // Mask used for taking elements modulo R. Since R is a power of 2,
    // x % R == x & MASK
    const MASK: u128 = Self::R - 1;

    // The length of R
    const SHIFT: u32 = Self::R.ilog2();

    // "Magic" number R' such that RR' is congruent to 1 mod N
    const MAGIC: u128 = {
        let mut r0 = Self::R;
        let mut r1 = N;
        let mut t0: i128 = 0;
        let mut t1: i128 = 1;
        let mut s0: i128 = 1;
        let mut s1: i128 = 0;
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (t0, t1) = (t1, t0 - (q as i128) * t1);
            (s0, s1) = (s1, s0 - (q as i128) * s1);
        }
        if t0 < 0 {
            (-t0) as u128
        } else {
            Self::R - (t0 as u128)
        }
    };

    const fn redc(x: u128) -> Montgomery<N> {
        let m = ((x & Self::MASK) * Self::MAGIC) & Self::MASK;
        let t = (x + m * N) >> Self::SHIFT;
        if t >= N {
            Montgomery(t - N)
        } else {
            Montgomery(t)
        }
    }

    const fn redc2((hi, lo): (u128, u128)) -> Montgomery<N> {
        let (_, m) = carrying_mul(lo & Self::MASK, Self::MAGIC);
        let m = m & Self::MASK;
        let (c1, t) = carrying_mul(m, N);
        let (c2, t) = carrying_add(t, lo);
        let t = shrd(t, c1 + c2 + hi, Self::SHIFT as usize);
        if t >= N {
            Montgomery(t - N)
        } else {
            Montgomery(t)
        }
    }

    // These next several methods are constant implementations of the usual std::ops methods.
    // Once const impls lands in Rust, they will be removed, and the impl std::ops will be made
    // const.
    pub(crate) const fn const_mul(&self, rhs: &Montgomery<N>) -> Montgomery<N> {
        let (hi, lo) = carrying_mul(self.0, rhs.0);
        Montgomery::<N>::redc2((hi, lo))
    }

    /// Converts a `u128` into its Montgomery representation.
    /// This operation is expensive.
    pub(crate) const fn from_u128(src: u128) -> Montgomery<N> {
        let r2 = long_multiply::<N>(Self::R, Self::R);
        Montgomery::<N>::redc2(carrying_mul(src, r2))
    }

    pub(crate) const fn const_pow(self, mut n: u128) -> Montgomery<N> {
        let mut x = self;
        let mut y = Montgomery::<N>::from_u128(1);
        if n == 0 {
            return y;
        }
        while n > 1 {
            if n % 2 == 1 {
                y = y.const_mul(&x);
            }
            x = x.const_mul(&x);
            n >>= 1;
        }
        x.const_mul(&y)
    }
}

impl<const N: u128> From<u128> for Montgomery<N> {
    fn from(src: u128) -> Montgomery<N> {
        Montgomery::<N>::from_u128(src)
    }
}

impl<const N: u128> From<Montgomery<N>> for u128 {
    fn from(src: Montgomery<N>) -> u128 {
        Montgomery::<N>::redc2((0, src.0)).0
    }
}

impl<const N: u128> Add<Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn add(self, rhs: Montgomery<N>) -> Montgomery<N> {
        let sum = self.0 + rhs.0;
        if sum >= N {
            Montgomery(sum - N)
        } else {
            Montgomery(sum)
        }
    }
}

impl<const N: u128> Add<Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn add(self, rhs: Montgomery<N>) -> Montgomery<N> {
        let sum = self.0 + rhs.0;
        if sum >= N {
            Montgomery(sum - N)
        } else {
            Montgomery(sum)
        }
    }
}

impl<const N: u128> Add<&Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn add(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        let sum = self.0 + rhs.0;
        if sum >= N {
            Montgomery(sum - N)
        } else {
            Montgomery(sum)
        }
    }
}

impl<const N: u128> Add<&Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn add(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        let sum = self.0 + rhs.0;
        if sum >= N {
            Montgomery(sum - N)
        } else {
            Montgomery(sum)
        }
    }
}

impl<const N: u128> Sub<Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn sub(self, rhs: Montgomery<N>) -> Montgomery<N> {
        if self.0 >= rhs.0 {
            Montgomery(self.0 - rhs.0)
        } else {
            Montgomery(N + self.0 - rhs.0)
        }
    }
}

impl<const N: u128> Sub<Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn sub(self, rhs: Montgomery<N>) -> Montgomery<N> {
        if self.0 >= rhs.0 {
            Montgomery(self.0 - rhs.0)
        } else {
            Montgomery(N + self.0 - rhs.0)
        }
    }
}

impl<const N: u128> Sub<&Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn sub(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        if self.0 >= rhs.0 {
            Montgomery(self.0 - rhs.0)
        } else {
            Montgomery(N + self.0 - rhs.0)
        }
    }
}

impl<const N: u128> Sub<&Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn sub(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        if self.0 >= rhs.0 {
            Montgomery(self.0 - rhs.0)
        } else {
            Montgomery(N + self.0 - rhs.0)
        }
    }
}

impl<const N: u128> Mul<Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn mul(self, rhs: Montgomery<N>) -> Montgomery<N> {
        self.const_mul(&rhs)
    }
}

impl<const N: u128> Mul<Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn mul(self, rhs: Montgomery<N>) -> Montgomery<N> {
        self.const_mul(&rhs)
    }
}

impl<const N: u128> Mul<&Montgomery<N>> for Montgomery<N> {
    type Output = Montgomery<N>;
    fn mul(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        self.const_mul(rhs)
    }
}

impl<const N: u128> Mul<&Montgomery<N>> for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn mul(self, rhs: &Montgomery<N>) -> Montgomery<N> {
        self.const_mul(rhs)
    }
}

impl<const N: u128> Neg for Montgomery<N> {
    type Output = Montgomery<N>;
    fn neg(self) -> Montgomery<N> {
        if self.0 == 0 {
            Montgomery(0)
        } else {
            Montgomery(N - self.0)
        }
    }
}

impl<const N: u128> Neg for &Montgomery<N> {
    type Output = Montgomery<N>;
    fn neg(self) -> Montgomery<N> {
        if self.0 == 0 {
            Montgomery(0)
        } else {
            Montgomery(N - self.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    fn finds_magic() {
        assert_eq!(Montgomery::<13>::MAGIC, 11);
    }

    #[test]
    fn test_multiply() {
        let five = Montgomery::<13>::from(5);
        let six = Montgomery::<13>::from(6);
        let thirty = &five * &six;
        assert_eq!(u128::from(thirty), 30 % 13);
    }

    #[test]
    fn test_multiply_big() {
        let a = 0xABCDEF0123456789;
        let b = 0x9876543210ABCDEF;
        let am = Montgomery::<BIG_P>::from(a);
        let bm = Montgomery::<BIG_P>::from(b);
        let res = am * bm;
        assert_eq!(u128::from(res), (a * b) % BIG_P);
    }

    #[test]
    fn test_pow() {
        let a = 0xABCDEF;
        let n = 0xFEDCBA;
        let res = Montgomery::<BIG_P>::from(a).const_pow(n);
        assert_eq!(u128::from(res), intpow::<BIG_P>(a, n));
    }
}
