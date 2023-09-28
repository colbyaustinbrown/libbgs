//! Various number theory utility methods used throughout the libbgs crate.

const fn gcd(mut a: u128, mut b: u128) -> u128 {
    let mut t;
    while b != 0 {
        t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Returns `x` to the power of `n`, modulo `m`.
pub const fn intpow<const M: u128>(mut x: u128, mut n: u128) -> u128 {
    if n == 0 {
        return 1;
    }
    let mut y = 1;
    while n > 1 {
        if n % 2 == 1 {
            y = if M == 0 {
                y * x
            } else {
                long_multiply::<M>(x, y)
            };
        }
        x = if M == 0 {
            x * x
        } else {
            long_multiply::<M>(x, x)
        };
        n >>= 1;
    }
    if M == 0 {
        y * x
    } else {
        long_multiply::<M>(y, x)
    }
}


/// Returns a pseudo-random integer modulo `q`, unique for every `i` between `0` and `q`.
/// This acts suitably well as a random number generator for several modular arithmetic operations,
/// including randomly searching for quadratic (non) residues.
pub const fn standard_affine_shift(q: u128, i: u128) -> u128 {
    let mut m = 4 * q / 5;
    while gcd(m, q) != 1 {
        m -= 1;
    }
    let a = 2 * q / 3;
    (m * i + a) % q
}

/// Returns the product of `a` and `b` modulo `m`.
/// This function will panic if `m >= 2^127`.
/// Otherwise, it is guarenteed that there will not be integer overflow.
pub const fn long_multiply<const M: u128>(mut a: u128, mut b: u128) -> u128 {
    if M == 0 {
        return a * b;
    }

    a %= M;
    b %= M;

    let mut res = 0;
    while b > 0 {
        if b & 1 == 1 {
            res += a;
            // Note: this is significantly faster (~40%)
            // than res %= m on benchmarking
            if res >= M {
                res -= M;
            }
        }
        a *= 2;
        // see above comment
        if a >= M {
            a -= M;
        }
        b /= 2;
    }
    res
}

/// Returns a quadratic non-residue modulo `p`.
/// That is, it returns an integer $a \in \mathbb{Z} / p\mathbb{Z}$ such that there is no $x$
/// satisfying $x^2 = a \mod p$.
pub fn find_nonresidue<const P: u128>() -> u128 {
    if P % 4 == 3 {
        P - 1
    } else if P % 8 == 3 || P % 8 == 5 {
        2
    } else {
        let mut res = 0;
        for i in 0..P {
            let a = standard_affine_shift(P, i);
            let half_pow = intpow::<P>(a % P, (P - 1) / 2);
            if half_pow == P - 1 {
                res = a;
                break;
            }
        }
        res
    }
}

pub const fn carrying_mul(a: u128, b: u128) -> (u128, u128) {
    let a_lo = a & 0xFF_FF_FF_FF_FF_FF_FF_FF;
    let a_hi = a >> 64;
    let b_lo = b & 0xFF_FF_FF_FF_FF_FF_FF_FF;
    let b_hi = b >> 64;
    
    let cross = a_hi * b_lo + a_lo * b_hi;
    let (c, res_lo) = carrying_add(a_lo * b_lo, (cross & 0xFF_FF_FF_FF_FF_FF_FF_FF) << 64);
    let res_hi = a_hi * b_hi + (cross >> 64) + c;
    (res_hi, res_lo)
}

pub const fn carrying_add(a: u128, b: u128) -> (u128, u128) {
    let lo = (a & 0xFF_FF_FF_FF_FF_FF_FF_FF) + (b & 0xFF_FF_FF_FF_FF_FF_FF_FF);
    let hi = (a >> 64) + (b >> 64) + (lo >> 64);
    (hi >> 64, (hi << 64) | (lo & 0xFF_FF_FF_FF_FF_FF_FF_FF))
}

pub const fn shrd(dst: u128, src: u128, n: usize) -> u128 {
    (dst >> n) | ((src & ((1 << n) - 1)) << (128 - n))
}

#[cfg(test)]
pub mod tests {
    use super::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    fn test_long_multiply() {
        let a = 109_9511_627_777 % BIG_P;
        let b = 846_719_626_338_931_482_199_954 % BIG_P;
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 641287801218053764509);
    }

    #[test]
    fn test_long_multiply_2() {
        let a = 109_9511_627_777 % BIG_P;
        let b = 590_295_810_358_705_651_713 % BIG_P;
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 443668638203145304518);
    }

    #[test]
    fn test_long_multiply_3() {
        let a = 590_295_810_358_705_651_713 % BIG_P;
        let b = 1_267_650_600_818_525_211_855_408_857_088 % BIG_P;
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 917612211409095050983);
    }

    #[test]
    fn test_long_multiply_4() {
        let a = (1 << 63) % BIG_P;
        let res = long_multiply::<BIG_P>(a, a);
        assert_eq!(res, 283_147_111_208_487_425_469);


        let a = ((1 << 127) - 1) % BIG_P;
        let res = long_multiply::<BIG_P>(a, a);
        assert_eq!(res, 327_682_438_216_164_803_859);
    }

    #[test]
    fn test_long_multiply_5() {
        let a: u128 = 538_744_077_496_950_347_511;
        let b: u128 = 10_022_347_072_413_323_143;
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 1);
    }

    #[test]
    fn test_long_multiply_6() {
        assert_eq!(long_multiply::<0>(100, 100), 10_000);
    }

    #[test]
    fn carrying_multiply() {
        let (hi, lo) = carrying_mul(BIG_P, BIG_P);
        assert_eq!(hi, 0xb7a);
        assert_eq!(lo, 0xbc6270503128ac1b3edc190389237521);
    }

    #[test]
    fn carrying_addition() {
        let (hi, lo) = carrying_add(BIG_P << 58, BIG_P << 58);
        assert_eq!(hi, 1);
        assert_eq!(lo, 0xb1ae4d6e2ef50f2f7800000000000000);
    }
}
