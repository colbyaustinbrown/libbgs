//! Various number theory utility methods used throughout the libbgs crate.
use std::arch::asm;

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
pub const fn intpow_const<const M: u128>(mut x: u128, mut n: u128) -> u128 {
    if n == 0 {
        return 1;
    }
    let mut y = 1;
    while n > 1 {
        if n % 2 == 1 {
            y = if M == 0 {
                y * x
            } else {
                long_multiply_const::<M>(x, y)
            };
        }
        x = if M == 0 {
            x * x
        } else {
            long_multiply_const::<M>(x, x)
        };
        n >>= 1;
    }
    if M == 0 {
        y * x
    } else {
        long_multiply_const::<M>(y, x)
    }
}

/// Returns `x` to the power of `n`, modulo `m`.
pub fn intpow<const M: u128>(mut x: u128, mut n: u128) -> u128 {
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
pub fn long_multiply<const M: u128>(mut a: u128, mut b: u128) -> u128 {
    a += 0;
    b += 0;
    let res_lo: u64;
    let res_hi: u64;
    unsafe {
        asm!(
            "mov {res_lo}, 0",
            "mov {res_hi}, 0",
            "mov {tmp}, 0",
            "mov {tmp2}, 0",
            "mov {msk}, 0",
            "2:",
            "mov {tmp}, {b_lo}",
            "or {tmp}, {b_hi}",
            "jz 2f",

            // b /= 2 (only the low bits for now)
            "shr {b_lo}, 1",

            // if the removed bit of b was 1,
            // add a to the res
            "mov {msk}, 0",
            "setnc {msk:l}",
            "sub {msk}, 1",
            "mov {tmp}, {msk}",
            "and {tmp}, {a_lo}",
            "and {msk}, {a_hi}",
            "add {res_lo}, {tmp}",
            "adc {res_hi}, {msk}",

            // b /= 2 (high bits now)
            // and a *= 2
            "mov {tmp}, 0",
            "shr {b_hi}, 1",
            "setc {tmp:l}",
            "shl {tmp}, 63",
            "xor {b_lo}, {tmp}",
            "shl {a_hi}, 1",
            "shl {a_lo}, 1",
            "adc {a_hi}, 0",

            // if a >= m,
            "mov {tmp}, 0",
            "mov {tmp2}, 0",
            "mov {msk}, 0",
            "cmp {a_lo}, {m_lo}",
            "setb {msk:l}",
            "cmp {a_hi}, {m_hi}",
            "setne {tmp:l}",
            "setbe {tmp2:l}",
            "or {msk}, {tmp}",
            "and {msk}, {tmp2}",
            "sub {msk}, 1",
            // a -= m
            "mov {tmp}, {msk}",
            "and {tmp}, {m_lo}",
            "and {msk}, {m_hi}",
            "sub {a_lo}, {tmp}",
            "adc {msk}, 0",
            "sub {a_hi}, {msk}",

            // if res >= m,
            "mov {tmp}, 0",
            "mov {tmp2}, 0",
            "cmp {res_lo}, {m_lo}",
            "setb {msk:l}",
            "cmp {res_hi}, {m_hi}",
            "setne {tmp:l}",
            "setbe {tmp2:l}",
            "or {msk}, {tmp}",
            "and {msk}, {tmp2}",
            "sub {msk}, 1",
            // res -= m
            "mov {tmp}, {msk}",
            "and {tmp}, {m_lo}",
            "and {msk}, {m_hi}",
            "sub {res_lo}, {tmp}",
            "adc {msk}, 0",
            "sub {res_hi}, {msk}",

            "jmp 2b",
            "2:",
        
            a_hi = inout(reg) (a >> 64) as u64 => _,
            a_lo = inout(reg) (a & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            b_hi = inout(reg) (b >> 64) as u64 => _,
            b_lo = inout(reg) (b & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            m_lo = inout(reg) (M & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            m_hi = in(reg) (M >> 64) as u64,
            res_hi = out(reg) res_hi,
            res_lo = out(reg) res_lo,
            tmp = out(reg) _,
            tmp2 = out(reg) _,
            msk = out(reg) _,
        );
    }
    ((res_hi as u128) << 64) | (res_lo as u128)
}

/// Returns the product of `a` and `b` modulo `m`.
/// This function will panic if `m >= 2^127`.
/// Otherwise, it is guarenteed that there will not be integer overflow.
pub const fn long_multiply_const<const M: u128>(mut a: u128, mut b: u128) -> u128 {
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

/// Returns the Legendre symbol of `a` modulo `p`, i.e.,
/// $$\left(\frac{a}{p}\right)_L = a^{\frac{p - 1}{2}} \mod p$$.
pub fn legendre<const P: u128>(a: u128) -> u128 {
    intpow::<P>(a, (P - 1) / 2)
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
            if intpow::<P>(a, (P - 1) / 2) == P - 1 {
                res = a;
                break;
            }
        }
        res
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    const BIG_P: u128 = 1_000_000_000_000_000_124_399;

    #[test]
    fn test_long_multiply() {
        let a = 109_9511_627_777 % BIG_P;
        let b = 846_719_626_338_931_482_199_954 % BIG_P;
        let res = long_multiply_const::<BIG_P>(a, b);
        assert_eq!(res, 641287801218053764509);
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 641287801218053764509);
    }

    #[test]
    fn test_long_multiply_2() {
        let a = 109_9511_627_777 % BIG_P;
        let b = 590_295_810_358_705_651_713 % BIG_P;
        let res = long_multiply_const::<BIG_P>(a, b);
        assert_eq!(res, 443668638203145304518);
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 443668638203145304518);
    }

    #[test]
    fn test_long_multiply_3() {
        let a = 590_295_810_358_705_651_713 % BIG_P;
        let b = 1_267_650_600_818_525_211_855_408_857_088 % BIG_P;
        let res = long_multiply_const::<BIG_P>(a, b);
        assert_eq!(res, 917612211409095050983);
        let res = long_multiply::<BIG_P>(a, b);
        assert_eq!(res, 917612211409095050983);
    }

    #[test]
    fn test_long_multiply_4() {
        let a = (1 << 63) % BIG_P;
        let res = long_multiply_const::<BIG_P>(a, a);
        assert_eq!(res, 283_147_111_208_487_425_469);
        let res = long_multiply::<BIG_P>(a, a);
        assert_eq!(res, 283_147_111_208_487_425_469);


        let a = ((1 << 127) - 1) % BIG_P;
        let res = long_multiply_const::<BIG_P>(a, a);
        assert_eq!(res, 327_682_438_216_164_803_859);
        let res = long_multiply::<BIG_P>(a, a);
        assert_eq!(res, 327_682_438_216_164_803_859);
    }
}
