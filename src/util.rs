//! Various number theory utility methods used throughout the libbgs crate.
use std::arch::{asm, global_asm};

/*
 * Registers used:
 * rdi, rsi -- a lo and hi (input, clobbered)
 * rcx, rdx -- b lo and hi (input, clobbered)
 * r8, r9 -- m lo and hi (input, clobbered)
 * r10, r11 -- output registers
 * rax, r13 -- temporary values
 */
global_asm!(
    ".globl long_multiply",
    "long_multiply:",
    "xor r10, r10",
    "xor r11, r11",

    "2:",
    // if b % 2 == 1
    // b /= 2
    "xor rax, rax",
    "shrd rdi, rsi, 1",
    "setnc al",
    "shr rsi, 1",
    "sub rax, 1",
    
    // res += a
    "mov r13, rax",
    "and r13, rcx",
    "and rax, rdx",
    "add r10, r13",
    "adc r11, rax",

    // if res >= m,
    "mov rax, r8",
    "or rax, r9",
    "jz 5f",
    "xor rax, rax",
    "cmp r11, r9",
    "jne 5f",
    "cmp r10, r8",
    "5:",
    "setbe al",
    "sub rax, 1",
    // res -= m
    "mov r13, rax",
    "and r13, r8",
    "and rax, r9",
    "sub r10, r13",
    "sbb r11, rax",

    // a *= 2
    "shld rdx, rcx, 1",
    "shl rcx, 1",

    // if a >= m,
    "mov rax, r8",
    "or rax, r9",
    "jz 5f",
    "xor rax, rax",
    "cmp rdx, r9",
    "jne 5f",
    "cmp rcx, r8",
    "5:",
    "setbe al",
    "sub rax, 1",
    // a -= m
    "mov r13, rax",
    "and r13, r8",
    "and rax, r9",
    "sub rcx, r13",
    "sbb rdx, rax",

    "mov r13, rdi",
    "or r13, rsi",
    "jnz 2b",
    "ret",
);

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
pub fn intpow<const M: u128>(x: u128, n: u128) -> u128 {
    let res_lo: u64;
    let res_hi: u64;
    unsafe {
        asm!(
            "mov r15, 0",       // y_hi = 0
            "mov r14, 1",       // y_lo = 1 
            "mov rcx, rdi",     // b_lo = x_lo
            "mov rdx, rsi",     // b_hi = y_hi
            "xor r10, r10",     // r_lo = 0
            "xor r11, r11",     // r_hi = 0
            "cmp {n_hi}, 0",
            "jne 2f",
            "cmp {n_lo}, 1",
            "ja 2f",
            "jb 9f",
            "mov r10, rdi",
            "mov r11, rsi",
            "jmp 4f",
            "9:",
            "mov r11, 0",
            "mov r10, 1",
            "jmp 4f",

            "2:",
            "shrd {n_lo}, {n_hi}, 1", // n /= 2
            "setc al",
            "shr {n_hi}, 1",
            "test rax, 1",          // if n was odd
            "jz 3f",                // {

            "xchg rcx, r14",
            "xchg rdx, r15",
            "call long_multiply",   //      y = long_multiply(x, y, m)
            "mov rdi, r14",
            "mov rsi, r15",         
            "mov r14, r10",
            "mov r15, r11",

            "3:",                   // }
            "mov rcx, rdi",         // 
            "mov rdx, rsi",
            "call long_multiply",
            "mov rdi, r10",
            "mov rsi, r11",
            "mov rcx, rdi",
            "mov rdx, rsi",

            "cmp {n_hi}, 0",
            "jne 2b",
            "cmp {n_lo}, 1",
            "jne 2b",

            "mov rcx, r14",
            "mov rdx, r15",
            "call long_multiply",
            "4:",

            in("r8") (M & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64,
            in("r9") (M >> 64) as u64,
            inout("rdi") (x & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            inout("rsi") (x >> 64) as u64 => _,
            out("rcx") _,
            out("rdx") _,
            n_lo = inout(reg) (n & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            n_hi = inout(reg) (n >> 64) as u64 => _,
            out("rax") _,
            out("r13") _,
            out("r10") res_lo,
            out("r11") res_hi,
            out("r14") _,
            out("r15") _,
            options(pure, nomem, nostack),
        );
    }
    ((res_hi as u128) << 64) | (res_lo as u128)
    /*
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
    */
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
pub fn long_multiply<const M: u128>(a: u128, b: u128) -> u128 {
    let res_lo: u64;
    let res_hi: u64;
    unsafe {
        asm!(
            "xor r10, r10",
            "xor r11, r11",

            "2:",
            // if b % 2 == 1
            // b /= 2
            "xor rax, rax",
            "shrd rdi, rsi, 1",
            "setnc al",
            "shr rsi, 1",
            "sub rax, 1",
            
            // res += a
            "mov r13, rax",
            "and r13, rcx",
            "and rax, rdx",
            "add r10, r13",
            "adc r11, rax",

            // a *= 2
            "shld rdx, rcx, 1",
            "shl rcx, 1",

            // "mov rax, r9",
            // "or rax, r8",
            // "jz 3f",

            // if res >= m,
            "xor rax, rax",
            "cmp r11, r9",
            "jne 5f",
            "cmp r10, r8",
            "5:",
            "setbe al",
            "sub rax, 1",
            // res -= m
            "mov r13, rax",
            "and r13, r8",
            "and rax, r9",
            "sub r10, r13",
            "sbb r11, rax",

            // if a >= m,
            "xor rax, rax",
            "cmp rdx, r9",
            "jne 5f",
            "cmp rcx, r8",
            "5:",
            "setbe al",
            "sub rax, 1",
            // a -= m
            "mov r13, rax",
            "and r13, r8",
            "and rax, r9",
            "sub rcx, r13",
            "sbb rdx, rax",

            "3:",
            "mov r13, rdi",
            "or r13, rsi",
            "jnz 2b",
            inout("rcx") (a & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            inout("rdx") (a >> 64) as u64 => _,
            inout("rdi") (b & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64 => _,
            inout("rsi") (b >> 64) as u64 => _,
            in("r8") (M & 0xFF_FF_FF_FF_FF_FF_FF_FF) as u64,
            in("r9") (M >> 64) as u64,
            out("r10") res_lo,
            out("r11") res_hi,
            out("rax") _,
            out("r13") _,
            options(pure, nomem, nostack),
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
    fn test_long_multiply_1() {
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
    fn test_intpow_1() {
        assert_eq!(intpow::<13>(10, 5), 4);
    }

    #[test]
    fn test_intpow_2() {
        let a: u128 = 123_456_789_123_456_789;
        let b: u128 = 123_456;
        let res = intpow::<BIG_P>(a, b);
        assert_eq!(res, 748_418_297_532_596_964_525);
    }

    #[test]
    fn test_intpow_3() {
        assert_eq!(intpow::<0>(3, 1), 3);
    }

    #[test]
    fn test_intpow_4() {
        assert_eq!(intpow::<0>(999, 0), 1);
    }
}
