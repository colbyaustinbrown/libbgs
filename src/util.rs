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
pub const fn intpow(mut x: u128, mut n: u128, m: u128) -> u128 {
    if n == 0 {
        return 1;
    }
    let mut y = 1;
    while n > 1 {
        if n % 2 == 1 {
            y = if m == 0 {
                y * x
            } else {
                long_multiply(x, y, m)
            };
        }
        x = if m == 0 {
            x * x
        } else {
            long_multiply(x, x, m)
        };
        n >>= 1;
    }
    if m == 0 {
        y * x
    } else {
        long_multiply(y, x, m)
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
pub const fn long_multiply(mut a: u128, mut b: u128, m: u128) -> u128 {
    a %= m;
    b %= m;

    if a < (1 << 63) && b < (1 << 63) {
        return (a * b) % m;
    } else if m >= (1 << 126) {
        panic!("p too large! Greater than 1 << 126.");
    }

    let mut res = 0;
    while b > 0 {
        if b & 1 == 1 {
            res += a;
            // Note: this is significantly faster (~40%)
            // than res %= m on benchmarking
            if res >= m {
                res -= m;
            }
        }
        a *= 2;
        // see above comment
        if a >= m {
            a -= m;
        }
        b /= 2;
    }
    res
}

/// Returns the Legendre symbol of `a` modulo `p`, i.e.,
/// $$\left(\frac{a}{p}\right)_L = a^{\frac{p - 1}{2}} \mod p$$.
pub fn legendre(a: u128, p: u128) -> u128 {
    intpow(a, (p - 1) / 2, p)
}

/// Returns a quadratic non-residue modulo `p`.
/// That is, it returns an integer $a \in \mathbb{Z} / p\mathbb{Z}$ such that there is no $x$
/// satisfying $x^2 = a \mod p$.
pub fn find_nonresidue(p: u128) -> u128 {
    if p % 4 == 3 {
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
    }
}
