fn gcd(mut a: u128, mut b: u128) -> u128 {
    let mut t;
    while b != 0 {
        t = a % b;
        a = b;
        b = t;
    }
    a
}

pub fn intpow(mut x: u128, mut n: u128, m: u128) -> u128 {
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

pub fn standard_affine_shift(q: u128, i: u128) -> u128 {
    let mut m = 4 * q / 5;
    while gcd(m, q) != 1 {
        m -= 1;
    }
    let a = 2 * q / 3;
    (m * i + a) % q
}

pub fn long_multiply(mut a: u128, mut b: u128, m: u128) -> u128 {
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

pub fn legendre(a: u128, p: u128) -> u128 {
    intpow(a, (p - 1) / 2, p)
}

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
