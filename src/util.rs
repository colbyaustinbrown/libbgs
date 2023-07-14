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
            y *= x;
            if m != 0 {
                y %= m;
            }
        }
        x *= x;
        if m != 0 {
            y %= m;
        }
        n >>= 1;
    }
    if m == 0 {y*x} else {y*x % m}
}

pub fn standard_affine_shift(q: u128, i: u128) -> u128 {
    let mut m = 4*q / 5;
    while gcd(m, q) != 1 {
        m = m - 1;
    }
    let a = 2*q / 3;
    (m * i + a) % q
}

