use crate::util::*;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Factorization {
    pub value: u128,
    pub prime_powers: Vec<(u128, u128)>,
}

impl Factorization {
    pub fn factor(&self, i: usize) -> u128 {
        intpow(self.prime_powers[i].0, self.prime_powers[i].1, 0)
    }

    pub fn new(mut n: u128) {
        // TODO replace with an actual factorization algorithm
        let mut res = Factorization {
            value: n,
            prime_powers: Vec::new()
        };
        let mut i = 2;
        while n > 1 {
            let mut pow = 0;
            while n % i == 0 {
                pow += 1;
                n /= i;
            }
            if pow > 0 {
                res.prime_powers.push((n, pow as u128));
            }
            i += 1;
        }
    }

    pub fn len(&self) -> usize {
        self.prime_powers.len()
    }
}

