use std::collections::{HashSet, HashMap};

use itertools::*;

use crate::markoff::disjoint::Disjoint;
use crate::numbers::fp::*;
use crate::util::*;

pub struct OrbitTester<'a, const P: u128> {
    f: &'a FpStar<P>,
    targets: HashSet<u128>
}

pub struct OrbitTesterResults {
    results: HashMap<u128, Disjoint<u128, bool>>
}

impl<'a, const P: u128> OrbitTester<'a, P> {
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new(true, |x,y| *x && *y));
        }

        let mut inv2 = self.f.from_int(2);
        inv2.invert(self.f);

        for (x,y) in self.targets.iter()
            .cartesian_product(self.targets.iter())
        {
            let x = *x;
            let y = *y;
            let x2 = intpow(x, 2, P);
            let y2 = intpow(y, 2, P);
            let tmp = self.f.from_int(long_multiply(4, x2 + y2, P));
            let mut disc = self.f.from_int(long_multiply(9, long_multiply(x, y, P), P));
            disc.sub(&tmp);
            let mut z = self.f.from_int(long_multiply(3, long_multiply(x, y, P), P));
            let mut candidates = Vec::new();
            match disc.int_sqrt() {
                Some(FpNum{value: 0}) => {
                    z.multiply(&inv2, self.f);
                    candidates.push(z);
                },
                Some(x) => {
                    let mut z1 = z.clone();
                    z1.add(&x);
                    z1.multiply(&inv2, self.f);
                    candidates.push(z1);

                    z.sub(&x);
                    z.multiply(&inv2, self.f);
                    candidates.push(z);
                },
                None => {}
            }

            let it: Vec<(&FpNum<P>, bool)> = candidates.iter()
                .map(|z| {
                    (z, results.contains_key(&z.value))
                }).collect();
            let Some(disjoint) = results.get_mut(&x) else { continue; };
            for (z, pred) in it {
                if pred {
                    disjoint.associate(&y, &z.value);
                } else {
                    disjoint.update(&y, false);
                }
            }
        }

        OrbitTesterResults { results }
    }

    pub fn new(f: &FpStar<P>) -> OrbitTester<P> {
        OrbitTester {
            f,
            targets: HashSet::new()
        }
    }

    pub fn add_target(mut self, t: u128) -> OrbitTester<'a, P> {
        self.targets.insert(t);
        self
    }
}

impl OrbitTesterResults {
    pub fn results(&self) -> impl Iterator<Item = (&u128, &Disjoint<u128, bool>)> {
        self.results.iter()
    }
}
