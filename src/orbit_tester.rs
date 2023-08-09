use std::collections::{HashSet, HashMap};

use itertools::*;

use crate::markoff::disjoint::Disjoint;
use crate::numbers::fp::*;
use crate::util::*;

pub struct OrbitTester<'a> {
    f: &'a FpStar,
    targets: HashSet<u128>
}

pub struct OrbitTesterResults {
    results: HashMap<u128, Disjoint<u128, bool>>
}

impl<'a> OrbitTester<'a> {
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new(true, |x,y| *x && *y));
        }

        let p = self.f.p();
        let mut inv2 = self.f.from_int(2);
        inv2.invert(self.f);

        for (x,y) in self.targets.iter()
            .cartesian_product(self.targets.iter())
        {
            let x = *x;
            let y = *y;
            let x2 = intpow(x, 2, p);
            let y2 = intpow(y, 2, p);
            let tmp = self.f.from_int(long_multiply(4, x2 + y2, p));
            let mut disc = self.f.from_int(long_multiply(9, long_multiply(x, y, p), p));
            disc.sub(&tmp, self.f);
            let mut z = self.f.from_int(long_multiply(3, long_multiply(x, y, p), p));
            let mut candidates = Vec::new();
            match disc.int_sqrt(self.f) {
                Some(FpNum{value: 0}) => {
                    z.multiply(&inv2, self.f);
                    candidates.push(z);
                },
                Some(x) => {
                    let mut z1 = z.clone();
                    z1.add(&x, self.f);
                    z1.multiply(&inv2, self.f);
                    candidates.push(z1);

                    z.sub(&x, self.f);
                    z.multiply(&inv2, self.f);
                    candidates.push(z);
                },
                None => {}
            }

            let it: Vec<(&FpNum, bool)> = candidates.iter()
                .map(|z| {
                    (z, results.contains_key(&z.value()))
                }).collect();
            let Some(disjoint) = results.get_mut(&x) else { continue; };
            for (z, pred) in it {
                if pred {
                    disjoint.associate(&y, &z.value());
                } else {
                    disjoint.update(&y, false);
                }
            }
        }

        OrbitTesterResults { results }
    }

    pub fn new(f: &FpStar) -> OrbitTester {
        OrbitTester {
            f,
            targets: HashSet::new()
        }
    }

    pub fn add_target(mut self, t: u128) -> OrbitTester<'a> {
        self.targets.insert(t);
        self
    }
}

impl OrbitTesterResults {
    pub fn results(&self) -> impl Iterator<Item = (&u128, &Disjoint<u128, bool>)> {
        self.results.iter()
    }
}
