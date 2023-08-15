use std::collections::{HashSet, HashMap};

use itertools::*;

use crate::markoff::disjoint::Disjoint;
use crate::numbers::fp::*;

pub struct OrbitTester<'a, const P: u128> {
    f: &'a FpStar<P>,
    targets: HashSet<u128>
}

pub struct OrbitTesterResults {
    failures: u128,
    results: HashMap<u128, Disjoint<u128, bool>>
}

impl<'a, const P: u128> OrbitTester<'a, P> {
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new(true, |x,y| *x && *y));
        }

        let mut inv2 = self.f.from_int(2);
        inv2 = inv2.invert(self.f);

        let mut failures = 0;

        for (x,y) in self.targets.iter()
            .combinations_with_replacement(2)
            .map(|v| (v[0], v[1]))
        {
            let x = self.f.from_int(*x);
            let y = self.f.from_int(*y);
            let disc = 9 * x * y - 4 * (x * x + y * y);
            let z = 3 * x * y;
            let mut candidates = Vec::new();

            match disc.int_sqrt() {
                Some(FpNum{value: 0}) => {
                    candidates.push(z * inv2);
                },
                Some(x) => {
                    candidates.push((z + x) * inv2);
                    candidates.push((z - x) * inv2);
                },
                None => {
                    failures += 1;
                }
            }

            let it: Vec<(&FpNum<P>, bool)> = candidates.iter()
                .map(|z| {
                    (z, results.contains_key(&z.value))
                }).collect();
            let Some(disjoint) = results.get_mut(&x.value) else { continue; };
            for (z, pred) in it {
                if pred {
                    disjoint.associate(&x.value, &z.value);
                    disjoint.associate(&y.value, &z.value);
                } else {
                    disjoint.update(&x.value, false);
                    disjoint.update(&y.value, false);
                }
            }
        }

        OrbitTesterResults { 
            failures,
            results 
        }
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

    pub fn failures(&self) -> u128 {
        self.failures
    }
}
