use std::collections::{HashMap, HashSet};

use itertools::*;

use crate::markoff::disjoint::Disjoint;
use crate::numbers::fp::*;

pub struct OrbitTester<'a, const P: u128> {
    f: &'a FpStar<P>,
    targets: HashSet<u128>,
}

pub struct OrbitTesterResults {
    failures: u128,
    results: HashMap<u128, Disjoint<u128>>,
}

impl<'a, const P: u128> OrbitTester<'a, P> {
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new());
        }

        let mut inv2 = FpNum::from(2);
        inv2 = inv2.invert(self.f);

        let mut failures = 0;

        for (x, y) in self
            .targets
            .iter()
            .combinations_with_replacement(2)
            .map(|v| (v[0], v[1]))
        {
            let x = FpNum::from(*x);
            let y = FpNum::from(*y);
            let disc = 9 * x * y - 4 * (x * x + y * y);
            let neg_b = 3 * x * y;
            let mut candidates = Vec::new();

            match disc.int_sqrt() {
                Some(FpNum(0)) => {
                    candidates.push(neg_b * inv2);
                }
                Some(root_disc) => {
                    candidates.push((neg_b + root_disc) * inv2);
                    candidates.push((neg_b - root_disc) * inv2);
                }
                None => {
                    failures += 1;
                }
            }

            let it: Vec<(&FpNum<P>, bool)> = candidates
                .iter()
                .map(|z| (z, results.contains_key(&z.0)))
                .collect();
            let Some(disjoint) = results.get_mut(&x.0) else { continue; };
            for (z, pred) in it {
                if pred {
                    disjoint.associate(x.0, z.0);
                    disjoint.associate(y.0, z.0);
                }
            }
        }

        OrbitTesterResults { failures, results }
    }

    pub fn new(f: &FpStar<P>) -> OrbitTester<P> {
        OrbitTester {
            f,
            targets: HashSet::new(),
        }
    }

    pub fn add_target(mut self, t: u128) -> OrbitTester<'a, P> {
        self.targets.insert(t);
        self
    }
}

impl OrbitTesterResults {
    pub fn results(&self) -> impl Iterator<Item = (&u128, &Disjoint<u128>)> {
        self.results.iter()
    }

    pub fn failures(&self) -> u128 {
        self.failures
    }
}
