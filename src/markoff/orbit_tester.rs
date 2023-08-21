use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;

use itertools::*;
use rayon::prelude::*;

use crate::markoff::disjoint::Disjoint;
use crate::numbers::fp::*;

pub struct OrbitTester<'a, const P: u128> {
    f: &'a FpStar<P>,
    targets: HashSet<u128>,
}

pub struct OrbitTesterResults {
    failures: u64,
    results: HashMap<u128, Disjoint<u128>>,
}

type Msg = (u128, u128, u128);

impl<'a, const P: u128> OrbitTester<'a, P> {
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new());
        }

        let mut inv2 = FpNum::from(2);
        inv2 = inv2.invert(self.f);

        // TOOD: is it a problem that this is a u64 and not a u128?
        let failures = AtomicU64::new(0);
        let (tx, rx) = std::sync::mpsc::sync_channel::<Msg>(1024);

        let handle = thread::spawn(move || {
            for (x, y, z) in rx.iter() {
                if results.contains_key(&z) {
                    results.get_mut(&x).map(|disjoint| {
                        disjoint.associate(y, z);
                    });
                    results.get_mut(&y).map(|disjoint| {
                        disjoint.associate(x, z);
                    });
                }
            }

            results
        });

        self.targets
            .iter()
            .combinations_with_replacement(2)
            .map(|v| (v[0], v[1]))
            .par_bridge()
            .for_each(|(x, y)| {
                let x = FpNum::from(*x);
                let y = FpNum::from(*y);
                let disc = 9 * x * y - 4 * (x * x + y * y);
                let neg_b = 3 * x * y;

                match disc.int_sqrt() {
                    Some(FpNum(0)) => {
                        let z = neg_b * inv2;
                        _ = tx.send((x.0, y.0, z.0));
                    }
                    Some(root_disc) => {
                        let z = (neg_b + root_disc) * inv2;
                        _ = tx.send((x.0, y.0, z.0));
                        let z = (neg_b - root_disc) * inv2;
                        _ = tx.send((x.0, y.0, z.0));
                    }
                    None => {
                        failures.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
        drop(tx);

        let results = handle.join().unwrap();

        OrbitTesterResults {
            failures: failures.into_inner(),
            results,
        }
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

    pub fn failures(&self) -> u64 {
        self.failures
    }
}
