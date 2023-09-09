use std::collections::{HashMap, HashSet};
use std::thread;

use itertools::*;
use rayon::prelude::*;

use crate::markoff::Disjoint;
use crate::numbers::{FpNum, GroupElem};

/// Configures tests to be run on orbits of the Markoff graph modulo `P`.
pub struct OrbitTester<const P: u128> {
    targets: HashSet<u128>,
}

/// The results of a successfully run `OrbitTester`.
pub struct OrbitTesterResults {
    results: HashMap<u128, Disjoint<u128>>,
}

type Msg = (u128, u128, u128);

impl<const P: u128> OrbitTester<P> {
    /// Consume and run this `OrbitTester`, blocking until completion, and returning the results.
    /// This method may spawn multiple worker threads, which are guarenteed to be joined before
    /// `run` returns.
    pub fn run(self) -> OrbitTesterResults {
        let mut results = HashMap::with_capacity(self.targets.len());
        for x in &self.targets {
            results.insert(*x, Disjoint::new());
        }

        let mut inv2 = FpNum::<P>::from(2);
        inv2 = inv2.inverse();

        let (tx, rx) = std::sync::mpsc::sync_channel::<Msg>(1024);

        let handle = thread::spawn(move || {
            for (x, y, z) in rx.iter() {
                if results.contains_key(&z) {
                    if let Some(disjoint) = results.get_mut(&x) {
                        disjoint.associate(y, y);
                    }
                    if let Some(disjoint) = results.get_mut(&y) {
                        disjoint.associate(x, z);
                    }
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

                // We use the non-normalized equation: x^2 + y^2 + z^2 - xyz = 0
                let disc = x * y - 4 * (x * x + y * y);
                let neg_b = x * y;

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
                    None => {}
                }
            });
        drop(tx);

        let results = handle.join().unwrap();

        OrbitTesterResults {
            results,
        }
    }

    /// Creates a new `OrbetTester` with default settings and no targets.
    pub fn new() -> OrbitTester<P> {
        OrbitTester {
            targets: HashSet::new(),
        }
    }

    /// Adds a target order to the list of orders to be tested.
    pub fn add_target(&mut self, t: u128) {
        self.targets.insert(t);
    }
}

impl OrbitTesterResults {
    /// The results of the test, as an iterator yielding each coordinate of a target order, along
    /// with the partitioning of the target orders into disjoint sets, which are subsets of the
    /// orbits under the fixed first coordinate.
    pub fn results(&self) -> impl Iterator<Item = (&u128, &Disjoint<u128>)> {
        self.results.iter()
    }
}
