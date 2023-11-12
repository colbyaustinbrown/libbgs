use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use libbgs::markoff::*;
use libbgs::numbers::*;
use libbgs::streams::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

struct Stock {}

impl_factors!(Stock, 1_000_000_000_000_000_124_399);

fn main() {
    let fp_decomp = SylowDecomp::<Stock, 7, FpNum<BIG_P>>::new();
    let fp2_decomp = SylowDecomp::<Stock, 11, QuadNum<BIG_P>>::new();

    const LIMIT: u128 = 10_000;

    let mut fp_stream_builder = SylowStreamBuilder::new()
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    let mut fp2_stream_builder = SylowStreamBuilder::new()
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    println!("Adding these targets from F_p: ");
    for d in FpNum::<BIG_P>::FACTORS.maximal_divisors(LIMIT) {
        println!("\t{d:?}");
        fp_stream_builder = fp_stream_builder.add_target(&d);
    }
    println!("Adding these targets from F_p^2: ");
    for d in QuadNum::<BIG_P>::FACTORS.maximal_divisors::<{ QuadNum::<BIG_P>::LENGTH }>(LIMIT) {
        println!("\t{d:?}");
        fp2_stream_builder = fp2_stream_builder.add_target(&d);
    }

    let tester = Mutex::new(OrbitTester::<BIG_P>::new());
    let count = AtomicUsize::new(0);
    println!("Loading coordinates into the Orbit Tester.");
    fp_stream_builder
        .into_iter()
        .map(|(x, _)| Coord::from_chi(&x, &fp_decomp))
        .chain(
            fp2_stream_builder
                .into_iter()
                .map(|(x, _)| Coord::from_chi(&x, &fp2_decomp)),
        )
        .for_each(|x| {
            count.fetch_add(1, Ordering::Relaxed);
            tester.lock().unwrap().add_target(u128::from(x));
        });

    println!(
        "Loaded {} coordinates into the Orbit Tester.",
        count.into_inner()
    );

    println!("Running the Orbit Tester.");
    let results = tester.into_inner().unwrap().run();
    println!("Testing complete.");

    let mut repless_count = 0;
    for (x, disjoint) in results.results() {
        let mut orbits = disjoint.get_sets().peekable();

        if orbits.peek().is_none() {
            repless_count += 1;
            continue;
        }

        for (key, ord) in orbits {
            println!("For coordinate {x}: Representative: {key} is an orbit of order {ord}.");
        }
    }
    println!("{repless_count} coordinates had no representative.");
}
