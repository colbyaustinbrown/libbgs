use libbgs::markoff::*;
use libbgs::numbers::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

#[derive(PartialEq, Eq)]
struct Phantom {}

impl Factored<Phantom, 7> for FpNum<BIG_P>{
    const FACTORS: Factorization<7> = Factorization::new([
        (2, 1),
        (7, 1),
        (13, 1),
        (29, 2),
        (43, 1),
        (705737, 1),
        (215288719, 1),
    ]);
}

impl Factored<Phantom, 11> for QuadNum<BIG_P> {
    const FACTORS: Factorization<11> = Factorization::new([
        (2, 4),
        (3, 1),
        (5, 2),
        (11, 2),
        (17, 1),
        (19, 1),
        (23, 1),
        (97, 1),
        (757, 1),
        (1453, 1),
        (8689, 1),
    ]);
}

fn main() {
    let fp = Factorization::new([
        (2, 1),
        (7, 1),
        (13, 1),
        (29, 2),
        (43, 1),
        (705737, 1),
        (215288719, 1),
    ]);
    let fp2_fact = Factorization::new([
        (2, 4),
        (3, 1),
        (5, 2),
        (11, 2),
        (17, 1),
        (19, 1),
        (23, 1),
        (97, 1),
        (757, 1),
        (1453, 1),
        (8689, 1),
    ]);
    let fp_decomp = SylowDecomp::new();
    let fp2_decomp = SylowDecomp::new();

    const LIMIT: u128 = 100;

    let mut fp_stream_builder = SylowStreamBuilder::new(&fp_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    let mut fp2_stream_builder = SylowStreamBuilder::new(&fp2_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    println!("Adding these targets from F_p: ");
    for d in fp.maximal_divisors(LIMIT) {
        println!("\t{d:?}");
        fp_stream_builder = fp_stream_builder.add_target(d);
    }
    println!("Adding these targets from F_p^2: ");
    for d in fp2_fact.maximal_divisors(LIMIT) {
        println!("\t{d:?}");
        fp2_stream_builder = fp2_stream_builder.add_target(d);
    }
    let stream = fp_stream_builder
        .build()
        .map(|x| Coord::from_chi_fp(&x, &fp_decomp))
        .chain(
            fp2_stream_builder
                .build()
                .map(|x| Coord::from_chi_quad(&x, &fp2_decomp)),
        );

    let mut tester = OrbitTester::<BIG_P>::new();
    let mut count = 0;
    println!("Loading coordinates into the Orbit Tester.");
    for x in stream {
        count += 1;
        tester = tester.add_target(u128::from(x));
    }
    println!("Loaded {count} coordinates into the Orbit Tester.");

    println!("Running the Orbit Tester.");
    let results = tester.run();
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
    println!("There were {} pairs of coordinates which did not form the first two coordinates in the solution of the Markoff equation.", results.failures());
}
