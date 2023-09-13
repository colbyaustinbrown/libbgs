use libbgs::numbers::*;

use std::time::*;

use rayon::iter::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;
// const BIG_P_SQRT: u128 = 31_622_776_602;

struct Phantom {}

impl Factor<Phantom, 7> for FpNum<BIG_P> {
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

fn main() {
    // let decomp = SylowDecomp::<Phantom, 7, FpNum<BIG_P>>::new();
    let mut builder = SylowStreamBuilder::<Phantom, 7, FpNum<BIG_P>>::new()
        .add_flag(flags::LEQ)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::NO_UPPER_HALF);

    builder = FactorStream::<7>::new(
            FpNum::<BIG_P>::FACTORS.prime_powers(),
            10_000_000,
            true
        )
        .into_iter()
        .inspect(|x| println!("adding target {x:?}"))
        .fold(builder, |b, x| b.add_target(x));

    let stream = builder
        .into_par_iter();

    println!("starting the stream");
    let now = Instant::now();
    let count = stream.count();
    let elapsed = now.elapsed();

    println!("stream count is {count}");
    println!("timer says {} seconds", elapsed.as_secs());
}