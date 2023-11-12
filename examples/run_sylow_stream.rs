use libbgs::numbers::*;
use libbgs::streams::*;

use std::time::*;

use rayon::iter::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

struct Stock {}
impl_factors!(Stock, 1_000_000_000_000_000_124_399);

fn main() {
    // let decomp = SylowDecomp::<Phantom, 7, FpNum<BIG_P>>::new();
    let factors = FactorStream::new(FpNum::<BIG_P>::FACTORS.prime_powers(), 10_000_000, true);

    let builder = SylowStreamBuilder::<Stock, 7, FpNum<BIG_P>>::new()
        .add_flag(flags::LEQ)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::NO_UPPER_HALF)
        .add_targets_from_factors(factors);

    let stream = builder.into_par_iter();

    println!("starting the stream");
    let now = Instant::now();
    let count = stream.count();
    let elapsed = now.elapsed();

    println!("stream count is {count}");
    println!("timer says {} seconds", elapsed.as_secs());
}
