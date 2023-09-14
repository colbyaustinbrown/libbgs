use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::*;

const BIG_P: u128 = 1_000_000_000_000_000_124_399;

fn find_inverses() {
    for i in 1..100 {
        let x = FpNum::<BIG_P>::from(i);
        let y = QuadNum::<BIG_P>::from((0, i));
        x.inverse();
        y.inverse();
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("find inverses", |b| b.iter(|| find_inverses()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
