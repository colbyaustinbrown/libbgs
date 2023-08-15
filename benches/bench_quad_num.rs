use criterion::{criterion_group, criterion_main, Criterion};

use libbgs::numbers::factorization::*;
use libbgs::numbers::quad_field::*;

fn find_sqrt() {
    let fp = QuadField::<1_000_000_000_000_000_124_399>::make();
    for i in 3..103 {
        let mut x = fp.int_sqrt(i);
        let y = x.clone();
        assert_ne!(x, i);
        x = x.multiply(&y, &fp);
        assert_eq!(x, i);
    }
}

fn criterion_benchmark(c: &mut Criterion ){
    c.bench_function("find sqrt", |b| b.iter(|| find_sqrt()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
