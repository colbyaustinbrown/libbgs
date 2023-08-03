use criterion::{black_box, criterion_group, criterion_main, Criterion};

use libbgs::numbers::factorization::*;
use libbgs::numbers::quad_field::*;

fn find_sqrt(n: u128) {
    let fp = QuadField::make(
        Factorization::new(
            vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]
        ),
        Factorization::new(
            vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]
        )
    );
    for i in 3..103 {
        let mut x = fp.int_sqrt(i);
        let y = x.clone();
        assert_ne!(x, i);
        x.multiply(&y, &fp);
        assert_eq!(x, i);
    }
}

fn criterion_benchmark(c: &mut Criterion ){
    c.bench_function("find sqrt", |b| b.iter(|| find_sqrt(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
