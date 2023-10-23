//! Stock implementations of the `Factor` trait for some commonly used primes.
use crate::numbers::*;
use crate::impl_factors;

/// Empty marker struct for stock implementations.
#[derive(Clone)]
pub struct Stock {}

impl_factors!(Stock, 100_003);
impl_factors!(Stock, 1_000_033..1_000_100);

impl Factor<Stock> for FpNum<1_234_577> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 4),
        (7, 1),
        (73, 1),
        (151, 1),
    ]);
}

impl Factor<Stock> for QuadNum<1_234_577> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (3, 1),
        (205_763, 1),
    ]);
}

impl Factor<Stock> for FpNum<1_500_007> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (3, 1),
        (53, 2),
        (89, 1),
    ]);
}

impl Factor<Stock> for QuadNum<1_500_007> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 3),
        (97, 1),
        (1933, 1),
    ]);
}

impl Factor<Stock> for FpNum<4_999_999> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (3, 1),
        (191, 1),
        (4363, 1),
    ]);
}

impl Factor<Stock> for QuadNum<4_999_999> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 6),
        (5, 7),
    ]);
}

impl Factor<Stock> for FpNum<10_000_139> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (37, 1),
        (337, 1),
        (401, 1),
    ]);
}

impl Factor<Stock> for QuadNum<10_000_139> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 2),
        (3, 1),
        (5, 1),
        (166_669, 1),
    ]);
}

impl Factor<Stock> for FpNum<101_010_079> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (3, 4),
        (13, 1),
        (47963, 1),
    ]);
}

impl Factor<Stock> for QuadNum<101_010_079> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 5),
        (5, 1),
        (19, 1),
        (149, 1),
        (223, 1),
    ]);
}

impl Factor<Stock> for FpNum<1_000_000_000_000_000_124_399> {
    const FACTORS: Factorization = Factorization::new(&[
        (2, 1),
        (7, 1),
        (13, 1),
        (29, 2),
        (43, 1),
        (705737, 1),
        (215288719, 1),
    ]);
}

impl Factor<Stock> for QuadNum<1_000_000_000_000_000_124_399> {
    const FACTORS: Factorization = Factorization::new(&[
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

