//! Stock implementations of the `Factor` trait for some commonly used primes.
use crate::numbers::*;

/// Empty marker struct for stock implementations.
pub struct Stock {}

impl Factor<Stock, 7> for FpNum<1_000_000_000_000_000_124_399> {
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

impl Factor<Stock, 11> for QuadNum<1_000_000_000_000_000_124_399> {
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

