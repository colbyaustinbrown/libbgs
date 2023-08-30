#![allow(dead_code)]
#![warn(missing_docs)]

//! A library for the manipulation and investigation of Markoff numbers and Markoff graphs modulo
//! primes, specifically with a focus on the Bourgain, Gamburd, and Sarnak algorithm.
pub mod markoff;
pub mod numbers;
pub mod util;

#[cfg(test)]
mod tests {
    use crate::numbers::*;

    #[test]
    fn fact_works() {
        let fact = Factorization::new(vec![(2, 3)]);
        assert_eq!(fact.value(), 8);
    }
}
