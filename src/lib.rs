#![allow(dead_code)]
#![warn(missing_docs)]
pub mod markoff;
pub mod numbers;
pub mod util;

#[cfg(test)]
mod tests {
    use crate::numbers::factorization::*;

    #[test]
    fn fact_works() {
        let fact = Factorization::new(vec![(2, 3)]);
        assert_eq!(fact.value(), 8);
    }
}
