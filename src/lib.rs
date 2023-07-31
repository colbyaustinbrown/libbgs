#![allow(dead_code)]
pub mod numbers;
pub mod util;
pub mod coord;
pub mod triple;

#[cfg(test)]
mod tests {
    use crate::numbers::factorization::*;

    #[test]
    fn fact_works() {
        let fact = Factorization::new(
            vec![(2, 3)]
        );
        assert_eq!(fact.value(), 8);
    }
}
