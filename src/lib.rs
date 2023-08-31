#![allow(dead_code)]
pub mod markoff;
pub mod numbers;
pub mod util;

#[cfg(test)]
mod tests {
    use crate::numbers::*;

    #[test]
    fn fact_works() {
        let fact = Factorization::new([(2, 3)]);
        assert_eq!(fact.value(), 8);
    }
}
