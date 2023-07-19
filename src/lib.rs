#![allow(dead_code)]
mod factorization;
mod quad_field;
mod util;
mod semigroup;
mod sylow;
mod fp;

#[cfg(test)]
mod tests {
    use crate::factorization::*;

    #[test]
    fn fact_works() {
        let fact = Factorization {
            value: 8,
            factors: vec![8],
            prime_powers: vec![(2, 3)]
        };
        assert_eq!(fact.value, 8);
    }
}
