#![allow(dead_code)]
pub mod factorization;
pub mod quad_field;
pub mod util;
pub mod semigroup;
pub mod sylow;
pub mod fp;
pub mod sylow_factory;
pub mod markoff_triple;
pub mod group;

#[cfg(test)]
mod tests {
    use crate::factorization::*;

    #[test]
    fn fact_works() {
        let fact = Factorization::new(
            vec![(2, 3)]
        );
        assert_eq!(fact.value(), 8);
    }
}
