mod factorization;
mod quad_field;
mod util;
mod semigroup;
mod sylow;

#[cfg(test)]
mod tests {
    use crate::factorization::*;
    use crate::quad_field::*;
    use crate::sylow::*;

    #[test]
    fn fact_works() {
        let fact = Factorization {
            value: 8,
            primepowers: vec![(2, 3, 8)]
        };
        assert_eq!(fact.value, 8);
    }
}
