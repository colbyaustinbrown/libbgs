use crate::semigroup::*;

pub trait Group: Semigroup {}

pub trait GroupElem: SemigroupElem {
    fn invert(&mut self);
}

