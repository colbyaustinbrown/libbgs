use crate::numbers::semigroup::*;

pub trait Group: Semigroup {}

pub trait GroupElem: SemigroupElem {
    fn invert(&mut self, g: &Self::Group);
}

