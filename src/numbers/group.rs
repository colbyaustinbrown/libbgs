use crate::numbers::semigroup::*;

pub trait Group: Semigroup {}

pub trait GroupElem: SemigroupElem {
    fn invert(&self, g: &Self::Group) -> Self;
}
