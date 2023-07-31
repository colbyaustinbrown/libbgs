use std::rc::Rc;

use libbgs::triple::*;
use libbgs::numbers::quad_field::*;
use libbgs::numbers::factorization::*;

fn main() {
    println!("hello, world?");
    let six = Factorization::new(vec![(2, 1), (3, 1)]);
    let eight = Factorization::new(vec![(2, 3)]);
    let f = Rc::new(QuadFieldExt::make(six, eight));
    let x = MarkoffTriple::make(&f, 1, 1, 1);
    println!("x looks like {:?}", x);
    println!("order of x is {}", x.ord(Pos::A));
}
