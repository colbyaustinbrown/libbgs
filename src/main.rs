use libbgs::triple::*;
use libbgs::numbers::quad_field::*;
use libbgs::numbers::factorization::*;

fn main() {
        let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
        let fp2 = QuadField::make(
            fp,
            Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
        );
    let x = MarkoffTriple::make(1, 1, 1, &fp2);
    println!("x looks like {:?}", x);
    println!("order of x is {}", x.ord(Pos::A, &fp2));
}
