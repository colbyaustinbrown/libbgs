use libbgs::markoff::coord::*;
use libbgs::markoff::triple::*;
use libbgs::numbers::factorization::*;
use libbgs::numbers::group::*;
use libbgs::numbers::quad_field::*;

fn main() {
    let p = 13;
    let fp = Factorization::new(vec![(2, 2), (3, 1)]);
    let fp2 = QuadField::make(fp, Factorization::new(vec![(2, 1), (7, 1)]));
    for a in 0..p {
        for b in 0..p {
            for c in 0..p {
                if (a*a + b*b + c*c + 3*a*b*c*p - 3*a*b*c) % p != 0 {
                    continue;
                }
                let x = MarkoffTriple::make(a, b, c, &fp2);
                for m in x.rot(Pos::A, &fp2) {
                    println!(
                        "{} {} {} ({})",
                        m.a.v(),
                        m.b.v(),
                        m.c.v(),
                        m.ord(Pos::C, &fp2)
                    );
                }
                println!();
            }
        }
    }
}
