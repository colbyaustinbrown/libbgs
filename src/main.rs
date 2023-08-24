use libbgs::markoff::triple::*;
use libbgs::numbers::factorization::*;
use libbgs::numbers::quad_field::*;

const P: u128 = 13;

fn main() {
    let fp_facts = Factorization::new(vec![(2, 2), (3, 1)]);
    let fp2 = QuadField::<P>::make();
    let fp2_facts = Factorization::new(vec![(2, 1), (7, 1)]);
    for a in 0..P {
        for b in 0..P {
            for c in 0..P {
                if (a*a + b*b + c*c + 3*a*b*c*P - 3*a*b*c) % P != 0 {
                    continue;
                }
                let x = MarkoffTriple::make(a, b, c, &fp2);
                for m in x.rot(Pos::A, &fp2) {
                    println!(
                        "{} {} {} ({})",
                        m.a.v(),
                        m.b.v(),
                        m.c.v(),
                        m.ord(Pos::C, &fp2, &fp_facts, &fp2_facts)
                    );
                }
                println!();
            }
        }
    }
}
