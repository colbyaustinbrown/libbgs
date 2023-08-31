use libbgs::markoff::*;
use libbgs::numbers::*;

const P: u128 = 13;

fn main() {
    let fp_facts = Factorization::new(vec![(2, 2), (3, 1)]);
    let fp2 = QuadField::<P> {};
    let fp2_facts = Factorization::new(vec![(2, 1), (7, 1)]);
    for a in 0..P {
        for b in 0..P {
            for c in 0..P {
                if (a*a + b*b + c*c + 3*a*b*c*P - 3*a*b*c) % P != 0 {
                    continue;
                }
                for (y, z) in Coord::from(a).rot(Coord::from(b), Coord::from(c)) {
                    println!(
                        "{a} {} {} ({})",
                        y,
                        z,
                        z.get_ord(&fp2, &fp_facts, &fp2_facts).value()
                    );
                }
                println!();
            }
        }
    }
}
