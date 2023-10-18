use libbgs::markoff::*;
use libbgs::numbers::*;

const P: u128 = 13;

#[derive(PartialEq, Eq)]
struct Phantom {}

fn main() {
    for a in 0..P {
        for b in 0..P {
            for c in 0..P {
                if (a*a + b*b + c*c + 3*a*b*c*P - 3*a*b*c) % P != 0 {
                    continue;
                }
                for (y, z) in Coord::<13>::from(a).rot(Coord::from(b), Coord::from(c)) {
                    println!(
                        "{a} {} {} ({})",
                        u128::from(y),
                        u128::from(z),
                        z.rot_order()
                    );
                }
                println!();
            }
        }
    }
}

impl Factor<Phantom> for FpNum<13> {
    const FACTORS: Factorization = Factorization::new(&[(2, 2), (3, 1)]);
}

impl Factor<Phantom> for QuadNum<13> {
    const FACTORS: Factorization = Factorization::new(&[(2, 1), (7, 1)]);
}
