use either::Either::*;

use libbgs::triple::*;
use libbgs::numbers::quad_field::*;
use libbgs::numbers::factorization::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::sylow_factory::*;
use libbgs::numbers::group::*;
use libbgs::coord::*;

fn main() {
    /*
    let p = 13;
    let fp = Factorization::new(vec![(2, 2), (3, 1)]);
    let fp2 = QuadField::make(
        fp,
        Factorization::new(vec![(2, 1), (7, 1)])
    );
    for a in 0..p {
        for b in 0..p {
            for c in 0..p {
                if (a*a + b*b + c*c + 3*a*b*c*p - 3*a*b*c) % p != 0 { continue; }
                let x = MarkoffTriple::make(a, b, c, &fp2);
                for m in x.rot(Pos::A, &fp2) {
                    println!("{} {} {} ({})", m.a.v(), m.b.v(), m.c.v(), m.ord(Pos::C, &fp2));
                }
                println!();
            }
        }
    }
    */
    let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
    let fp2 = QuadField::make(
        fp.clone(),
        Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
    );
    println!("p is {}", fp2.p());
    let decomp = SylowDecomp::new(&fp);
    //let factory = sylow_factory(&decomp, &vec![1, 1, 1, 0, 0, 0, 0], Mode::LEQ);
    let factory = SylowFactory::new(&decomp, 3, 1, flags::NO_UPPER_HALF);
    let mut count = 0;
    for y in factory {
        let mut yin = y.clone();
        yin.invert(&decomp);
        let a = Coord::from_chi(Right((y.to_product(&decomp), yin.to_product(&decomp))), &fp2);
        println!("{} has order {}", a.v(), a.get_ord(&fp2).value());
        count += 1;
    }
    println!("total: {count}");
}
