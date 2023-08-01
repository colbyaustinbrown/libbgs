use libbgs::triple::*;
use libbgs::numbers::quad_field::*;
use libbgs::numbers::factorization::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::sylow_factory::*;

fn main() {
    let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
    let fp2 = QuadField::make(
        fp.clone(),
        Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
    );
    // let p = 13;
    /*
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
    let decomp = SylowDecomp::new(&fp);
    let factory = SylowFactory::new(&decomp, 1, 1);
    for a in factory.map(|v| v.to_product(&decomp)) {
        println!("a is {a:?}");
        let x = MarkoffTriple::make(a.value, 1, 1, &fp2);
        println!("ord is {}", x.ord(Pos::A, &fp2));
    }
}
