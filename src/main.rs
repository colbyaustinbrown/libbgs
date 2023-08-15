use either::Either::*;

use libbgs::numbers::fp::*;
use libbgs::numbers::quad_field::*;
use libbgs::numbers::factorization::*;
use libbgs::numbers::sylow::*;
use libbgs::numbers::sylow_stream::*;
// use libbgs::numbers::group::*;
use libbgs::markoff::coord::*;
// use libbgs::markoff::triple::*;
use libbgs::orbit_tester::*;
    
const BIG_P: u128 = 1_000_000_000_000_000_124_399;

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
    /*
    let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
    let fp2 = QuadField::make(
        fp.clone(),
        Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)])
    );
    println!("p is {}", fp2.p());
    let decomp = SylowDecomp::new(&fp);
    //let stream = sylow_stream(&decomp, &vec![1, 1, 1, 0, 0, 0, 0], Mode::LEQ);
    let mut count = 0;
    let builder = SylowStreamBuilder::new(&decomp)
        .add_target(vec![0, 1, 1, 0, 0, 0, 0])
        .add_target(vec![0, 0, 1, 1, 0, 0, 0])
        // .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::LEQ);

    let tester_builder = OrbitTester::new(&fp);
    for x in builder.build() { 

    }
    for y in builder.build() {
        let mut yin = y.clone();
        // println!("{:?}", y.coords);
        yin.invert(&decomp);
        let a = Coord::from_chi(Right((y.to_product(&decomp), yin.to_product(&decomp))), &fp2);
        println!("{} has order {}", a.v(), a.get_ord(&fp2).value());
        count += 1;
    }
    println!("total: {count}");
    */

    let fp = Factorization::new(vec![(2, 1), (7, 1), (13, 1), (29, 2), (43, 1), (705737, 1), (215288719, 1)]);
    let fp2_fact = Factorization::new(vec![(2, 4), (3, 1), (5, 2), (11, 2), (17, 1), (19, 1), (23, 1), (97, 1), (757, 1), (1453, 1), (8689, 1)]);
    let fp2 = QuadField::<BIG_P>::make();
    let fp_decomp = SylowDecomp::new(&FpStar::<BIG_P> {}, fp.clone());
    let fp2_decomp = SylowDecomp::new(&fp2, fp2_fact.clone());

    const LIMIT: u128 = 500;

    let mut fp_stream_builder = SylowStreamBuilder::new(&fp_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    let mut fp2_stream_builder = SylowStreamBuilder::new(&fp2_decomp)
        .add_flag(flags::NO_UPPER_HALF)
        .add_flag(flags::NO_PARABOLIC)
        .add_flag(flags::LEQ);
    println!("Adding these targets from F_p: ");
    for d in fp.maximal_divisors(LIMIT) {
        println!("\t{d:?}");
        fp_stream_builder = fp_stream_builder.add_target(d);
    }
    println!("Adding these targets from F_p^2: ");
    for d in fp2_fact.maximal_divisors(LIMIT) {
        println!("\t{d:?}");
        fp2_stream_builder = fp2_stream_builder.add_target(d);
    }
    let stream = fp_stream_builder.build()
        .map(|x| Right((x.clone(), x)))
        .chain(fp2_stream_builder.build()
            .map(|x| Left((x.clone(), x)))
        )
        .map(|e| e.either(
            |(l, mut l_i)| { 
                l_i = l_i.invert(&fp2_decomp);
                Left((
                    l.to_product(&fp2_decomp), 
                    l_i.to_product(&fp2_decomp)
                ))
            },
            |(r, mut r_i)| { 
                r_i = r_i.invert(&fp_decomp);
                Right((
                    r.to_product(&fp_decomp), 
                    r_i.to_product(&fp_decomp)
                ))
            }
        ))
        .map(Coord::from_chi);

    let mut tester = OrbitTester::new(&FpStar::<BIG_P> {});
    let mut count = 0;
    println!("Loading coordinates into the Orbit Tester.");
    for x in stream {
        count += 1;
        tester = tester.add_target(x.v());
    }
    println!("Loaded {count} coordinates into the Orbit Tester.");

    println!("Running the Orbit Tester.");
    let results = tester.run();
    println!("Testing complete.");

    let mut repless_count = 0;
    let mut big_count = 0;
    for (x, disjoint) in results.results() {
        let mut orbits = disjoint.get_orbits().peekable();

        if orbits.peek().is_none() {
            repless_count += 1;
            continue;
        }

        for (key, set) in orbits {
            if set.data {
                println!("For coordinate {x}: Representative: {key} is an orbit of small coordinates."); 
            } else {
                big_count += 1;
            }
        }
    }
    println!("{repless_count} coordinates had no representative.");
    println!("{big_count} orbits had big coordinates.");
    println!("There were {} pairs of coordinates which did not form the first two coordinates in the solution of the Markoff equation.", results.failures());
}

