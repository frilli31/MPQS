#![allow(non_snake_case)]
#![allow(unused_must_use)]

use std::fmt;
use std::time::Instant;

use rug::Integer;

use rabin_miller::is_rabin_miller_prime;

pub mod algebra;
//pub mod memory_shared_MPQS;
//pub mod message_MPQS;
pub mod rabin_miller;
pub mod serial_MPQS;
pub mod tonelli_shanks;

pub fn main() {
    //std::env::set_var("RUST_LOG", "Factorization::algebra=INFO");
    //env_logger::init();

    let _s1 = "676292275716558246502605230897191366469551764092181362779759";
    let _s2 = "2736300383840445596906210796102273501547527150973747";

    let _p1 = "1201121312171223122912311237".parse::<Integer>().unwrap();
    let _p2 = "3023706637809542222940030043".parse::<Integer>().unwrap();

    let p1 = "98761037 7233144895 5342113853".parse::<Integer>().unwrap();
    let p2 = "1001446553 1244957205 9845328443"
        .parse::<Integer>()
        .unwrap();
    //let p3 = "4825641527 1247992250 9389813571".parse::<Integer>().unwrap();

    let _n = _s2.parse::<Integer>().unwrap();

    let _r = time(|| serial_MPQS::mpqs(&(_p1.clone() * &_p2)));
    //    let _r = time(|| message_MPQS::mpqs(&(_p1.clone() * &_p2)));
    //    let _r = time(|| memory_shared_MPQS::mpqs(&(_p1.clone() * &_p2)));
}

pub fn check_is_divisor(n: Integer, qs: Option<Integer>) {
    match qs {
        Some(qs) => {
            let (q, r) = n.clone().div_rem(qs.clone());
            assert_eq!(r, Integer::new());
            println!("{} = {} * {}", n, qs, q);
        }
        None => {
            assert!(is_rabin_miller_prime(&n));
            println!("{} is Rabin Miller prime", n);
        }
    }
}

pub fn time<F, O: fmt::Debug>(f: F) -> O
where
    F: Fn() -> O,
{
    let t = Instant::now();
    let ris = f();
    println!("R: {:?} in {:?}", ris, t.elapsed());
    ris
}

// Results: Factorization di 1201121312171223122912311237 * 3023706637809542222940030043
// Python - 524 s
// Rust Serial - 39 s
// Rust Parallel - 20 s
//
// >>> 524 / 39
// 13.435897435897436
// >>> 524 / 20
// 26.2
// >>>
