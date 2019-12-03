#![allow(non_snake_case)]
#![allow(unused_must_use)]

use std::fmt;
use std::time::Instant;

use rug::Integer;

use crate::rabin_miller::is_rabin_miller_prime;

pub mod algebra;
pub mod memory_shared_MPQS;
pub mod message_MPQS;
pub mod rabin_miller;
pub mod serial_MPQS;
pub mod tonelli_shanks;

pub fn modular_inv(a0: Integer, m0: Integer) -> Integer {
    if m0 == 1 {
        return Integer::from(1);
    }

    let (mut a, mut m, mut x0, mut inv) = (a0, m0.clone(), Integer::new(), Integer::from(2));

    while a > 1 {
        inv = inv - (a.clone() / &m) * &x0;
        a %= &m;
        std::mem::swap(&mut a, &mut m);
        std::mem::swap(&mut x0, &mut inv);
    }
    if inv < 0 {
        inv += m0;
    }
    inv
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
