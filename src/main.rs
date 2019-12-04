#![allow(non_snake_case)]

use clap::{App, Arg};
use rug::Integer;
use rug::integer::IsPrime;

use MPQS::*;

pub fn main() {
    let app = App::new("MPQS")
        .version("1.0")
        .author("Luca Allegro <luca.all1996@gmail.com>")
        .about("Rust implementation of Multi-polinomial Quadratic Sieve, an algorithm to factorize numbers")
        .arg(Arg::with_name("algorithm")
            .short("a")
            .long("algorithm")
            .value_name("S or M or A")
            .help("S for serial, M for Memory Shared or A for Actor(message passing)")
            .required(true)
            .validator(|v| if v == "S" || v == "M" || v == "A" { Ok(()) } else { Err("Algorithm should be one of S, M or A".to_owned()) })
            .takes_value(true))
        .arg(Arg::with_name("number")
            .short("n")
            .long("number")
            .value_name("N")
            .validator(|v| if v.chars().all(|c| c.is_ascii_digit()) { Ok(()) } else { Err("Number accepts only digits".to_owned()) })
            .required_unless("product")
            .help("The number to factorize")
            .takes_value(true))
        .arg(Arg::with_name("product")
            .short("p")
            .long("product")
            .value_name("N N")
            .help("Insert a pair of number to factorize their product")
            .validator(|v| if v.chars().all(|c| c.is_ascii_digit()) { Ok(()) } else { Err("Number accepts only digits".to_owned()) })
            .number_of_values(2)
            .takes_value(true))
        .get_matches();

    let n: Integer = app
        .value_of("number")
        .map(|n| n.parse::<Integer>().unwrap())
        .unwrap_or_else(|| {
            let numbers: Vec<&str> = app.values_of("product").unwrap().collect();
            numbers[0].parse::<Integer>().unwrap() * numbers[1].parse::<Integer>().unwrap()
        });

    if rabin_miller::is_rabin_miller_prime(&n) || n.is_probably_prime(15) == IsPrime::Probably {
        println!(
            "{} is probably prime. Don't waste time trying to factorize it ;)",
            n
        );
    } else {
        let r = match app.value_of("algorithm").unwrap() {
            "S" => time(|| serial_MPQS::mpqs(&n)),
            "M" => time(|| memory_shared_MPQS::mpqs(&n)),
            "A" => time(|| message_MPQS::mpqs(&n)),
            _ => panic!(""),
        };
        check_is_divisor(n, r);
    }
}

// Results: Factorization di 1201121312171223122912311237 * 3023706637809542222940030043
// Python - 524 s
