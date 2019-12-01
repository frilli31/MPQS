use std::cmp::{max, min};
use std::collections::HashMap;

use log::info;
use primal_sieve;
use rug::Integer;

use crate::algebra;
use crate::tonelli_shanks::tonelli_shanks;

pub fn mpqs(n: &Integer) -> Option<Integer> {
    let InitResult {
        mut roota,
        factorbase,
        tsqrt,
        xmax,
        tlog,
        thresh,
        min_prime,
    } = initialize_qs(n);

    let mut smooths = Vec::new();
    let mut partials: HashMap<Integer, (Integer, (Integer, Integer))> = HashMap::new();
    let sievesize = 1_i64 << 15;

    loop {
        loop {
            roota.next_prime_mut();
            if n.legendre(&roota) == 1 {
                break;
            }
        }
        info!("Loop 1, roota: {}", roota);
        let a = roota.clone() * &roota;
        let b = tonelli_shanks(&n, &roota);

        let int2: Integer = b.clone() * 2;
        let intermediate = int2.invert(&roota).unwrap();

        let b = (b.clone() - (b.clone() * b - n) * intermediate) % &a;

        let c = (b.clone() * &b - n) / &a;

        info!("a={} \t b={} \t c={}", a, b, c);

        let mut s1: HashMap<u64, i64> = HashMap::new();
        let mut s2: HashMap<u64, i64> = HashMap::new();

        for (i, p) in factorbase.iter().enumerate() {
            let ainv = a
                .clone()
                .pow_mod(&Integer::from(p - 2), &Integer::from(*p))
                .unwrap();
            let mut sol1 = (tsqrt[i].clone() - &b) * &ainv % p;
            let mut sol2 = (-tsqrt[i].clone() - &b) * &ainv % p;
            sol1 -= ((sol1.clone() + xmax) / p) * p;
            sol2 -= ((sol2.clone() + xmax) / p) * p;

            s1.insert(*p, (sol1 + xmax).to_i64().unwrap());
            s2.insert(*p, (sol2 + xmax).to_i64().unwrap());
        }

        for low in (-xmax..xmax + 1).step_by(sievesize as usize + 1) {
            let high = min(xmax, low + sievesize);
            let size = high - low;
            let size_plus_1 = size + 1;

            let mut S = vec![0_f64; size_plus_1 as usize];

            for (i, p) in factorbase.iter().enumerate() {
                if *p < min_prime {
                    continue;
                }
                let mut sol1 = s1[p];
                let mut sol2 = s2[p];
                let logp = tlog[i];

                let p_i64 = *p as i64;
                while sol1 <= size || sol2 <= size {
                    if sol1 <= size {
                        S[sol1 as usize] += logp;
                        sol1 += p_i64;
                    }
                    if sol2 <= size {
                        S[sol2 as usize] += logp;
                        sol2 += p_i64;
                    }
                }
                s1.insert(*p, sol1 - size_plus_1);
                s2.insert(*p, sol2 - size_plus_1);
            }

            for i in 0..size_plus_1 {
                if S[i as usize] > thresh {
                    let x = i + low;
                    let tofact: Integer = a.clone() * x.pow(2) + b.clone() * x * 2 + &c;
                    let mut nf = tofact.clone().abs();

                    for p in factorbase.iter() {
                        while nf.clone() % p == 0 {
                            nf /= p;
                        }
                    }
                    if nf == 1 {
                        smooths.push((a.clone() * x + &b, (tofact, roota.clone())));
                    } else {
                        match partials.remove(&nf) {
                            Some((pairv, pairvals)) => {
                                smooths.push((
                                    pairv * (a.clone() * x + &b),
                                    (tofact * pairvals.0, pairvals.1 * &roota * nf),
                                ));
                            }
                            None => {
                                partials.insert(nf, (a.clone() * x + &b, (tofact, roota.clone())));
                            }
                        }
                    }
                }
            }
        }
        if smooths.len() > factorbase.len() {
            break;
        }
        info!("{} relations found using {}", smooths.len(), roota);
    }
    info!("{:?}", smooths);
    algebra::algebra(factorbase, smooths, n)
}

pub struct InitResult {
    pub roota: Integer,
    pub factorbase: Vec<u64>,
    pub tsqrt: Vec<Integer>,
    pub xmax: i64,
    pub tlog: Vec<f64>,
    pub thresh: f64,
    pub min_prime: u64,
}

pub fn initialize_qs(n: &Integer) -> InitResult {
    let _root2n: Integer = (n * Integer::from(2)).sqrt();

    info!("Isqrt is {}", _root2n);

    let bound: usize = (n.to_f64().log10().powi(2) * 5_f64) as usize;

    info!("Bound is {}", bound);

    let factorbase: Vec<u64> = primal_sieve::Sieve::new(bound)
        .primes_from(2)
        .take_while(|x| x <= &bound)
        .filter(|x| n.legendre(&Integer::from(*x as u64)) == 1 || *x == 2)
        .map(|x| x as u64)
        .collect();

    info!(
        "Largest prime used is {:?}",
        factorbase[factorbase.len() - 1]
    );
    info!("Factorbase len {:?}", factorbase.len());

    let (mut tsqrt, tlog): (Vec<Integer>, Vec<f64>) = factorbase
        .iter()
        .map(|p| (tonelli_shanks(&n, &Integer::from(*p)), (*p as f64).log10()))
        .unzip();
    tsqrt[0] = Integer::new();

    info!("{:?}", tsqrt);
    info!("{:?}", tlog);

    let xmax: i64 = factorbase.len() as i64 * 60 * 4;
    let mval: Integer = (_root2n.clone() * xmax) >> 1;
    let thresh = mval.to_f64().log10() * 0.735;

    info!("XMAX is {:?}", xmax);
    info!("mval is {:?}", mval);
    info!("thresh is {:?}", thresh);

    let min_prime = (thresh * 3_f64) as u64;
    let fudge: f64 = factorbase
        .iter()
        .take_while(|p| p < &&(min_prime))
        .enumerate()
        .map(|(i, _)| tlog[i])
        .sum::<f64>()
        / 4_f64;

    info!("min_prime is: {}", min_prime);
    info!("thresh is: {}", thresh);
    info!("Fudge is: {}", fudge);

    let thresh = thresh - fudge;

    let mut roota = (_root2n / xmax).sqrt();
    if roota.is_divisible_2pow(1) {
        roota += 1;
    }

    info!("ROOTA is {}", roota);
    let roota: Integer = max(roota, Integer::from(3));
    info!("ROOTA is {}", roota);

    InitResult {
        roota,
        factorbase,
        tsqrt,
        xmax,
        tlog,
        thresh,
        min_prime,
    }
}

#[cfg(test)]
mod tests {
    use rug::Integer;

    use crate::check_is_divisor;
    use crate::serial_MPQS::mpqs;

    #[test]
    fn test_qs() {
        let n = "523022617466601111760007224100074291200000001"
            .parse::<Integer>()
            .unwrap();
        let ris = "14029308060317546154181".parse::<Integer>().unwrap();
        check_is_divisor(n, Some(ris));
    }

    #[test]
    fn test_qs_2() {
        let n = "9986801107".parse::<Integer>().unwrap();
        check_is_divisor(n.clone(), mpqs(&n));
    }
}
