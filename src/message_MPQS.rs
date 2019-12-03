use std::cmp::min;
use std::collections::HashMap;

use log::info;
use rug::ops::Pow;
use rug::Integer;

use crate::algebra;
use crate::serial_MPQS::{initialize_qs, InitResult};
use crate::tonelli_shanks::tonelli_shanks;

/// Nothing Shared
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

    // Multi Producer - Single Consumer
    let (result_sender, result_receiver) = std::sync::mpsc::sync_channel(12);

    for _ in 0..num_cpus::get() {
        let z = n.clone();
        let factorbase = factorbase.clone();
        let result_sender = result_sender.clone();
        let tsqrt = tsqrt.clone();
        let tlog = tlog.clone();

        std::thread::spawn(move || {
            sieve_actor(
                z,
                factorbase,
                result_sender,
                tsqrt,
                tlog,
                xmax,
                min_prime,
                thresh,
            )
        });
    }

    let mut smooths = Vec::with_capacity(factorbase.len() + 100);
    let mut partials: HashMap<Integer, (Integer, (Integer, Integer))> = HashMap::new();

    while smooths.len() < factorbase.len() + 1 {
        let (mut sm, part, sender) = result_receiver.recv().unwrap();
        smooths.append(&mut sm);

        for (key, (pairv2, pairvals2)) in part {
            match partials.remove(&key) {
                Some((pairv, pairvals)) => {
                    smooths.push((
                        pairv * pairv2,
                        (pairvals2.0 * pairvals.0, pairvals.1 * pairvals2.1 * key),
                    ));
                }
                None => {
                    partials.insert(key, (pairv2, pairvals2));
                }
            }
        }
        roota.next_prime_mut();
        while n.legendre(&roota) != 1 {
            roota.next_prime_mut();
        }
        sender.send(roota.clone());
    }
    std::mem::drop(result_receiver);
    algebra::algebra(factorbase, smooths, n)
}

fn sieve_actor(
    n: Integer,
    factorbase: Vec<u64>,
    sender: std::sync::mpsc::SyncSender<(
        Vec<(Integer, (Integer, Integer))>,
        HashMap<Integer, (Integer, (Integer, Integer))>,
        std::sync::mpsc::SyncSender<Integer>,
    )>,
    tsqrt: Vec<Integer>,
    tlog: Vec<f64>,
    xmax: i64,
    min_prime: u64,
    thresh: f64,
) {
    let sievesize = 1_i64 << 15;
    let mut partials: HashMap<Integer, (Integer, (Integer, Integer))> = HashMap::new();
    let mut smooths: Vec<(Integer, (Integer, Integer))> = Vec::new();

    let (roota_sender, roota_receiver) = std::sync::mpsc::sync_channel(2);
    sender.send((Vec::new(), HashMap::new(), roota_sender.clone()));
    sender.send((Vec::new(), HashMap::new(), roota_sender.clone()));

    loop {
        let roota = roota_receiver.recv().unwrap();

        info!("Loop 1, roota: {}, n: {}", roota, n);
        let a = roota.clone().pow(2);
        let b = tonelli_shanks(&n, &roota);

        let int2: Integer = b.clone() * 2;
        let intermediate = int2.invert(&roota).expect("Inverse does not exist");
        let b = (-(b.clone() * &b - &n) * intermediate + &b) % &a;

        let c = (b.clone() * &b - &n) / &a;

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
                while sol1 < size_plus_1 || sol2 < size_plus_1 {
                    if sol1 < size_plus_1 {
                        S[sol1 as usize] += logp;
                        sol1 += p_i64;
                    }
                    if sol2 < size_plus_1 {
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
        if let Err(_) = sender.send((
            smooths.drain(..).collect(),
            partials.drain().collect(),
            roota_sender.clone(),
        )) {
            return;
        };
    }
}

#[cfg(test)]
mod tests {
    use rug::Integer;

    use crate::check_is_divisor;

    use super::mpqs;

    #[test]
    fn test_qs() {
        let n = "523022617466601111760007224100074291200000001"
            .parse::<Integer>()
            .unwrap();

        check_is_divisor(n.clone(), mpqs(&n));
    }

    #[test]
    fn test_qs_2() {
        let n = "9986801107".parse::<Integer>().unwrap();

        check_is_divisor(n.clone(), mpqs(&n));
    }

    #[test]
    #[ignore]
    fn test_qs_3() {
        let n = "2736300383840445596906210796102273501547527150973747"
            .parse::<Integer>()
            .unwrap();
        check_is_divisor(n.clone(), mpqs(&n));
    }

    #[test]
    #[ignore]
    fn test_qs_4() {
        let n = "676292275716558246502605230897191366469551764092181362779759"
            .parse::<Integer>()
            .unwrap();
        check_is_divisor(n.clone(), mpqs(&n));
    }
}
