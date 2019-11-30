use std::cmp::min;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use chashmap::CHashMap;
use crossbeam::queue::ArrayQueue;
use log::info;
use rug::ops::Pow;
use rug::Integer;

use crate::algebra;
use crate::serial_MPQS::{initialize_qs, InitResult};
use crate::tonelli_shanks::tonelli_shanks;
use std::sync::mpsc::SyncSender;

pub fn mpqs(n: &Integer) -> Option<Integer> {
    let InitResult {
        roota,
        factorbase,
        tsqrt,
        xmax,
        tlog,
        thresh,
        min_prime,
    } = initialize_qs(n);

    let smooths = ArrayQueue::new(factorbase.len() + 100);

    let (sender, receiver) = std::sync::mpsc::sync_channel(num_cpus::get());
    let roota = Arc::new(Mutex::new(roota));

    let arc_smooths = Arc::new(smooths);
    let partials = Arc::new(CHashMap::default());

    for _ in 0..num_cpus::get() {
        let z = n.clone();
        let factorbase = factorbase.clone();
        let sender = sender.clone();
        let tsqrt = tsqrt.clone();
        let tlog = tlog.clone();
        let arc_smooths = arc_smooths.clone();
        let roota = roota.clone();
        let partials = partials.clone();

        std::thread::spawn(move || {
            thread_loop(
                z,
                factorbase,
                arc_smooths,
                sender,
                roota,
                tsqrt,
                tlog,
                xmax,
                min_prime,
                thresh,
                partials,
            )
        });
    }

    let _ = receiver.recv();

    let mut new_smooth: Vec<_> = Vec::with_capacity(arc_smooths.len());
    for _ in 0..arc_smooths.len() {
        new_smooth.push(arc_smooths.pop().unwrap());
    }
    algebra::algebra(factorbase, new_smooth, n)
}

fn thread_loop(
    n: Integer,
    factorbase: Vec<Integer>,
    smooths: Arc<ArrayQueue<(Integer, (Integer, Integer))>>,
    sender: SyncSender<()>,
    roota: Arc<Mutex<Integer>>,
    tsqrt: Vec<Integer>,
    tlog: Vec<f64>,
    xmax: i64,
    min_prime: u64,
    thresh: f64,
    partials: Arc<CHashMap<Integer, (Integer, (Integer, Integer))>>,
) {
    let sievesize = 1_i64 << 15;

    loop {
        let my_roota: Integer = {
            let mut aq_roota = roota.lock().unwrap();
            aq_roota.next_prime_mut();
            while n.legendre(&aq_roota) != 1 {
                aq_roota.next_prime_mut();
            }
            aq_roota.clone()
        };
        info!("Loop 1, roota: {}, n: {}", my_roota, n);
        let a = my_roota.clone().pow(2);
        let b = tonelli_shanks(&n, &my_roota);

        let int2: Integer = b.clone() * 2;
        let intermediate = int2.invert(&my_roota).expect("Inverse does not exist");
        let b = (-(b.clone() * &b - &n) * intermediate + &b) % &a;

        let c = (b.clone() * &b - &n) / &a;

        info!("a={} \t b={} \t c={}", a, b, c);

        let mut s1: HashMap<Integer, Integer> = HashMap::new();
        let mut s2: HashMap<Integer, Integer> = HashMap::new();

        for (i, p) in factorbase.iter().enumerate() {
            let p_minus_2 = p.clone() - 2;
            let ainv = a.clone().pow_mod(&p_minus_2, p).unwrap();
            let mut sol1 = (tsqrt[i].clone() - &b) * &ainv % p;
            let mut sol2 = (-tsqrt[i].clone() - &b) * &ainv % p;
            sol1 -= ((sol1.clone() + xmax) / p) * p;
            sol2 -= ((sol2.clone() + xmax) / p) * p;

            s1.insert(p.clone(), sol1 + xmax);
            s2.insert(p.clone(), sol2 + xmax);
        }

        for low in (-xmax..xmax + 1).step_by(sievesize as usize + 1) {
            let high = min(xmax, low + sievesize);
            let size = high - low;
            let size_plus_1 = size + 1;

            let mut S = vec![0_f64; size_plus_1 as usize];

            for (i, p_i) in factorbase.iter().enumerate() {
                if *p_i < min_prime {
                    continue;
                }
                let p = p_i.to_i64().unwrap();
                let mut sol1 = s1[p_i].to_i64().unwrap();
                let mut sol2 = s2[p_i].to_i64().unwrap();
                let logp = tlog[i];

                while sol1 <= size || sol2 <= size {
                    if sol1 <= size {
                        S[sol1 as usize] += logp;
                        sol1 += p;
                    }
                    if sol2 <= size {
                        S[sol2 as usize] += logp;
                        sol2 += p;
                    }
                }
                s1.insert(p_i.clone(), Integer::from(sol1 - size_plus_1));
                s2.insert(p_i.clone(), Integer::from(sol2 - size_plus_1));
            }

            for i in 0..size_plus_1 {
                if S[i as usize] > thresh {
                    let x = i + low;
                    let tofact: Integer = a.clone() * x.pow(2) + b.clone() * x * 2 + &c;
                    let mut nf = tofact.clone().abs();

                    for p in factorbase.iter() {
                        while nf.is_divisible(p) {
                            nf /= p;
                        }
                    }

                    if nf == 1 {
                        smooths.push((a.clone() * x + &b, (tofact, my_roota.clone())));
                    } else {
                        match partials.remove(&nf) {
                            Some((pairv, pairvals)) => {
                                smooths.push((
                                    pairv * (a.clone() * x + &b),
                                    (tofact * pairvals.0, pairvals.1 * &my_roota * nf),
                                ));
                            }
                            None => {
                                partials
                                    .insert(nf, (a.clone() * x + &b, (tofact, my_roota.clone())));
                            }
                        }
                    }
                }
            }
        }

        if smooths.len() > factorbase.len() {
            sender.send(());
            return;
        }
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
