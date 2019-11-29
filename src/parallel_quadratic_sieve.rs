use chashmap::CHashMap;
use crossbeam::queue::ArrayQueue;
use log::info;
use primal_sieve;
use rug::ops::Pow;
use rug::Integer;
use std::cmp::{max, min};
use std::collections::HashMap;

use std::sync::{mpsc::Sender, Arc, Mutex};

use crate::algebra;
use crate::tonelli_shank::tonelli_shank;

/// Memory shared: smooth and maybe partial
pub fn parallel_qs(n: &Integer) -> Option<Integer> {
    let _root2n: Integer = (n * Integer::from(2)).sqrt();

    info!("Isqrt is {}", _root2n);

    let bound: usize = (n.to_f64().log10().powf(2_f64) * 5_f64) as usize;

    info!("Bound is {}", bound);

    let factorbase: Vec<Integer> = {
        let mut v = vec![Integer::from(2)];
        v.extend(
            primal_sieve::Sieve::new(bound)
                .primes_from(3)
                .take_while(|x| x <= &bound)
                .map(|x| Integer::from(x))
                .filter(|x| n.legendre(x) == 1),
        );
        v
    };
    info!(
        "Largest prime used is {:?}",
        factorbase[factorbase.len() - 1]
    );
    info!("Factorbase len {:?}", factorbase.len());

    let mut tsqrt: Vec<Integer> = Vec::with_capacity(factorbase.len());
    let mut tlog: Vec<f64> = Vec::with_capacity(factorbase.len());

    for p in factorbase.iter() {
        //info!("Starting Tonelli Shank on {} - {}", n, p);
        let f1 = tonelli_shank(&n, &p);
        info!("Tonelli Shank on {} - {} = \t {}", n, p, f1);
        tsqrt.push(f1);
        tlog.push(p.to_f64().log10());
    }
    tsqrt[0] = Integer::new();

    info!("{:?}", tsqrt);
    info!("{:?}", tlog);

    let xmax: i64 = factorbase.len() as i64 * 60 * 4;
    let mval: Integer = (_root2n.clone() * xmax) >> 1;
    let thresh = mval.to_f64().log10() * 0.735;
    let min_prime = (thresh * 3_f64) as u64;
    let mut fudge: f64 = factorbase
        .iter()
        .take_while(|p| *p < &min_prime)
        .enumerate()
        .map(|(i, _)| tlog[i])
        .sum();
    fudge /= 4_f64;

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

    let smooths = ArrayQueue::new(factorbase.len() + 100);

    let (sender, receiver) = std::sync::mpsc::channel();
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
    algebra::algebra(factorbase.clone(), new_smooth.clone(), n)
}

fn thread_loop(
    n: Integer,
    factorbase: Vec<Integer>,
    smooths: Arc<ArrayQueue<(Integer, (Integer, Integer))>>,
    sender: Sender<()>,
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
        let b = tonelli_shank(&n, &my_roota);

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

        for low in (0 - xmax..xmax + 1).step_by(sievesize as usize + 1) {
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

            for i in 0..size + 1 {
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

    use super::parallel_qs;
    use crate::check_is_divisor;
    use rug::Integer;

    #[test]
    fn test_qs() {
        let n = "523022617466601111760007224100074291200000001"
            .parse::<Integer>()
            .unwrap();

        check_is_divisor(n.clone(), parallel_qs(&n));
    }

    #[test]
    fn test_qs_2() {
        let n = "9986801107".parse::<Integer>().unwrap();

        check_is_divisor(n.clone(), parallel_qs(&n));
    }

    #[test]
    fn test_qs_3() {
        let n = "2736300383840445596906210796102273501547527150973747"
            .parse::<Integer>()
            .unwrap();
        check_is_divisor(n.clone(), parallel_qs(&n));
    }

    #[test]
    #[ignore]
    fn test_qs_4() {
        let n = "676292275716558246502605230897191366469551764092181362779759"
            .parse::<Integer>()
            .unwrap();
        check_is_divisor(n.clone(), parallel_qs(&n));
    }
}
