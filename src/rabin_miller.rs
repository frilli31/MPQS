use rug::Integer;

pub fn is_rabin_miller_prime(n: &Integer) -> bool {
    if *n == 2 || *n == 3 {
        return true;
    } else if *n <= 1 || (n.clone() % 2) == 0 {
        return false;
    }

    let mut s: Integer = Integer::from(0);
    let mut d: Integer = n.clone() - 1;
    while (d.clone() % 2) == 0 {
        d >>= 1;
        s = s + 1;
    }

    for a in primes_to_test(&n)
        .iter()
        .map(|a| Integer::from(*a))
        .filter(|a| a.clone() + 1 < *n)
    {
        if try_composite(&a, &d, n, &s) {
            return false;
        }
    }
    true
}

fn try_composite(a: &Integer, d: &Integer, n: &Integer, s: &Integer) -> bool {
    let two = Integer::from(2);
    let mut a_to_n = a.clone().pow_mod(d, n).unwrap();

    if a_to_n == 1 {
        false
    } else {
        let mut i = Integer::new();

        while i < s.clone() - 1 {
            if a_to_n == n.clone() - 1 {
                return false;
            }
            a_to_n = a_to_n.pow_mod(&two, n).unwrap();
            i = i + 1;
        }
        a_to_n != n.clone() - 1
    }
}

fn primes_to_test(n: &Integer) -> Vec<u64> {
    let n128 = n.to_u128();

    match n128 {
        Some(n) => match n {
            0..=1_373_653 => [2, 3].to_vec(),
            1_373_653..=25_326_001 => [2, 3, 5].to_vec(),
            25_326_001..=118670087467 => [2, 3, 5, 7].to_vec(),
            118670087467..=2152302898747 => [2, 3, 5, 7, 11].to_vec(),
            2152302898747..=3474749660383 => [2, 3, 5, 7, 11, 13].to_vec(),
            3474749660383..=341550071728321 => [2, 3, 5, 7, 11, 13, 17].to_vec(),
            341550071728321..=3825123056546413051 => [2, 3, 5, 7, 11, 13, 17, 19, 23].to_vec(),
            3825123056546413051..=318665857834031151167461 => {
                [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37].to_vec()
            }
            318665857834031151167461..=3317044064679887385961981 => {
                [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41].to_vec()
            }
            _ => [
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
            ]
            .to_vec(),
        },
        None => [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
        ]
        .to_vec(),
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_is_prime() {
        let pairs: Vec<(u64, bool)> = vec![
            (1, false),
            (2, true),
            (3, true),
            (4, false),
            (5, true),
            (123123423467, false),
            (4373, true),
            (1048576, false),
        ];

        for (n, r) in pairs {
            assert_eq!(is_rabin_miller_prime(&Integer::from(n)), r);
        }
    }
}