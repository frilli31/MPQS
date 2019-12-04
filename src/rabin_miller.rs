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
        s += 1;
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
            i += 1;
        }
        a_to_n != n.clone() - 1
    }
}

fn primes_to_test(n: &Integer) -> Vec<u64> {
    let n128 = n.to_u128();

    match n128 {
        Some(n) => match n {
            0..=1_373_653 => [2, 3].to_vec(),
            1_373_654..=25_326_001 => [2, 3, 5].to_vec(),
            25_326_002..=118_670_087_467 => [2, 3, 5, 7].to_vec(),
            118_670_087_468..=2_152_302_898_747 => [2, 3, 5, 7, 11].to_vec(),
            2_152_302_898_748..=3_474_749_660_383 => [2, 3, 5, 7, 11, 13].to_vec(),
            3_474_749_660_384..=341_550_071_728_321 => [2, 3, 5, 7, 11, 13, 17].to_vec(),
            341_550_071_728_322..=3_825_123_056_546_413_051 => {
                [2, 3, 5, 7, 11, 13, 17, 19, 23].to_vec()
            }
            3_825_123_056_546_413_052..=318_665_857_834_031_151_167_461 => {
                [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37].to_vec()
            }
            318_665_857_834_031_151_167_462..=3_317_044_064_679_887_385_961_981 => {
                [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41].to_vec()
            }
            _ => [
                2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
                83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
            ]
            .to_vec(),
        },
        None => [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
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
