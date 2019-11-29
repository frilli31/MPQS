use log::info;
use rug::Integer;

/// Tonelli Shank Algorithm
pub fn tonelli_shank(n: &Integer, p: &Integer) -> Integer {
    if n.legendre(p) != 1 || *n == 0 {
        Integer::new()
    } else if *p == 2 {
        p.clone()
    } else if p.mod_u(4) == 3 {
        n.clone().pow_mod(&((p.clone() + 1) / 4), p).unwrap()
    } else {
        let (mut s, mut e): (Integer, _) = ((p.clone() - 1), 0);
        while s.is_divisible_2pow(1) {
            s /= 2;
            e += 1;
        }

        let mut z: Integer = 2.into();

        while z.legendre(p) != -1 {
            z += 1;
        }

        let mut x = n.clone().pow_mod(&((s.clone() + 1) / 2), p).unwrap();
        let mut b = n.clone().pow_mod(&s, p).unwrap();
        let mut g = z.pow_mod(&s, p).unwrap();
        let mut r = e;

        info!("x = {} b={} g={} r={}", x, b, g, r);

        let TWO = Integer::from(2);
        loop {
            let mut t = b.clone();

            let mut m = 0;

            while m < r && t != 1 {
                t = t.pow_mod(&TWO, p).unwrap();

                m += 1;
            }
            if m == 0 {
                return x;
            }

            let gs = g
                .clone()
                .pow_mod(&Integer::from(Integer::u_pow_u(2, r - m - 1)), p)
                .unwrap();

            g = gs.clone() * &gs % p;
            x = x * gs % p;
            b = b * &g % p;
            r = m;

            info!("g={} x={} b={} r={}", g, x, b, r);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tonelli_shank() {
        let (n, p) = (Integer::from(23479349), Integer::from(23));
        assert_eq!(tonelli_shank(&n, &p), Integer::from(12));
    }
}
