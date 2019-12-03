use rug::Integer;

/// Tonelli Shanks Algorithm
pub fn tonelli_shanks(n: &Integer, p: &Integer) -> Integer {
    if n.legendre(p) != 1 || *n == 0 {
        Integer::new()
    } else if *p == 2 {
        p.clone()
    } else if p.is_congruent_u(3, 4) {
        n.clone().pow_mod(&((p.clone() + 1) / 4), p).unwrap()
    } else {
        let mut s: Integer = p.clone() - 1;
        let e = s.find_one(0).unwrap();
        s >>= e;

        let mut z: Integer = Integer::from(2);
        while z.legendre(p) != -1 {
            z += 1;
        }

        let mut b = n.clone().pow_mod(&s, p).unwrap();
        let mut g = z.pow_mod(&s, p).unwrap();
        let mut x = n.clone().pow_mod(&((s + 1) / 2), p).unwrap();
        let mut r = e;

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
                .pow_mod(&(Integer::from(1) << (r - m - 1)), p)
                .unwrap();

            g = gs.clone().square() % p;
            x = x * gs % p;
            b = b * &g % p;
            r = m;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tonelli_shank() {
        let (n, p) = (Integer::from(23479349), Integer::from(23));
        assert_eq!(tonelli_shanks(&n, &p), Integer::from(12));
    }
}
