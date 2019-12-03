use rug::Integer;
use rug::ops::Pow;

pub fn algebra(
    factorbase: &[u64],
    smooths: &[(Integer, (Integer, Integer))],
    settings: &Integer,
) -> Option<Integer> {
    let n = settings;
    let mut m_vector: Vec<Integer> = smooths
        .iter()
        .map(|(_, (el, _))| create_vector(el, &factorbase))
        .collect();

    let factorbase_new = {
        let mut temp = vec![Integer::from(-1)];
        factorbase.iter().for_each(|i| temp.push(Integer::from(*i)));
        temp
    };

    let mut h_vector: Vec<Integer> = (0..m_vector.len())
        .map(|i| Integer::from(1) << i as u32)
        .collect();

    reduce_row_echelon_form(&mut m_vector, &mut h_vector, factorbase_new.len());

    let nul_cols: Vec<Integer> = h_vector
        .into_iter()
        .enumerate()
        .filter(|(i, _)| m_vector[*i] == 0)
        .map(|(_, e)| e)
        .collect();

    for nc in nul_cols {
        let mut lhs = Integer::from(1);
        let mut rhs = vec![Integer::new(); factorbase_new.len()];
        let mut rhspr = Integer::from(1);

        for (index, (lh, (rh, ra))) in smooths.iter().enumerate() {
            if (Integer::from(1) << index as u32) & &nc > 0 {
                lhs *= lh;
                rhspr *= ra;
                if *rh < 0 {
                    rhs[0] += 1;
                }
                let mut rh = rh.clone();
                for j in 1..factorbase_new.len() {
                    while rh.is_divisible(&factorbase_new[j]) {
                        rh /= &factorbase_new[j];
                        rhs[j] += 1;
                    }
                }
            }
        }
        for (j, factor) in factorbase_new.iter().enumerate() {
            rhspr *= factor.clone().pow(rhs[j].to_u32().unwrap() >> 1);
        }
        let g = Integer::from(rhspr - lhs).gcd(n);
        if g != 1 && g != *n {
            return Some(g);
        }
    }
    None
}

fn create_vector(n: &Integer, factor_base: &[u64]) -> Integer {
    let mut n = n.clone();
    let mut a = Integer::new();
    let lg = factor_base.len() - 1;
    if n < 0 {
        a |= Integer::from(2) << lg as u32;
        n = -n;
    }
    for (i, p) in factor_base.iter().enumerate() {
        if n.clone() % p == 0 {
            let mut c = 0;
            while n.clone() % p == 0 {
                n /= p;
                c += 1;
            }
            if c & 1 > 0 {
                a |= Integer::from(1) << (lg - i) as u32;
            }
        }
    }
    a
}

fn reduce_row_echelon_form(m: &mut Vec<Integer>, h: &mut Vec<Integer>, column_count: usize) {
    if m.is_empty() {
        return;
    }
    let mut lead = 0;
    let row_count = m.len();
    for r in 0..row_count {
        if lead >= column_count {
            return;
        }
        let mut i = r;
        while (Integer::from(1) << lead as u32) & &m[i] == 0 {
            i += 1;
            if i == row_count {
                i = r;
                lead += 1;
                if column_count == lead {
                    return;
                }
            }
        }
        m.swap(i, r);
        h.swap(i, r);

        for index in 0..row_count {
            if index != r && (Integer::from(1) << lead as u32) & &m[index] != 0 {
                let temp1 = m[r].clone();
                let temp2 = h[r].clone();
                m[index] ^= temp1;
                h[index] ^= temp2;
            }
        }
        lead += 1;
    }
}
