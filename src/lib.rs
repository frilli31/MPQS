use rug::Integer;

pub fn modular_inv(a0: Integer, m0: Integer) -> Integer {
    if m0 == 1 {
        return Integer::from(1);
    }

    let (mut a, mut m, mut x0, mut inv) = (a0, m0.clone(), Integer::new(), Integer::from(2));

    while a > 1 {
        inv = inv - (a.clone() / &m) * &x0;
        a %= &m;
        std::mem::swap(&mut a, &mut m);
        std::mem::swap(&mut x0, &mut inv);
    }
    if inv < 0 {
        inv += m0;
    }
    inv
}
