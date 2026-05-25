use rand::{Rng, RngExt};
use rand_distr::StandardNormal;

#[derive(Clone, PartialEq)]
pub enum Ternary {
    Zero,
    PlusOne,
    MinusOne,
}

pub fn sample_ternary(n: usize, rng: &mut impl Rng, h: usize) -> Vec<Ternary> {
    let mut coefs = vec![Ternary::Zero; n];
    let mut n_set: usize = 0;

    while n_set < h {
        let idx = rng.random_range(0..n);
        if coefs[idx] != Ternary::Zero {
            continue;
        }

        coefs[idx] = match rng.random_bool(0.5) {
            true => Ternary::PlusOne,
            false => Ternary::MinusOne,
        };
        n_set += 1;
    }

    coefs
}

pub fn sample_gaussian(n: usize, rng: &mut impl Rng, sigma: f64) -> Vec<i64> {
    (0..n)
        .map(|_| (rng.sample::<f64, _>(StandardNormal) * sigma).round() as i64)
        .collect()
}
