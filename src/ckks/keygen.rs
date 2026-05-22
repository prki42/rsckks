use rand::RngExt;
use rand_distr::StandardNormal;

use crate::ckks::CkksContext;
use crate::ckks::types::{GaloisKeys, PublicKey, RelinKey, SecretKey};
use crate::rns::{CoeffForm, NttForm, Poly};

pub struct KeyGenerator<'a> {
    ctx: &'a CkksContext,
}

impl<'a> KeyGenerator<'a> {
    pub fn new(ctx: &'a CkksContext) -> Self {
        KeyGenerator { ctx }
    }

    /// Samples a new secret key
    pub fn secret_key(&self) -> SecretKey {
        let mut rng = rand::rng();
        let n = self.ctx.n();

        // TODO: temp
        let h = n / 2;

        #[derive(Clone, PartialEq)]
        enum Ternary {
            Zero,
            PlusOne,
            MinusOne,
        }

        let mut coefs = vec![Ternary::Zero; n];
        let mut n_set: usize = 0;

        while n_set < h {
            let idx = rng.random_range(0..n);
            if coefs[idx] != Ternary::Zero {
                continue;
            }

            match rng.random_bool(0.5) {
                true => {
                    coefs[idx] = Ternary::PlusOne;
                }
                false => {
                    coefs[idx] = Ternary::MinusOne;
                }
            }
            n_set += 1;
        }

        let q_moduli_count = self.ctx.ring_q.num_moduli();
        let mut sk_coefs = vec![vec![0u64; n]; q_moduli_count];

        for coef_idx in 0..n {
            for (mod_idx, limb) in sk_coefs.iter_mut().enumerate() {
                let modulus = self.ctx.ring_q().modulus(mod_idx);
                limb[coef_idx] = match coefs[coef_idx] {
                    Ternary::PlusOne => 1,
                    Ternary::MinusOne => modulus - 1,
                    Ternary::Zero => 0,
                }
            }
        }

        SecretKey {
            sk: self.ctx.ring_q.ntt(Poly::<CoeffForm>::new(sk_coefs)),
        }
    }

    pub fn public_key(&self, sk: &SecretKey) -> PublicKey {
        let ring_q = self.ctx.ring_q();
        let mut rng = rand::rng();

        let a = {
            let a_limbs: Vec<Vec<u64>> = (0..ring_q.num_moduli())
                .map(|i| {
                    let qi = ring_q.modulus(i);
                    (0..self.ctx.n()).map(|_| rng.random_range(0..qi)).collect()
                })
                .collect();

            Poly::<NttForm>::new(a_limbs)
        };

        // TODO: temp, should be revised later
        let e = {
            let errors: Vec<i64> = (0..self.ctx.n())
                .map(|_| (rng.sample::<f64, _>(StandardNormal) * 3.2).round() as i64)
                .collect();

            let e_limbs: Vec<Vec<u64>> = (0..ring_q.num_moduli())
                .map(|i| {
                    let q = ring_q.modulus(i) as i64;
                    errors
                        .iter()
                        .map(|&e| ((e % q) + q) as u64 % q as u64)
                        .collect()
                })
                .collect();

            ring_q.ntt(Poly::<CoeffForm>::new(e_limbs))
        };

        // b = e - a*s
        let mut b = e;
        ring_q.sub_inplace(&mut b, &ring_q.mul(&a, &sk.sk));

        PublicKey { a, b }
    }

    pub fn relin_key(&self, _sk: &SecretKey) -> RelinKey {
        todo!()
    }

    pub fn galois_keys(&self, _sk: &SecretKey, _rotations: &[i32]) -> GaloisKeys {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::test_utils::make_test_ctx;

    // TODO: proper tests

    #[test]
    fn secretkey_gen() {
        let ctx = make_test_ctx();
        let k = KeyGenerator::new(&ctx);
        k.secret_key();
    }

    #[test]
    fn publickey_gen() {
        let ctx = make_test_ctx();
        let k = KeyGenerator::new(&ctx);
        k.public_key(&k.secret_key());
    }
}
