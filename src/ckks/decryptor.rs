use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, Plaintext, SecretKey};

pub struct Decryptor<'a> {
    ctx: &'a CkksContext,
    sk: &'a SecretKey,
}

impl<'a> Decryptor<'a> {
    pub fn new(ctx: &'a CkksContext, sk: &'a SecretKey) -> Self {
        Decryptor { ctx, sk }
    }

    pub fn decrypt(&self, ct: &Ciphertext) -> Plaintext {
        let mut p = self.ctx.ring_q.mul(&ct.c1, &self.sk.sk);
        self.ctx.ring_q.add_inplace(&mut p, &ct.c0);

        Plaintext {
            data: p,
            scale: ct.scale,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::{
        encoder::Encoder,
        encryptor::Encryptor,
        keygen::KeyGenerator,
        params::{CkksParams, gen_ckks_context},
        test_utils::make_test_ctx,
    };
    use num_complex::Complex64;
    use proptest::prelude::*;

    #[test]
    fn decrypt_test() {
        let ctx = make_test_ctx();
        let k = KeyGenerator::new(&ctx);
        let sk = k.secret_key();
        let e = Decryptor::new(&ctx, &sk);

        let mut rng = rand::rng();

        let ct = Ciphertext {
            c0: ctx.ring_q.sample_uniform(&mut rng),
            c1: ctx.ring_q.sample_uniform(&mut rng),
            level: 3,
            scale: 10.0,
        };

        e.decrypt(&ct);
    }

    fn complex_vec(max_val: f64, max_len: usize) -> impl Strategy<Value = Vec<Complex64>> {
        proptest::collection::vec((-max_val..max_val, -max_val..max_val), 1..=max_len).prop_map(
            |v| {
                v.into_iter()
                    .map(|(re, im)| Complex64::new(re, im))
                    .collect()
            },
        )
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(30))]

        #[test]
        fn encrypt_decrypt_roundtrip(values in complex_vec(10.0, 8)) {

            // TODO: use make_test_ctx when you figure out proper examples and bounds
            let ctx = gen_ckks_context(&CkksParams {
                first_mod_size: 61,
                scaling_size: 55,
                mul_depth: 15,
                ring_dim: 16,
            })
            .unwrap();

            let keygen = KeyGenerator::new(&ctx);
            let sk = keygen.secret_key();
            let pk = keygen.public_key(&sk);
            let enc = Encryptor::new(&ctx, &pk);
            let dec = Decryptor::new(&ctx, &sk);
            let encoder = Encoder::new(&ctx);

            let pt = encoder.encode(&values, &ctx);
            let ct = enc.encrypt(&pt);
            let pt_dec = dec.decrypt(&ct);
            let result = encoder.decode(&pt_dec, &ctx);

            for (i, (&expected, &actual)) in values.iter().zip(result.iter()).enumerate() {
                let diff = (expected - actual).norm();
                prop_assert!(
                    // TODO: change to a proper bound
                    diff < 1.0,
                    "slot {i}: expected {expected}, got {actual}, diff {diff}"
                );
            }
        }
    }
}
