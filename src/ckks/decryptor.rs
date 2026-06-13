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
        let mut p = self.ctx.ring_q.mul(&ct.c1, &self.sk.sk_q);
        self.ctx.ring_q.add_inplace(&mut p, &ct.c0);

        Plaintext {
            data: p,
            scale: ct.scale,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::test_utils::*;
    use proptest::prelude::*;

    proptest! {
        #![proptest_config(no_persist(30))]

        #[test]
        fn encrypt_decrypt_roundtrip(
            (tc, values) in with_test_env(|tc| complex_vec(1.0, tc.slots()))
        ) {
            let ct = tc.encrypt(&values);
            let result = tc.decrypt_decode(&ct);
            // TODO: change to a proper bound
            assert_slots_approx(&values, &result, 1.0)?;
        }
    }
}
