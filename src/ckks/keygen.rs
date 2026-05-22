use crate::ckks::CkksContext;
use crate::ckks::types::{GaloisKeys, PublicKey, RelinKey, SecretKey};

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
        SecretKey {
            sk: self
                .ctx
                .ring_q
                .sample_ternary(&mut rng, self.ctx.ring_q.n() / 2),
        }
    }

    /// Generates a public key from the secret key
    pub fn public_key(&self, sk: &SecretKey) -> PublicKey {
        let ring_q = self.ctx.ring_q();
        let mut rng = rand::rng();

        let a = self.ctx.ring_q.sample_uniform(&mut rng);
        let e = self.ctx.ring_q.sample_gaussian(&mut rng, 3.2);

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
