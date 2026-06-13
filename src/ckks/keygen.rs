use crate::ckks::CkksContext;
use crate::ckks::types::{PublicKey, RelinKey, SecretKey};
use crate::sampling::{sample_gaussian, sample_ternary};

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
        let s = sample_ternary(self.ctx.n(), &mut rng, self.ctx.n() / 2);

        SecretKey {
            sk_q: self.ctx.ring_q.poly_from_ternary(&s),
            sk_p: self.ctx.ring_p.poly_from_ternary(&s),
        }
    }

    /// Generates a public key from the secret key
    pub fn public_key(&self, sk: &SecretKey) -> PublicKey {
        let ring_q = self.ctx.ring_q();
        let mut rng = rand::rng();

        let a = self.ctx.ring_q.sample_uniform(&mut rng);
        let e = self
            .ctx
            .ring_q
            .poly_from_i64(&sample_gaussian(self.ctx.n(), &mut rng, 3.2));

        // b = e - a*s
        let mut b = e;
        ring_q.sub_inplace(&mut b, &ring_q.mul(&a, &sk.sk_q));

        PublicKey { a, b }
    }

    pub fn relin_key(&self, sk: &SecretKey) -> RelinKey {
        let mut rng = rand::rng();
        let mut s2_q = self.ctx.ring_q.mul(&sk.sk_q, &sk.sk_q);

        self.ctx
            .ring_q
            .mul_const(&mut s2_q, &self.ctx.cross_ring.p_mod_q);

        let a_q = self.ctx.ring_q.sample_uniform(&mut rng);
        let a_p = self.ctx.ring_p.sample_uniform(&mut rng);

        let err = sample_gaussian(self.ctx.n(), &mut rng, 3.2);

        let mut err_q = self.ctx.ring_q.poly_from_i64(&err);
        let mut err_p = self.ctx.ring_p.poly_from_i64(&err);

        self.ctx
            .ring_p
            .sub_inplace(&mut err_p, &self.ctx.ring_p.mul(&a_p, &sk.sk_p));
        let b_p = err_p;

        self.ctx
            .ring_q
            .sub_inplace(&mut err_q, &self.ctx.ring_q.mul(&a_q, &sk.sk_q));
        self.ctx.ring_q.add_inplace(&mut err_q, &s2_q);
        let b_q = err_q;

        RelinKey {
            mod_q: (b_q, a_q),
            mod_p: (b_p, a_p),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::CKKS_TEST_ENVS;

    #[test]
    fn secretkey_gen() {
        let tc = &CKKS_TEST_ENVS[0];
        let k = KeyGenerator::new(&tc.ctx);
        k.secret_key();
    }

    #[test]
    fn publickey_gen() {
        let tc = &CKKS_TEST_ENVS[0];
        let k = KeyGenerator::new(&tc.ctx);
        k.public_key(&k.secret_key());
    }

    #[test]
    fn relinkey_gen() {
        let tc = &CKKS_TEST_ENVS[0];
        let k = KeyGenerator::new(&tc.ctx);
        let sk = k.secret_key();
        k.relin_key(&sk);
    }
}
