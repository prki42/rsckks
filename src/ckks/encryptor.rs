use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, Plaintext, PublicKey};
use crate::sampling::{sample_gaussian, sample_ternary};

pub struct Encryptor<'a> {
    ctx: &'a CkksContext,
    pk: &'a PublicKey,
}

impl<'a> Encryptor<'a> {
    pub fn new(ctx: &'a CkksContext, pk: &'a PublicKey) -> Self {
        Encryptor { ctx, pk }
    }

    pub fn encrypt(&self, pt: &Plaintext) -> Ciphertext {
        let mut rng = rand::rng();

        let e0 = self
            .ctx
            .ring_q
            .poly_from_i64(&sample_gaussian(self.ctx.n(), &mut rng, 3.2));
        let e1 = self
            .ctx
            .ring_q
            .poly_from_i64(&sample_gaussian(self.ctx.n(), &mut rng, 3.2));

        let v = self.ctx.ring_q.poly_from_ternary(&sample_ternary(
            self.ctx.n(),
            &mut rng,
            self.ctx.n() / 2,
        ));

        // c0 = (v * b) + p + e0
        let mut c0 = self.ctx.ring_q.mul(&v, &self.pk.b);
        self.ctx.ring_q.add_inplace(&mut c0, &pt.data);
        self.ctx.ring_q.add_inplace(&mut c0, &e0);

        // c1 = (v * a) + e1
        let mut c1 = v;
        self.ctx.ring_q.mul_inplace(&mut c1, &self.pk.a);
        self.ctx.ring_q.add_inplace(&mut c1, &e1);

        Ciphertext {
            c0,
            c1,
            scale: pt.scale,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ckks::{keygen::KeyGenerator, test_utils::make_test_ctx};

    #[test]
    fn encrypt_test() {
        let ctx = make_test_ctx();
        let k = KeyGenerator::new(&ctx);
        let sk = k.secret_key();
        let pk = k.public_key(&sk);
        let e = Encryptor::new(&ctx, &pk);

        let mut rng = rand::rng();

        let pt = Plaintext {
            data: ctx.ring_q.sample_uniform(&mut rng),
            scale: 1.0,
        };

        e.encrypt(&pt);
    }
}
