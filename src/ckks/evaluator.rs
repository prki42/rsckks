use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, Plaintext, RelinKey};
use thiserror::Error;

pub struct Evaluator<'a> {
    ctx: &'a CkksContext,
}

#[derive(Error, Debug)]
pub enum EvalError {
    #[error("levels do not match, {0} != {1}")]
    MismatchedLevels(usize, usize),

    #[error("scales do not match, {0} != {1}")]
    MismatchedScales(f64, f64),
}

impl<'a> Evaluator<'a> {
    pub fn new(ctx: &'a CkksContext) -> Self {
        Evaluator { ctx }
    }

    pub fn add(&self, a: &Ciphertext, b: &Ciphertext) -> Result<Ciphertext, EvalError> {
        Self::check_same_level(a, b)?;
        Self::check_same_scale(a, b)?;

        let c0 = self.ctx.ring_q.add(&a.c0, &b.c0);
        let c1 = self.ctx.ring_q.add(&a.c1, &b.c1);

        Ok(Ciphertext {
            c0,
            c1,
            scale: a.scale,
        })
    }

    pub fn sub(&self, a: &Ciphertext, b: &Ciphertext) -> Result<Ciphertext, EvalError> {
        Self::check_same_level(a, b)?;
        Self::check_same_scale(a, b)?;

        let c0 = self.ctx.ring_q.sub(&a.c0, &b.c0);
        let c1 = self.ctx.ring_q.sub(&a.c1, &b.c1);

        Ok(Ciphertext {
            c0,
            c1,
            scale: a.scale,
        })
    }

    pub fn rescale(&self, _ct: &Ciphertext) -> Ciphertext {
        todo!()
    }

    pub fn mul(
        &self,
        a: &Ciphertext,
        b: &Ciphertext,
        evk: &RelinKey,
    ) -> Result<Ciphertext, EvalError> {
        Self::check_same_level(a, b)?;

        let d0 = self.ctx.ring_q.mul(&a.c0, &b.c0);

        let mut d1 = self.ctx.ring_q.mul(&a.c0, &b.c1);
        self.ctx
            .ring_q
            .add_inplace(&mut d1, &self.ctx.ring_q.mul(&a.c1, &b.c0));

        let mut d2_ntt_q = self.ctx.ring_q.mul(&a.c1, &b.c1);

        let d2_coef_q = self.ctx.ring_q.intt(d2_ntt_q.clone());
        let d2_coef_p = self.ctx.ring_p.mod_change(
            &d2_coef_q,
            &self.ctx.ring_q,
            &self.ctx.cross_ring.q_hat_mod_p,
        );

        let mut d2_ntt_p = self.ctx.ring_p.ntt(d2_coef_p);

        let mut c0_q = self
            .ctx
            .ring_q
            .intt(self.ctx.ring_q.mul(&d2_ntt_q, &evk.mod_q.0));
        self.ctx.ring_q.mul_inplace(&mut d2_ntt_q, &evk.mod_q.1);
        let mut c1_q = self.ctx.ring_q.intt(d2_ntt_q);

        let c0_p = self
            .ctx
            .ring_p
            .intt(self.ctx.ring_p.mul(&d2_ntt_p, &evk.mod_p.0));
        self.ctx.ring_p.mul_inplace(&mut d2_ntt_p, &evk.mod_p.1);
        let c1_p = self.ctx.ring_p.intt(d2_ntt_p);

        self.ctx.ring_q.mod_scale_down(
            &mut c0_q,
            &c0_p,
            &self.ctx.ring_p,
            &self.ctx.cross_ring.p_hat_mod_q,
            &self.ctx.cross_ring.p_inv_mod_q,
        );
        self.ctx.ring_q.mod_scale_down(
            &mut c1_q,
            &c1_p,
            &self.ctx.ring_p,
            &self.ctx.cross_ring.p_hat_mod_q,
            &self.ctx.cross_ring.p_inv_mod_q,
        );

        self.ctx
            .ring_q
            .add_inplace(&mut c0_q, &self.ctx.ring_q.intt(d0));
        self.ctx
            .ring_q
            .add_inplace(&mut c1_q, &self.ctx.ring_q.intt(d1));

        let ql = self.ctx.ring_q.modulus(c0_q.limbs.len() - 1) as f64;

        self.ctx.ring_q.rescale(&mut c0_q);
        self.ctx.ring_q.rescale(&mut c1_q);

        let c0 = self.ctx.ring_q.ntt(c0_q);
        let c1 = self.ctx.ring_q.ntt(c1_q);

        Ok(Ciphertext {
            c0,
            c1,
            scale: a.scale * b.scale / ql,
        })
    }

    pub fn mul_plain(&self, _ct: &Ciphertext, _pt: &Plaintext) -> Ciphertext {
        todo!()
    }

    pub fn add_plain(&self, _ct: &Ciphertext, _pt: &Plaintext) -> Ciphertext {
        todo!()
    }

    pub fn relinearize(&self, _ct: &Ciphertext, _rlk: &RelinKey) -> Ciphertext {
        todo!()
    }

    fn check_same_level(a: &Ciphertext, b: &Ciphertext) -> Result<(), EvalError> {
        if a.level() != b.level() {
            return Err(EvalError::MismatchedLevels(a.level(), b.level()));
        }
        Ok(())
    }

    // TODO: different check (if at all??)
    fn check_same_scale(a: &Ciphertext, b: &Ciphertext) -> Result<(), EvalError> {
        if (a.scale - b.scale).abs() > 1e-6 {
            return Err(EvalError::MismatchedScales(a.scale, b.scale));
        }
        Ok(())
    }
}
