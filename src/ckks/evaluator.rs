use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, GaloisKeys, Plaintext, RelinKey};
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
            level: a.level,
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
            level: a.level,
        })
    }

    pub fn mul(&self, _a: &Ciphertext, _b: &Ciphertext) -> Ciphertext {
        todo!()
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

    pub fn rescale(&self, _ct: &Ciphertext) -> Ciphertext {
        todo!()
    }

    pub fn rotate(&self, _ct: &Ciphertext, _steps: i32, _gk: &GaloisKeys) -> Ciphertext {
        todo!()
    }

    fn check_same_level(a: &Ciphertext, b: &Ciphertext) -> Result<(), EvalError> {
        if a.level != b.level {
            return Err(EvalError::MismatchedLevels(a.level, b.level));
        }
        Ok(())
    }

    fn check_same_scale(a: &Ciphertext, b: &Ciphertext) -> Result<(), EvalError> {
        if (a.scale - b.scale).abs() > 1e-6 {
            return Err(EvalError::MismatchedScales(a.scale, b.scale));
        }
        Ok(())
    }
}
