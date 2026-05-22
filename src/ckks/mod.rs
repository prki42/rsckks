pub mod decryptor;
pub mod encoder;
pub mod encryptor;
pub mod evaluator;
pub mod keygen;
pub mod params;
pub mod types;

use thiserror::Error;

use crate::rns::{RnsRing, RnsRingError};

#[derive(Debug)]
pub struct CkksContext {
    ring_q: RnsRing,
    ring_p: RnsRing,
    scale: f64,
}

#[derive(Error, Debug)]
pub enum CkksContextErr {
    #[error(transparent)]
    RnsRing(#[from] RnsRingError),
}

impl CkksContext {
    pub fn new(
        q_moduli: &[u64],
        p_moduli: &[u64],
        n: usize,
        scale: f64,
    ) -> Result<Self, CkksContextErr> {
        Ok(CkksContext {
            ring_q: RnsRing::new(q_moduli, n)?,
            ring_p: RnsRing::new(p_moduli, n)?,
            scale,
        })
    }

    pub fn ring_q(&self) -> &RnsRing {
        &self.ring_q
    }

    pub fn ring_p(&self) -> &RnsRing {
        &self.ring_p
    }

    pub fn n(&self) -> usize {
        self.ring_q.n()
    }

    pub fn num_levels(&self) -> usize {
        self.ring_q.num_moduli()
    }

    pub fn scale(&self) -> f64 {
        self.scale
    }
}

#[cfg(test)]
pub(crate) mod test_utils {
    use super::*;

    // TODO: better test ctx, also should replace other adhoc ctx creations in tests

    pub fn make_test_ctx() -> CkksContext {
        CkksContext::new(&[998244353, 985661441, 754974721], &[469762049], 256, 64.0).unwrap()
    }
}
