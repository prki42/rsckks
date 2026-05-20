pub mod encoder;
pub mod evaluator;
pub mod keygen;
pub mod params;
pub mod types;

use crate::rns::RnsRing;

pub struct CkksContext {
    ring_q: RnsRing,
    ring_p: RnsRing,
}

impl CkksContext {
    pub fn new(q_moduli: &[u64], p_moduli: &[u64], n: usize) -> Self {
        CkksContext {
            ring_q: RnsRing::new(q_moduli, n),
            ring_p: RnsRing::new(p_moduli, n),
        }
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
}
