use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, GaloisKeys, Plaintext, RelinKey};

pub struct Evaluator<'a> {
    ctx: &'a CkksContext,
}

impl<'a> Evaluator<'a> {
    pub fn new(ctx: &'a CkksContext) -> Self {
        Evaluator { ctx }
    }

    pub fn add(&self, _a: &Ciphertext, _b: &Ciphertext) -> Ciphertext {
        todo!()
    }

    pub fn sub(&self, _a: &Ciphertext, _b: &Ciphertext) -> Ciphertext {
        todo!()
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
}
