use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, GaloisKeys, Plaintext, PublicKey, RelinKey, SecretKey};

pub struct KeyGenerator<'a> {
    ctx: &'a CkksContext,
}

impl<'a> KeyGenerator<'a> {
    pub fn new(ctx: &'a CkksContext) -> Self {
        KeyGenerator { ctx }
    }

    pub fn secret_key(&self) -> SecretKey {
        todo!()
    }

    pub fn public_key(&self, _sk: &SecretKey) -> PublicKey {
        todo!()
    }

    pub fn relin_key(&self, _sk: &SecretKey) -> RelinKey {
        todo!()
    }

    pub fn galois_keys(&self, _sk: &SecretKey, _rotations: &[i32]) -> GaloisKeys {
        todo!()
    }

    pub fn encrypt(&self, _pk: &PublicKey, _pt: &Plaintext) -> Ciphertext {
        todo!()
    }

    pub fn decrypt(&self, _sk: &SecretKey, _ct: &Ciphertext) -> Plaintext {
        todo!()
    }
}
