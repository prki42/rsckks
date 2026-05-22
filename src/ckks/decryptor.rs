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

    pub fn decrypt(&self, _ct: &Ciphertext) -> Plaintext {
        todo!()
    }
}
