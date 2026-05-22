use crate::ckks::CkksContext;
use crate::ckks::types::{Ciphertext, Plaintext, PublicKey};

pub struct Encryptor<'a> {
    ctx: &'a CkksContext,
    pk: &'a PublicKey,
}

impl<'a> Encryptor<'a> {
    pub fn new(ctx: &'a CkksContext, pk: &'a PublicKey) -> Self {
        Encryptor { ctx, pk }
    }

    pub fn encrypt(&self, _pt: &Plaintext) -> Ciphertext {
        todo!()
    }
}
