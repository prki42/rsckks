use crate::rns::{NttForm, Poly};

#[derive(Debug)]
pub struct SecretKey {
    pub(crate) sk_q: Poly<NttForm>,
    pub(crate) sk_p: Poly<NttForm>,
}

pub struct PublicKey {
    pub(crate) a: Poly<NttForm>,
    pub(crate) b: Poly<NttForm>,
}

pub struct RelinKey {
    pub(crate) mod_q: (Poly<NttForm>, Poly<NttForm>),
    pub(crate) mod_p: (Poly<NttForm>, Poly<NttForm>),
}

pub struct Plaintext {
    pub(crate) data: Poly<NttForm>,
    pub(crate) scale: f64,
}

pub struct Ciphertext {
    pub(crate) c0: Poly<NttForm>,
    pub(crate) c1: Poly<NttForm>,
    pub(crate) scale: f64,
}

impl Ciphertext {
    pub fn scale(&self) -> f64 {
        self.scale
    }

    pub fn level(&self) -> usize {
        self.c0.limbs.len() - 1
    }
}

impl Plaintext {
    pub fn scale(&self) -> f64 {
        self.scale
    }

    pub fn level(&self) -> usize {
        self.data.limbs.len() - 1
    }
}
