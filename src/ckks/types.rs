use crate::rns::{NttForm, Poly};

#[derive(Debug)]
pub struct SecretKey {
    pub(crate) sk: Poly<NttForm>,
}

pub struct PublicKey {
    pub(crate) a: Poly<NttForm>,
    pub(crate) b: Poly<NttForm>,
}

pub struct RelinKey {
    pub(crate) keys: Vec<(Poly<NttForm>, Poly<NttForm>)>,
}

pub struct Plaintext {
    pub(crate) data: Poly<NttForm>,
    pub(crate) scale: f64,
}

pub struct Ciphertext {
    pub(crate) c0: Poly<NttForm>,
    pub(crate) c1: Poly<NttForm>,
    pub(crate) scale: f64,
    pub(crate) level: usize,
}

impl Ciphertext {
    pub fn scale(&self) -> f64 {
        self.scale
    }

    pub fn level(&self) -> usize {
        self.level
    }
}

impl Plaintext {
    pub fn scale(&self) -> f64 {
        self.scale
    }
}
