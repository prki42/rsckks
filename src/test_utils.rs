use std::sync::{Arc, LazyLock};

use num_complex::Complex64;
use proptest::prelude::*;
use proptest::test_runner::{FileFailurePersistence, TestCaseError};

use crate::arith::ModArith;
use crate::ckks::CkksContext;
use crate::ckks::decryptor::Decryptor;
use crate::ckks::encoder::Encoder;
use crate::ckks::encryptor::Encryptor;
use crate::ckks::keygen::KeyGenerator;
use crate::ckks::params::{CkksParams, gen_ckks_context};
use crate::ckks::types::{Ciphertext, SecretKey};

pub fn no_persist(cases: u32) -> ProptestConfig {
    ProptestConfig {
        cases,
        failure_persistence: Some(Box::new(FileFailurePersistence::Off)),
        ..ProptestConfig::default()
    }
}

pub struct CkksTestEnv {
    pub ctx: Arc<CkksContext>,
    pub encoder: Encoder,
    pub sk: SecretKey,
}

impl CkksTestEnv {
    pub fn slots(&self) -> usize {
        self.ctx.n() / 2
    }

    pub fn encrypt(&self, values: &[Complex64]) -> Ciphertext {
        let keygen = KeyGenerator::new(&self.ctx);
        let pk = keygen.public_key(&self.sk);
        let enc = Encryptor::new(&self.ctx, &pk);
        let pt = self.encoder.encode(values, &self.ctx);
        enc.encrypt(&pt)
    }

    pub fn decrypt_decode(&self, ct: &Ciphertext) -> Vec<Complex64> {
        let dec = Decryptor::new(&self.ctx, &self.sk);
        let pt = dec.decrypt(ct);
        self.encoder.decode(&pt, &self.ctx)
    }
}

impl From<Arc<CkksContext>> for CkksTestEnv {
    fn from(ctx: Arc<CkksContext>) -> Self {
        let encoder = Encoder::new(&ctx);
        let keygen = KeyGenerator::new(&ctx);
        let sk = keygen.secret_key();
        Self { ctx, encoder, sk }
    }
}

// NOTE: until I expose different ctx for encoder tests i must keep 2*scaling_size < first_mod_size
pub static TEST_PARAMS: &[CkksParams] = &[
    CkksParams {
        first_mod_size: 61,
        scaling_size: 27,
        mul_depth: 5,
        ring_dim: 16,
    },
    CkksParams {
        first_mod_size: 20,
        scaling_size: 9,
        mul_depth: 1,
        ring_dim: 16,
    },
];

pub static TEST_CONTEXTS: LazyLock<Vec<Arc<CkksContext>>> = LazyLock::new(|| {
    TEST_PARAMS
        .iter()
        .map(|params| gen_ckks_context(params).map(Arc::new))
        .collect::<Result<Vec<_>, _>>()
        .expect("TODO")
});

pub static CKKS_TEST_ENVS: LazyLock<Vec<CkksTestEnv>> = LazyLock::new(|| {
    TEST_CONTEXTS
        .iter()
        .map(|ctx| CkksTestEnv::from(Arc::clone(ctx)))
        .collect()
});

#[derive(Clone, Copy)]
pub struct TestEnvRef(usize);

impl std::fmt::Debug for TestEnvRef {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let p = &TEST_PARAMS[self.0];
        // TODO
        write!(f, "CkksCtx(N={}, depth={})", p.ring_dim, p.mul_depth)
    }
}

impl std::ops::Deref for TestEnvRef {
    type Target = CkksTestEnv;
    fn deref(&self) -> &CkksTestEnv {
        &CKKS_TEST_ENVS[self.0]
    }
}

pub fn with_test_env<S, T>(f: impl Fn(&CkksTestEnv) -> S) -> impl Strategy<Value = (TestEnvRef, T)>
where
    S: Strategy<Value = T>,
    T: std::fmt::Debug,
{
    (0..CKKS_TEST_ENVS.len()).prop_flat_map(move |i| (Just(TestEnvRef(i)), f(&CKKS_TEST_ENVS[i])))
}

pub fn complex_vec(max_val: f64, len: usize) -> impl Strategy<Value = Vec<Complex64>> {
    proptest::collection::vec((-max_val..max_val, -max_val..max_val), len..=len).prop_map(|v| {
        v.into_iter()
            .map(|(re, im)| Complex64::new(re, im))
            .collect()
    })
}

pub fn real_vec(max_val: f64, len: usize) -> impl Strategy<Value = Vec<Complex64>> {
    proptest::collection::vec(-max_val..max_val, len..=len)
        .prop_map(|v| v.into_iter().map(|re| Complex64::new(re, 0.0)).collect())
}

pub fn assert_slots_approx(
    expected: &[Complex64],
    actual: &[Complex64],
    bound: f64,
) -> Result<(), TestCaseError> {
    let len = expected.len().min(actual.len());
    for i in 0..len {
        let diff = (expected[i] - actual[i]).norm();
        if diff >= bound {
            return Err(TestCaseError::fail(format!(
                "slot {i}: expected {}, got {}, diff {diff} >= bound {bound}",
                expected[i], actual[i],
            )));
        }
    }
    Ok(())
}

pub fn naive_negacyclic(arith: &ModArith, a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len();
    let mut c = vec![0u64; n];
    for i in 0..n {
        for j in 0..n {
            if i + j < n {
                c[i + j] = arith.add(c[i + j], arith.mul(a[i], b[j]));
            } else {
                c[i + j - n] = arith.sub(c[i + j - n], arith.mul(a[i], b[j]));
            }
        }
    }
    c
}
