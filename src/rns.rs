use std::marker::PhantomData;

use rand::{Rng, RngExt};
use rand_distr::StandardNormal;
use thiserror::Error;

use crate::arith::{ArithError, ModArith};
use crate::ring::{Ring, RingError};

#[derive(Debug)]
pub struct CoeffForm;

#[derive(Debug)]
pub struct NttForm;

#[derive(Debug)]
pub struct Poly<F> {
    pub limbs: Vec<Vec<u64>>,
    _form: PhantomData<F>,
}

impl<F> Poly<F> {
    pub fn new(limbs: Vec<Vec<u64>>) -> Self {
        Poly {
            limbs,
            _form: PhantomData,
        }
    }
}

#[derive(Error, Debug)]
pub enum RnsRingError {
    #[error(transparent)]
    Arith(#[from] ArithError),
    #[error(transparent)]
    Ring(#[from] RingError),
}

#[derive(Debug)]
pub struct RnsRing {
    subrings: Vec<Ring>,
}

impl RnsRing {
    pub fn new(moduli: &[u64], n: usize) -> Result<Self, RnsRingError> {
        let subrings = moduli
            .iter()
            .map(|&q| -> Result<Ring, RnsRingError> {
                let arith = ModArith::new(q)?;
                Ok(Ring::new(arith, n)?)
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(RnsRing { subrings })
    }

    pub fn n(&self) -> usize {
        self.subrings[0].n()
    }

    pub fn num_moduli(&self) -> usize {
        self.subrings.len()
    }

    pub fn modulus(&self, i: usize) -> u64 {
        self.subrings[i].modulus()
    }

    pub fn ntt(&self, p: Poly<CoeffForm>) -> Poly<NttForm> {
        let mut limbs = p.limbs;
        for (limb, ring) in limbs.iter_mut().zip(&self.subrings) {
            ring.forward_ntt_ct(limb);
        }
        Poly::new(limbs)
    }

    pub fn intt(&self, p: Poly<NttForm>) -> Poly<CoeffForm> {
        let mut limbs = p.limbs;
        for (limb, ring) in limbs.iter_mut().zip(&self.subrings) {
            ring.inverse_ntt_gs(limb);
        }
        Poly::new(limbs)
    }

    pub fn add_inplace<F>(&self, a: &mut Poly<F>, b: &Poly<F>) {
        for ((la, lb), ring) in a.limbs.iter_mut().zip(&b.limbs).zip(&self.subrings) {
            ring.arith().add_vec(la, lb);
        }
    }

    pub fn sub_inplace<F>(&self, a: &mut Poly<F>, b: &Poly<F>) {
        for ((la, lb), ring) in a.limbs.iter_mut().zip(&b.limbs).zip(&self.subrings) {
            ring.arith().sub_vec(la, lb);
        }
    }

    pub fn mul_inplace(&self, a: &mut Poly<NttForm>, b: &Poly<NttForm>) {
        for ((la, lb), ring) in a.limbs.iter_mut().zip(&b.limbs).zip(&self.subrings) {
            ring.arith().mul_vec(la, lb);
        }
    }

    pub fn add<F>(&self, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        let mut out = Poly::new(a.limbs.clone());
        self.add_inplace(&mut out, b);
        out
    }

    pub fn sub<F>(&self, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        let mut out = Poly::new(a.limbs.clone());
        self.sub_inplace(&mut out, b);
        out
    }

    pub fn neg<F>(&self, a: &Poly<F>) -> Poly<F> {
        let limbs = a
            .limbs
            .iter()
            .zip(&self.subrings)
            .map(|(la, ring)| la.iter().map(|&x| ring.arith().neg(x)).collect())
            .collect();
        Poly::new(limbs)
    }

    pub fn mul(&self, a: &Poly<NttForm>, b: &Poly<NttForm>) -> Poly<NttForm> {
        let mut out = Poly::new(a.limbs.clone());
        self.mul_inplace(&mut out, b);
        out
    }

    pub fn sample_ternary(&self, rng: &mut impl Rng, h: usize) -> Poly<NttForm> {
        #[derive(Clone, PartialEq)]
        enum Ternary {
            Zero,
            PlusOne,
            MinusOne,
        }

        let n = self.n();
        let mut coefs = vec![Ternary::Zero; n];
        let mut n_set: usize = 0;

        while n_set < h {
            let idx = rng.random_range(0..n);
            if coefs[idx] != Ternary::Zero {
                continue;
            }

            match rng.random_bool(0.5) {
                true => {
                    coefs[idx] = Ternary::PlusOne;
                }
                false => {
                    coefs[idx] = Ternary::MinusOne;
                }
            }
            n_set += 1;
        }

        let q_moduli_count = self.num_moduli();
        let mut sk_coefs = vec![vec![0u64; n]; q_moduli_count];

        for coef_idx in 0..n {
            for (mod_idx, limb) in sk_coefs.iter_mut().enumerate() {
                let modulus = self.modulus(mod_idx);
                limb[coef_idx] = match coefs[coef_idx] {
                    Ternary::PlusOne => 1,
                    Ternary::MinusOne => modulus - 1,
                    Ternary::Zero => 0,
                }
            }
        }

        self.ntt(Poly::<CoeffForm>::new(sk_coefs))
    }

    pub fn sample_uniform(&self, rng: &mut impl Rng) -> Poly<NttForm> {
        let a_limbs: Vec<Vec<u64>> = (0..self.num_moduli())
            .map(|i| {
                let qi = self.modulus(i);
                (0..self.n()).map(|_| rng.random_range(0..qi)).collect()
            })
            .collect();

        Poly::<NttForm>::new(a_limbs)
    }

    // TODO: temp, should be revised later
    pub fn sample_gaussian(&self, rng: &mut impl Rng, sigma: f64) -> Poly<NttForm> {
        let errors: Vec<i64> = (0..self.n())
            .map(|_| (rng.sample::<f64, _>(StandardNormal) * sigma).round() as i64)
            .collect();

        let e_limbs: Vec<Vec<u64>> = (0..self.num_moduli())
            .map(|i| {
                let q = self.modulus(i) as i64;
                errors
                    .iter()
                    .map(|&e| ((e % q) + q) as u64 % q as u64)
                    .collect()
            })
            .collect();

        self.ntt(Poly::<CoeffForm>::new(e_limbs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const MODULI: &[u64] = &[7681, 12289, 40961];
    const N: usize = 256;

    fn make_rns_ring() -> RnsRing {
        RnsRing::new(MODULI, N).unwrap()
    }

    fn rns_poly_vec() -> impl Strategy<Value = Vec<Vec<u64>>> {
        MODULI
            .iter()
            .map(|&q| proptest::collection::vec(0..q, N))
            .collect::<Vec<_>>()
    }

    fn naive_negacyclic(arith: &ModArith, a: &[u64], b: &[u64]) -> Vec<u64> {
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

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(30))]

        #[test]
        fn ntt_roundtrip(limbs in rns_poly_vec()) {
            let rns = make_rns_ring();
            let original = limbs.clone();

            let ntt_form = rns.ntt(Poly::<CoeffForm>::new(limbs));
            let back = rns.intt(ntt_form);

            prop_assert_eq!(back.limbs, original);
        }

        #[test]
        fn add_commutative(a in rns_poly_vec(), b in rns_poly_vec()) {
            let rns = make_rns_ring();
            let a = Poly::<NttForm>::new(a);
            let b = Poly::<NttForm>::new(b);

            let ab = rns.add(&a, &b);
            let ba = rns.add(&b, &a);

            prop_assert_eq!(ab.limbs, ba.limbs);
        }

        #[test]
        fn add_sub_roundtrip(a in rns_poly_vec(), b in rns_poly_vec()) {
            let rns = make_rns_ring();
            let a = Poly::<NttForm>::new(a);
            let b = Poly::<NttForm>::new(b);

            let sum = rns.add(&a, &b);
            let back = rns.sub(&sum, &b);

            prop_assert_eq!(back.limbs, a.limbs);
        }

        #[test]
        fn neg_double_is_identity(a in rns_poly_vec()) {
            let rns = make_rns_ring();
            let a = Poly::<NttForm>::new(a);

            let neg_a = rns.neg(&a);
            let neg_neg_a = rns.neg(&neg_a);

            prop_assert_eq!(neg_neg_a.limbs, a.limbs);
        }

        #[test]
        fn add_neg_is_zero(a in rns_poly_vec()) {
            let rns = make_rns_ring();
            let a = Poly::<NttForm>::new(a);

            let neg_a = rns.neg(&a);
            let sum = rns.add(&a, &neg_a);

            for limb in &sum.limbs {
                prop_assert!(limb.iter().all(|&x| x == 0));
            }
        }

        #[test]
        fn mul_commutative(a in rns_poly_vec(), b in rns_poly_vec()) {
            let rns = make_rns_ring();
            let a = Poly::<NttForm>::new(a);
            let b = Poly::<NttForm>::new(b);

            let ab = rns.mul(&a, &b);
            let ba = rns.mul(&b, &a);

            prop_assert_eq!(ab.limbs, ba.limbs);
        }

        #[test]
        fn ntt_mul_matches_naive(a in rns_poly_vec(), b in rns_poly_vec()) {
            let rns = make_rns_ring();

            let c_ref: Vec<Vec<u64>> = a.iter()
                .zip(&b)
                .zip(MODULI)
                .map(|((la, lb), &q)| naive_negacyclic(&ModArith::new(q).unwrap(), la, lb))
                .collect();

            let a_ntt = rns.ntt(Poly::<CoeffForm>::new(a));
            let b_ntt = rns.ntt(Poly::<CoeffForm>::new(b));
            let c_ntt = rns.mul(&a_ntt, &b_ntt);
            let c = rns.intt(c_ntt);

            prop_assert_eq!(c.limbs, c_ref);
        }

        #[test]
        fn ntt_add_matches_coeff_add(a in rns_poly_vec(), b in rns_poly_vec()) {
            let rns = make_rns_ring();

            let c_ref: Vec<Vec<u64>> = a.iter()
                .zip(&b)
                .zip(MODULI)
                .map(|((la, lb), &q)| {
                    let arith = ModArith::new(q).unwrap();
                    let mut out = la.clone();
                    arith.add_vec(&mut out, lb);
                    out
                })
                .collect();

            let a_ntt = rns.ntt(Poly::<CoeffForm>::new(a));
            let b_ntt = rns.ntt(Poly::<CoeffForm>::new(b));
            let c_ntt = rns.add(&a_ntt, &b_ntt);
            let c = rns.intt(c_ntt);

            prop_assert_eq!(c.limbs, c_ref);
        }
    }
}
