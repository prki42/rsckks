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
pub struct CrossRingPrecomp {
    // (Q/q_i) mod p_k — for Q→P conversion
    pub q_hat_mod_p: Vec<Vec<u64>>,
    // (P/p_k) mod q_i — for P→Q conversion
    pub p_hat_mod_q: Vec<Vec<u64>>,
}

impl CrossRingPrecomp {
    pub fn new(ring_q: &RnsRing, ring_p: &RnsRing) -> Self {
        let q_hat_mod_p = Self::compute_hat_mod_target(ring_q, ring_p);
        let p_hat_mod_q = Self::compute_hat_mod_target(ring_p, ring_q);
        CrossRingPrecomp {
            q_hat_mod_p,
            p_hat_mod_q,
        }
    }

    // For each source modulus i, compute (product of all other source moduli) mod each target modulus
    fn compute_hat_mod_target(source: &RnsRing, target: &RnsRing) -> Vec<Vec<u64>> {
        let src_count = source.num_moduli();
        let tgt_count = target.num_moduli();
        let mut table = vec![vec![0u64; tgt_count]; src_count];

        for (i, table_row) in table.iter_mut().enumerate() {
            for (k, table_elem) in table_row.iter_mut().enumerate() {
                let ak = target.subrings[k].arith();
                let mut product = 1u64;
                for (j, _) in source.subrings.iter().enumerate() {
                    if i != j {
                        product = ak.mul(product, source.modulus(j) % ak.modulus());
                    }
                }
                *table_elem = product;
            }
        }

        table
    }
}

#[derive(Debug)]
pub struct RnsRing {
    subrings: Vec<Ring>,
    hat_inv: Vec<u64>,
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

        let hat_inv = Self::compute_hat_inv(&subrings, moduli);

        Ok(RnsRing { subrings, hat_inv })
    }

    fn compute_hat_inv(subrings: &[Ring], moduli: &[u64]) -> Vec<u64> {
        let k = moduli.len();
        let mut hat_inv = vec![0u64; k];
        for i in 0..k {
            let ai = subrings[i].arith();
            let mut product = 1u64;
            for (j, &qj) in moduli.iter().enumerate() {
                if i != j {
                    product = ai.mul(product, qj % ai.modulus());
                }
            }
            hat_inv[i] = ai.inv(product);
        }
        hat_inv
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

    /* TODO
      Regarding base conversion related code, once i figure out where to keep track of levels
      I should check the code below to make sure it complies with the decision.
    */

    pub fn crt_coefficients(&self, poly: &Poly<CoeffForm>) -> Vec<Vec<u64>> {
        poly.limbs
            .iter()
            .enumerate()
            .map(|(i, limb)| {
                let ai = self.subrings[i].arith();
                limb.iter().map(|&c| ai.mul(c, self.hat_inv[i])).collect()
            })
            .collect()
    }

    /// Approximate base conversion: result may differ from true value
    /// by the k*(product of source moduli) where k < (number of source moduli)
    pub fn base_convert(
        &self,
        crt_coeffs: &[Vec<u64>],
        cross_table: &[Vec<u64>],
    ) -> Poly<CoeffForm> {
        let n = crt_coeffs[0].len();
        let limbs = (0..self.num_moduli())
            .map(|k| {
                let ak = self.subrings[k].arith();
                (0..n)
                    .map(|j| {
                        let mut sum = 0u64;
                        for (i, crt_limb) in crt_coeffs.iter().enumerate() {
                            sum = ak.add(sum, ak.mul(crt_limb[j], cross_table[i][k]));
                        }
                        sum
                    })
                    .collect()
            })
            .collect();
        Poly::new(limbs)
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

            coefs[idx] = match rng.random_bool(0.5) {
                true => Ternary::PlusOne,
                false => Ternary::MinusOne,
            };
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

    const Q_MODULI: &[u64] = &[7681, 12289, 40961];
    const P_MODULI: &[u64] = &[65537, 786433];
    const N: usize = 256;

    fn make_rns_ring() -> RnsRing {
        RnsRing::new(Q_MODULI, N).unwrap()
    }

    fn rns_poly_vec() -> impl Strategy<Value = Vec<Vec<u64>>> {
        Q_MODULI
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
                .zip(Q_MODULI)
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
                .zip(Q_MODULI)
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

    fn make_q_ring() -> RnsRing {
        RnsRing::new(Q_MODULI, N).unwrap()
    }

    fn make_p_ring() -> RnsRing {
        RnsRing::new(P_MODULI, N).unwrap()
    }

    fn coeff_poly(coeff_max_val: u64) -> impl Strategy<Value = Vec<u64>> {
        proptest::collection::vec(0..=coeff_max_val, N)
    }

    fn coeff_poly_to_rns(coeffs: &[u64], ring: &RnsRing) -> Poly<CoeffForm> {
        let limbs = (0..ring.num_moduli())
            .map(|i| {
                let q = ring.modulus(i);
                coeffs.iter().map(|&c| c % q).collect()
            })
            .collect();
        Poly::new(limbs)
    }

    fn is_multiple_of_q_mod_p(diff: u64, q_mod_p: u64, p: u64, max_k: u64) -> bool {
        (0..max_k).any(|k| (k as u128 * q_mod_p as u128 % p as u128) as u64 == diff)
    }

    // Compute (product of all source moduli) mod each target modulus
    fn compute_product_mod_target(source: &RnsRing, target: &RnsRing) -> Vec<u64> {
        (0..target.num_moduli())
            .map(|k| {
                let ak = target.subrings[k].arith();
                let mut product = 1u64;
                for i in 0..source.num_moduli() {
                    product = ak.mul(product, source.modulus(i) % ak.modulus());
                }
                product
            })
            .collect()
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(30))]

        #[test]
        fn base_extend_q_to_p_error_bounded(coeffs in coeff_poly(3000)) {
            let ring_q = make_q_ring();
            let ring_p = make_p_ring();
            let cross = CrossRingPrecomp::new(&ring_q, &ring_p);
            let l = ring_q.num_moduli() as u64;

            let poly_q = coeff_poly_to_rns(&coeffs, &ring_q);
            let crt = ring_q.crt_coefficients(&poly_q);
            let poly_p = ring_p.base_convert(&crt, &cross.q_hat_mod_p);

            let q_mod_p = compute_product_mod_target(&ring_q, &ring_p);

            for k in 0..ring_p.num_moduli() {
                let pk = ring_p.modulus(k);
                for j in 0..N {
                    let diff = (poly_p.limbs[k][j] + pk - coeffs[j]) % pk;
                    prop_assert!(
                        is_multiple_of_q_mod_p(diff, q_mod_p[k], pk, l),
                        "coeff {} mod p_{}: diff {} not k*Q mod p for any k < {}", j, k, diff, l
                    );
                }
            }
        }
    }
}
