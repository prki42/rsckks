use std::marker::PhantomData;

use crate::arith::ModArith;
use crate::ring::Ring;

pub struct CoeffForm;
pub struct NttForm;

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

pub struct RnsRing {
    subrings: Vec<Ring>,
}

impl RnsRing {
    pub fn new(moduli: &[u64], n: usize) -> Self {
        RnsRing {
            subrings: moduli
                .iter()
                .map(|&q| Ring::new(ModArith::new(q), n))
                .collect(),
        }
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

    pub fn add<F>(&self, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        let limbs = a
            .limbs
            .iter()
            .zip(&b.limbs)
            .zip(&self.subrings)
            .map(|((la, lb), ring)| {
                let mut out = la.clone();
                ring.arith().add_vec(&mut out, lb);
                out
            })
            .collect();
        Poly::new(limbs)
    }

    pub fn sub<F>(&self, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        let limbs = a
            .limbs
            .iter()
            .zip(&b.limbs)
            .zip(&self.subrings)
            .map(|((la, lb), ring)| {
                let mut out = la.clone();
                ring.arith().sub_vec(&mut out, lb);
                out
            })
            .collect();
        Poly::new(limbs)
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
        let limbs = a
            .limbs
            .iter()
            .zip(&b.limbs)
            .zip(&self.subrings)
            .map(|((la, lb), ring)| {
                let mut out = la.clone();
                ring.arith().mul_vec(&mut out, lb);
                out
            })
            .collect();
        Poly::new(limbs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const MODULI: &[u64] = &[7681, 12289, 40961];
    const N: usize = 256;

    fn make_rns_ring() -> RnsRing {
        RnsRing::new(MODULI, N)
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
                .map(|((la, lb), &q)| naive_negacyclic(&ModArith::new(q), la, lb))
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
                    let arith = ModArith::new(q);
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
