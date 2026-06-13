use std::f64::consts::PI;

use num_complex::Complex64;

use crate::ckks::CkksContext;
use crate::ckks::types::Plaintext;
use crate::rns::{CoeffForm, NttForm, Poly};

// TODO: ponder upon my decision to not borrow context??

pub struct Encoder {
    n: usize,
    slots: usize,
    scale: f64,
    zetas: Vec<Complex64>,
    bit_rev: Vec<usize>,
}

impl Encoder {
    pub fn new(ctx: &CkksContext) -> Self {
        let n = ctx.n();
        let slots = n / 2;
        let m = 2 * n;

        let zetas: Vec<Complex64> = (0..m)
            .map(|k| {
                let angle = 2.0 * PI * k as f64 / m as f64;
                Complex64::new(angle.cos(), angle.sin())
            })
            .collect();

        let log_n = n.trailing_zeros();
        let bit_rev: Vec<usize> = (0..n)
            .map(|i| i.reverse_bits() >> (usize::BITS - log_n))
            .collect();

        Encoder {
            n,
            slots,
            scale: ctx.scale,
            zetas,
            bit_rev,
        }
    }

    pub fn slots(&self) -> usize {
        self.slots
    }

    pub fn encode(&self, values: &[Complex64], ctx: &CkksContext) -> Plaintext {
        assert!(values.len() <= self.slots);

        let mut v = vec![Complex64::new(0.0, 0.0); self.n];
        for (i, &z) in values.iter().enumerate() {
            v[i] = z;
            v[self.n - 1 - i] = z.conj();
        }

        self.ifft(&mut v);

        let ring_q = ctx.ring_q();
        let limbs: Vec<Vec<u64>> = (0..ring_q.num_moduli())
            .map(|l| {
                let q = ring_q.modulus(l);
                (0..self.n)
                    .map(|j| {
                        let coeff = (v[j] * self.zetas[j].conj()).re * self.scale;
                        Self::f64_to_mod(coeff, q)
                    })
                    .collect()
            })
            .collect();

        let poly_ntt = ring_q.ntt(Poly::<CoeffForm>::new(limbs));
        Plaintext {
            data: poly_ntt,
            scale: self.scale,
        }
    }

    pub fn decode(&self, pt: &Plaintext, ctx: &CkksContext) -> Vec<Complex64> {
        let ring_q = ctx.ring_q();
        let poly_coeff = ring_q.intt(Poly::<NttForm>::new(pt.data.limbs.clone()));

        let q = ring_q.modulus(0);
        let half_q = q / 2;

        let mut v: Vec<Complex64> = poly_coeff.limbs[0]
            .iter()
            .enumerate()
            .map(|(j, &c)| {
                let centered = if c > half_q {
                    -((q - c) as f64)
                } else {
                    c as f64
                };
                Complex64::new(centered / pt.scale, 0.0) * self.zetas[j]
            })
            .collect();

        self.fft(&mut v);

        v[..self.slots].to_vec()
    }

    fn fft(&self, a: &mut [Complex64]) {
        debug_assert_eq!(a.len(), 2 * self.slots, "invalid number of slots");
        self.bit_reverse_permute(a);

        let mut len = 2;
        while len <= self.n {
            let half = len / 2;
            let step = 2 * self.n / len;
            for start in (0..self.n).step_by(len) {
                for j in 0..half {
                    let w = self.zetas[j * step];
                    let u = a[start + j];
                    let t = a[start + j + half] * w;
                    a[start + j] = u + t;
                    a[start + j + half] = u - t;
                }
            }
            len *= 2;
        }
    }

    fn ifft(&self, a: &mut [Complex64]) {
        debug_assert_eq!(a.len(), 2 * self.slots, "invalid number of slots");
        self.bit_reverse_permute(a);

        let mut len = 2;
        while len <= self.n {
            let half = len / 2;
            let step = 2 * self.n / len;
            for start in (0..self.n).step_by(len) {
                for j in 0..half {
                    let w = self.zetas[j * step].conj();
                    let u = a[start + j];
                    let t = a[start + j + half] * w;
                    a[start + j] = u + t;
                    a[start + j + half] = u - t;
                }
            }
            len *= 2;
        }

        let inv_n = 1.0 / self.n as f64;
        for x in a.iter_mut() {
            *x *= inv_n;
        }
    }

    fn bit_reverse_permute(&self, a: &mut [Complex64]) {
        for i in 0..self.n {
            let j = self.bit_rev[i];
            if i < j {
                a.swap(i, j);
            }
        }
    }

    fn f64_to_mod(x: f64, q: u64) -> u64 {
        let r = x.round() as i128;
        let q = q as i128;
        ((r % q + q) % q) as u64
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::test_utils::*;
    use proptest::prelude::*;

    impl Encoder {
        pub fn encoding_err_upper_bound(&self) -> f64 {
            (self.n as f64) / (2.0 * self.scale)
        }

        pub fn mul_err_upper_bound(&self, inf_norm1: f64, inf_norm2: f64) -> f64 {
            (((self.n as f64) * (inf_norm1 + inf_norm2)) / (self.scale * 2.0))
                + ((self.n as f64).powf(2.0) / (4.0 * self.scale.powf(2.0)))
        }
    }

    proptest! {
        #![proptest_config(no_persist(30))]

        // TODO: use custom contexts for encoder (mul test has specific requirements)

        #[test]
        fn fft_ifft_roundtrip(
            (tc, (re, im)) in with_test_env(|tc| {
                let s = tc.slots();
                (proptest::collection::vec(-100.0..100.0, 2*s),
                    proptest::collection::vec(-100.0..100.0, 2*s))
            })
        ) {
            let original: Vec<Complex64> = re.into_iter()
                .zip(im)
                .map(|(r, i)| Complex64::new(r, i))
                .collect();
            let mut data = original.clone();

            tc.encoder.fft(&mut data);
            tc.encoder.ifft(&mut data);

            for (k, (&expected, &actual)) in original.iter().zip(data.iter()).enumerate() {
                let diff = (expected - actual).norm();
                prop_assert!(diff < 1e-6, "index {k}: expected {expected}, got {actual}");
            }
        }

        #[test]
        fn encode_decode_complex(
            (tc, values) in with_test_env(|tc| complex_vec(1.0, tc.slots()))
        ) {
            let pt = tc.encoder.encode(&values, &tc.ctx);
            let decoded = tc.encoder.decode(&pt, &tc.ctx);
            let bound = tc.encoder.encoding_err_upper_bound();
            assert_slots_approx(&values, &decoded, bound)?;
        }

        #[test]
        fn encode_decode_real(
            (tc, values) in with_test_env(|tc| real_vec(1.0, tc.slots()))
        ) {
            let pt = tc.encoder.encode(&values, &tc.ctx);
            let decoded = tc.encoder.decode(&pt, &tc.ctx);
            let bound = tc.encoder.encoding_err_upper_bound();
            assert_slots_approx(&values, &decoded, bound)?;
        }

        #[test]
        fn encode_add_decode(
            (tc, (a, b)) in with_test_env(|tc| {
                let s = tc.slots();
                (complex_vec(1.0, s), complex_vec(1.0, s))
            })
        ) {
            let ring_q = tc.ctx.ring_q();

            let pt_a = tc.encoder.encode(&a, &tc.ctx);
            let pt_b = tc.encoder.encode(&b, &tc.ctx);

            let sum_poly = ring_q.add(&pt_a.data, &pt_b.data);
            let pt_sum = Plaintext { data: sum_poly, scale: pt_a.scale };
            let decoded = tc.encoder.decode(&pt_sum, &tc.ctx);

            let expected: Vec<_> = a.iter().zip(&b).map(|(x, y)| x + y).collect();
            let bound = tc.encoder.encoding_err_upper_bound() * 2.0;
            assert_slots_approx(&expected, &decoded, bound)?;
        }

        #[test]
        fn encode_mul_decode(
            (tc, (a, b)) in with_test_env(|tc| {
                let s = tc.slots();
                (complex_vec(1.0, s), complex_vec(1.0, s))
            })
        ) {
            let ring_q = tc.ctx.ring_q();

            let pt_a = tc.encoder.encode(&a, &tc.ctx);
            let pt_b = tc.encoder.encode(&b, &tc.ctx);

            let prod_poly = ring_q.mul(&pt_a.data, &pt_b.data);
            let pt_prod = Plaintext { data: prod_poly, scale: pt_a.scale * pt_a.scale };
            let decoded = tc.encoder.decode(&pt_prod, &tc.ctx);

            let expected: Vec<_> = a.iter().zip(&b).map(|(x, y)| x * y).collect();
            let inf_norm = |z: &[Complex64]| -> f64 {
                z.iter().map(|z| z.norm()).fold(0.0_f64, f64::max)
            };
            let bound = tc.encoder.mul_err_upper_bound(inf_norm(&a), inf_norm(&b));
            assert_slots_approx(&expected, &decoded, bound)?;
        }
    }
}
