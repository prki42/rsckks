use std::f64::consts::PI;

use num_complex::Complex64;

use crate::ckks::CkksContext;
use crate::ckks::types::Plaintext;
use crate::rns::{CoeffForm, NttForm, Poly};

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
        debug_assert_eq!(a.len(), self.n);
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
        debug_assert_eq!(a.len(), self.n);
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
mod tests {
    use super::*;
    use proptest::prelude::*;

    const N: usize = 256;
    const SLOTS: usize = N / 2;

    fn make_ctx() -> CkksContext {
        CkksContext::new(&[998244353, 985661441, 754974721], &[469762049], N, 64.0).unwrap()
    }

    fn complex_vec(max_val: f64, max_len: usize) -> impl Strategy<Value = Vec<Complex64>> {
        proptest::collection::vec((-max_val..max_val, -max_val..max_val), 1..=max_len).prop_map(
            |v| {
                v.into_iter()
                    .map(|(re, im)| Complex64::new(re, im))
                    .collect()
            },
        )
    }

    fn real_vec(max_val: f64, max_len: usize) -> impl Strategy<Value = Vec<Complex64>> {
        proptest::collection::vec(-max_val..max_val, 1..=max_len)
            .prop_map(|v| v.into_iter().map(|re| Complex64::new(re, 0.0)).collect())
    }

    impl Encoder {
        fn encoding_error(&self) -> f64 {
            (self.n as f64) / (2.0 * self.scale)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(30))]

        #[test]
        fn fft_ifft_roundtrip(
            re in proptest::collection::vec(-100.0f64..100.0, N),
            im in proptest::collection::vec(-100.0f64..100.0, N),
        ) {
            let ctx = make_ctx();
            let encoder = Encoder::new(&ctx);

            let original: Vec<Complex64> = re.into_iter()
                .zip(im)
                .map(|(r, i)| Complex64::new(r, i))
                .collect();
            let mut data = original.clone();

            encoder.fft(&mut data);
            encoder.ifft(&mut data);

            for (k, (&expected, &actual)) in original.iter().zip(data.iter()).enumerate() {
                let diff = (expected - actual).norm();
                prop_assert!(diff < 1e-6, "index {k}: expected {expected}, got {actual}");
            }
        }

        // TODO: better tests for encoding, tighter/proper bounds

        #[test]
        fn encode_decode_complex_values(values in complex_vec(100.0, SLOTS)) {
            let ctx = make_ctx();
            let encoder = Encoder::new(&ctx);

            let pt = encoder.encode(&values, &ctx);
            let decoded = encoder.decode(&pt, &ctx);

            for (i, (&expected, &actual)) in values.iter().zip(decoded.iter()).enumerate() {
                let diff = (expected - actual).norm();
                prop_assert!(
                    diff < encoder.encoding_error(),
                    "slot {i}: expected {expected}, got {actual}, diff {diff}"
                );
            }
        }

        #[test]
        fn encode_decode_real_values(values in real_vec(100.0, SLOTS)) {
            let ctx = make_ctx();
            let encoder = Encoder::new(&ctx);

            let pt = encoder.encode(&values,  &ctx);
            let decoded = encoder.decode(&pt, &ctx);

            for (i, (&expected, &actual)) in values.iter().zip(decoded.iter()).enumerate() {
                let diff = (expected - actual).norm();
                prop_assert!(
                    diff < encoder.encoding_error(),
                    "slot {i}: expected {expected}, got {actual}, diff {diff}"
                );
            }
        }

        #[test]
        fn encode_add_decode(
            a in complex_vec(100.0, SLOTS),
            b in complex_vec(100.0, SLOTS),
        ) {
            let ctx = make_ctx();
            let encoder = Encoder::new(&ctx);
            let ring_q = ctx.ring_q();

            let len = a.len().min(b.len());
            let pt_a = encoder.encode(&a,  &ctx);
            let pt_b = encoder.encode(&b,  &ctx);

            let sum_poly = ring_q.add(&pt_a.data, &pt_b.data);
            let pt_sum = Plaintext { data: sum_poly, scale: pt_a.scale };
            let decoded = encoder.decode(&pt_sum, &ctx);

            for i in 0..len {
                let expected = a[i] + b[i];
                let diff = (expected - decoded[i]).norm();
                prop_assert!(
                    diff < encoder.encoding_error() * 2.0,
                    "slot {i}: expected {expected}, got {}, diff {diff}", decoded[i]
                );
            }
        }

        #[test]
        fn encode_mul_decode(
            a in complex_vec(3.0, SLOTS),
            b in complex_vec(3.0, SLOTS),
        ) {
            let ctx = make_ctx();
            let encoder = Encoder::new(&ctx);
            let ring_q = ctx.ring_q();

            let len = a.len().min(b.len());
            let pt_a = encoder.encode(&a, &ctx);
            let pt_b = encoder.encode(&b, &ctx);

            let prod_poly = ring_q.mul(&pt_a.data, &pt_b.data);
            let pt_prod = Plaintext { data: prod_poly, scale: pt_a.scale * pt_a.scale } ;
            let decoded = encoder.decode(&pt_prod, &ctx);

            for i in 0..len {
                let expected = a[i] * b[i];
                let diff = (expected - decoded[i]).norm();
                prop_assert!(
                    // TODO: temp
                    diff < 3.0,
                    "slot {i}: expected {expected}, got {}, diff {diff}", decoded[i]
                );
            }
        }
    }
}
