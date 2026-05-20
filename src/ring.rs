use crate::arith::ModArith;

pub struct Ring {
    arith: ModArith,
    n: usize,
    zetas: Vec<u64>,
    inv_zetas: Vec<u64>,
    n_inv: u64,
}

impl Ring {
    pub fn new(arith: ModArith, n: usize) -> Self {
        let psi = arith.primitive_root_of_unity(2 * n as u64);
        let log_n = n.trailing_zeros();

        let mut zetas = vec![0u64; n];
        let mut inv_zetas = vec![0u64; n];
        for k in 0..n {
            let exp = Self::bit_rev(k, log_n) as u64;
            zetas[k] = arith.pow(psi, exp);
            inv_zetas[k] = arith.inv(zetas[k]);
        }

        let n_inv = arith.inv(n as u64);

        Ring {
            arith,
            n,
            zetas,
            inv_zetas,
            n_inv,
        }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn modulus(&self) -> u64 {
        self.arith.modulus()
    }

    pub fn arith(&self) -> &ModArith {
        &self.arith
    }

    pub fn forward_ntt_ct(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.n);
        let mut k: usize = 1;
        let mut len = self.n / 2;
        while len >= 1 {
            let mut start = 0;
            while start < self.n {
                let zeta = self.zetas[k];
                k += 1;
                for j in start..(start + len) {
                    let t = self.arith.mul(zeta, a[j + len]);
                    a[j + len] = self.arith.sub(a[j], t);
                    a[j] = self.arith.add(a[j], t);
                }
                start += 2 * len;
            }
            len /= 2;
        }
    }

    pub fn inverse_ntt_gs(&self, a: &mut [u64]) {
        debug_assert_eq!(a.len(), self.n);
        let mut len = 1;
        while len < self.n {
            let k_base = self.n / (2 * len);
            let mut start = 0;
            let mut block = 0;
            while start < self.n {
                let zeta_inv = self.inv_zetas[k_base + block];
                for j in start..(start + len) {
                    let t = a[j];
                    a[j] = self.arith.add(t, a[j + len]);
                    a[j + len] = self.arith.mul(zeta_inv, self.arith.sub(t, a[j + len]));
                }
                start += 2 * len;
                block += 1;
            }
            len *= 2;
        }

        for x in a.iter_mut() {
            *x = self.arith.mul(*x, self.n_inv);
        }
    }

    fn bit_rev(x: usize, log_n: u32) -> usize {
        x.reverse_bits() >> (usize::BITS - log_n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const Q: u64 = 7681;
    const N: usize = 256;

    fn make_ring() -> Ring {
        Ring::new(ModArith::new(Q), N)
    }

    fn poly_vec() -> impl Strategy<Value = Vec<u64>> {
        proptest::collection::vec(0..Q, N)
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
        #![proptest_config(ProptestConfig::with_cases(50))]

        #[test]
        fn ntt_roundtrip(a in poly_vec()) {
            let ring = make_ring();
            let mut poly = a.clone();

            ring.forward_ntt_ct(&mut poly);
            ring.inverse_ntt_gs(&mut poly);

            prop_assert_eq!(poly, a);
        }

        #[test]
        fn ntt_mul_matches_naive(a in poly_vec(), b in poly_vec()) {
            let ring = make_ring();
            let c_ref = naive_negacyclic(ring.arith(), &a, &b);

            let mut a_ntt = a;
            let mut b_ntt = b;
            ring.forward_ntt_ct(&mut a_ntt);
            ring.forward_ntt_ct(&mut b_ntt);
            ring.arith().mul_vec(&mut a_ntt, &b_ntt);
            ring.inverse_ntt_gs(&mut a_ntt);

            prop_assert_eq!(a_ntt, c_ref);
        }

        #[test]
        fn ntt_add_matches_coeff_add(a in poly_vec(), b in poly_vec()) {
            let ring = make_ring();

            let mut c_coeff = a.clone();
            ring.arith().add_vec(&mut c_coeff, &b);

            let mut a_ntt = a;
            let mut b_ntt = b;
            ring.forward_ntt_ct(&mut a_ntt);
            ring.forward_ntt_ct(&mut b_ntt);
            ring.arith().add_vec(&mut a_ntt, &b_ntt);
            ring.inverse_ntt_gs(&mut a_ntt);

            prop_assert_eq!(a_ntt, c_coeff);
        }

        #[test]
        fn ntt_mul_identity(a in poly_vec()) {
            let ring = make_ring();

            let mut one = vec![0u64; N];
            one[0] = 1;
            let mut one_ntt = one;
            ring.forward_ntt_ct(&mut one_ntt);

            let mut a_ntt = a.clone();
            ring.forward_ntt_ct(&mut a_ntt);
            ring.arith().mul_vec(&mut a_ntt, &one_ntt);
            ring.inverse_ntt_gs(&mut a_ntt);

            prop_assert_eq!(a_ntt, a);
        }
    }
}
