use thiserror::Error;

#[derive(Debug)]
pub struct ModArith {
    modulus: u64,
}

#[derive(Error, Debug)]
pub enum ArithError {
    #[error("{0} is not a valid modulus, should be < 2^{1}")]
    InvalidModulusSize(u64, u32),
}

impl ModArith {
    pub fn new(modulus: u64) -> Result<Self, ArithError> {
        match modulus >> 61 {
            0 => Ok(ModArith { modulus }),
            _ => Err(ArithError::InvalidModulusSize(modulus, 61)),
        }
    }

    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    #[inline]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.modulus);
        debug_assert!(b < self.modulus);

        let s = a + b;
        if s >= self.modulus {
            s - self.modulus
        } else {
            s
        }
    }

    #[inline]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.modulus);
        debug_assert!(b < self.modulus);

        if a >= b { a - b } else { self.modulus - b + a }
    }

    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.modulus);
        debug_assert!(b < self.modulus);

        (a as u128 * b as u128 % self.modulus as u128) as u64
    }

    pub fn pow(&self, mut base: u64, mut exp: u64) -> u64 {
        let mut result = 1u64;
        base %= self.modulus;
        while exp > 0 {
            if exp & 1 == 1 {
                result = self.mul(result, base);
            }
            exp >>= 1;
            base = self.mul(base, base);
        }
        result
    }

    pub fn inv(&self, a: u64) -> u64 {
        self.pow(a, self.modulus - 2)
    }

    pub fn neg(&self, a: u64) -> u64 {
        if a == 0 { 0 } else { self.modulus - a }
    }

    /// Finds psi such that `order`^psi = 1 mod `self.modulus`
    ///
    /// Assumes `order` | `self.modulus` and `order` power of two
    pub fn primitive_root_of_unity(&self, order: u64) -> Option<u64> {
        debug_assert_eq!(
            (self.modulus - 1) % order,
            0,
            "order {order} does not divide q-1={}",
            self.modulus - 1
        );
        debug_assert!(
            order.is_power_of_two(),
            "order {order} is not a power of two"
        );

        let exp = (self.modulus - 1) / order;
        for g in 2..self.modulus {
            let root = self.pow(g, exp);
            if root != 1 && self.pow(root, order / 2) == self.modulus - 1 {
                return Some(root);
            }
        }
        None
    }

    pub fn add_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        for (lhs, &rhs) in a.iter_mut().zip(b) {
            *lhs = self.add(*lhs, rhs);
        }
    }

    pub fn sub_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        for (lhs, &rhs) in a.iter_mut().zip(b) {
            *lhs = self.sub(*lhs, rhs);
        }
    }

    pub fn mul_vec(&self, a: &mut [u64], b: &[u64]) {
        debug_assert_eq!(a.len(), b.len());
        for (lhs, &rhs) in a.iter_mut().zip(b) {
            *lhs = self.mul(*lhs, rhs);
        }
    }

    pub fn mul_vec_const(&self, a: &mut [u64], c: u64) {
        for el in a.iter_mut() {
            *el = self.mul(*el, c);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const Q: u64 = 7681;

    fn arith() -> ModArith {
        ModArith::new(Q).unwrap()
    }

    fn elem() -> impl Strategy<Value = u64> {
        0..Q
    }

    proptest! {
        #[test]
        fn add_commutative(a in elem(), b in elem()) {
            let m = arith();
            prop_assert_eq!(m.add(a, b), m.add(b, a));
        }

        #[test]
        fn add_identity(a in elem()) {
            let m = arith();
            prop_assert_eq!(m.add(a, 0), a);
        }

        #[test]
        fn add_sub_roundtrip(a in elem(), b in elem()) {
            let m = arith();
            prop_assert_eq!(m.sub(m.add(a, b), b), a);
        }

        #[test]
        fn mul_commutative(a in elem(), b in elem()) {
            let m = arith();
            prop_assert_eq!(m.mul(a, b), m.mul(b, a));
        }

        #[test]
        fn mul_identity(a in elem()) {
            let m = arith();
            prop_assert_eq!(m.mul(a, 1), a);
        }

        #[test]
        fn mul_zero(a in elem()) {
            let m = arith();
            prop_assert_eq!(m.mul(a, 0), 0);
        }

        #[test]
        fn mul_associative(a in elem(), b in elem(), c in elem()) {
            let m = arith();
            prop_assert_eq!(m.mul(m.mul(a, b), c), m.mul(a, m.mul(b, c)));
        }

        #[test]
        fn mul_distributive(a in elem(), b in elem(), c in elem()) {
            let m = arith();
            prop_assert_eq!(
                m.mul(a, m.add(b, c)),
                m.add(m.mul(a, b), m.mul(a, c))
            );
        }

        #[test]
        fn inv_roundtrip(a in 1..Q) {
            let m = arith();
            prop_assert_eq!(m.mul(a, m.inv(a)), 1);
        }

        #[test]
        fn result_in_range(a in elem(), b in elem()) {
            let m = arith();
            prop_assert!(m.add(a, b) < Q);
            prop_assert!(m.sub(a, b) < Q);
            prop_assert!(m.mul(a, b) < Q);
        }
    }
}
