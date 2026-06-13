use rand::RngExt;

/// Miller test (deterministic for `n` < 2^64)
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }

    const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if WITNESSES.contains(&n) {
        return true;
    }

    if n.is_multiple_of(2) {
        return false;
    }

    // n - 1 = d * 2^r
    let mut d = n - 1;
    let mut r = 0u32;
    while d.is_multiple_of(2) {
        d /= 2;
        r += 1;
    }

    'witness: for &a in &WITNESSES {
        // a >= n already handled above (n would be in WITNESSES)
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..(r - 1) {
            x = mod_pow(x, 2, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }

    true
}

/// Generates `size`-bit prime p such that 2*`n` | (p - 1)
pub fn random_ntt_prime(size: u32, n: usize) -> u64 {
    assert!((2..=64).contains(&size));

    let mut rng = rand::rng();
    let mask = if size == 64 {
        u64::MAX
    } else {
        (1u64 << size) - 1
    };
    let msb = 1u64 << (size - 1);

    loop {
        let p = (rng.random::<u64>() & mask) | msb | 1;
        if (p - 1).is_multiple_of(2 * n as u64) && is_prime(p) {
            break p;
        }
    }
}

pub fn distinct_random_ntt_primes(sizes: &[u32], n: usize) -> Vec<u64> {
    let mut res = vec![];
    for &size in sizes {
        loop {
            let p = random_ntt_prime(size, n);
            if !res.contains(&p) {
                res.push(p);
                break;
            }
        }
    }
    res
}

fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let m = modulus as u128;
    let mut result: u128 = 1;
    base %= modulus;
    let mut b = base as u128;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * b % m;
        }
        exp >>= 1;
        b = b * b % m;
    }
    result as u64
}

#[cfg(test)]
mod test {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn small_values() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(5));
    }

    #[test]
    fn known_primes() {
        let primes = [
            7,
            13,
            97,
            541,
            7919,
            104729,
            1299709,
            2147483647,     // Mersenne prime M31
            67280421310721, // prime factor of Fermat's number F6
        ];
        for p in primes {
            assert!(is_prime(p), "{p} should be prime");
        }
    }

    #[test]
    fn pseudoprimes() {
        let pseudoprimes = [
            9,
            15,
            21,
            25,
            33,
            35,
            49,
            85,
            91,
            121,
            133,
            169,
            217,
            221,
            341,
            545,
            781,
            1687,
            2047, // A298756 sequence ("Least strong pseudoprime to base n" for n<=37)
            25326001,
            161304001,
            960946321,
            1157839381,
            3215031751,
            3697278427,
            5764643587,
            6770862367,
            14386156093,
            15579919981,
            18459366157,
            19887974881,
            21276028621,
            27716349961,
            29118033181,
            37131467521,
            41752650241,
            42550716781,
            43536545821, // A056915 sequence ("Strong pseudoprimes to bases 2, 3 and 5")
        ];
        for c in pseudoprimes {
            assert!(!is_prime(c), "{c} should be composite");
        }
    }

    #[test]
    fn large_u64() {
        // Largest prime below 2^64
        assert!(is_prime(u64::MAX - 58));
        assert!(!is_prime(u64::MAX)); // 2^64 - 1 = 3 * 5 * 17 * ...
    }

    #[test]
    fn distinct_prime_generator() {
        let count = 10;
        let sizes = vec![55u32; count];
        let primes = distinct_random_ntt_primes(&sizes, 256);
        assert_eq!(
            primes.len(),
            count,
            "generated {} number of primes instead of {count}",
            primes.len()
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(40))]

        #[test]
        fn random_composite(a in (2..u32::MAX), b in (2..u32::MAX)) {
            assert!(!is_prime((a as u64) * (b as u64)))
        }

        // TODO: test related to ntt primes
    }
}
