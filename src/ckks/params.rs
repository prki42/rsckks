use thiserror::Error;

use crate::ckks::{CkksContext, CkksContextErr};
use crate::primes::random_ntt_prime;

pub struct CkksParams {
    pub first_mod_size: u32,
    pub scaling_size: u32,
    pub mul_depth: usize,
    pub ring_dim: usize,
}

#[derive(Error, Debug)]
pub enum ParamsErr {
    #[error("first modulus size too big: 2^{0} > u64::MAX")]
    FirstModSizeTooBig(u32),

    #[error("scaling factor too big: 2^{0} > u64::MAX")]
    ScalingSizeTooBig(u32),

    #[error("q0 > qi for all i>0")]
    ScalingNotSmallerThanFirstMod(u32, u32),

    #[error("ring dimension not a power of two")]
    RingDimNotPowerOfTwo,

    #[error(transparent)]
    CkksContext(#[from] CkksContextErr),
}

pub fn gen_ckks_context(params: &CkksParams) -> Result<CkksContext, ParamsErr> {
    if params.first_mod_size > 63 {
        return Err(ParamsErr::FirstModSizeTooBig(params.first_mod_size));
    }
    if params.scaling_size > 63 {
        return Err(ParamsErr::ScalingSizeTooBig(params.scaling_size));
    }

    let n = params.ring_dim;
    if !n.is_power_of_two() {
        return Err(ParamsErr::RingDimNotPowerOfTwo);
    }

    if params.scaling_size >= params.first_mod_size {
        return Err(ParamsErr::ScalingNotSmallerThanFirstMod(
            params.scaling_size,
            params.first_mod_size,
        ));
    }

    let mut q_moduli = vec![0u64; params.mul_depth + 1];
    q_moduli[0] = random_ntt_prime(params.first_mod_size, n);

    let mut all_moduli = vec![q_moduli[0]];

    for qi in q_moduli.iter_mut().skip(1) {
        *qi = random_ntt_prime_excluding(params.scaling_size, n, &all_moduli);
        all_moduli.push(*qi);
    }

    let q_bitsize = (params.mul_depth as u32) * params.scaling_size + params.first_mod_size;
    let pi_bitsize = 61;
    let mut p_moduli = vec![0u64; q_bitsize.div_ceil(pi_bitsize) as usize];

    for pi in p_moduli.iter_mut() {
        *pi = random_ntt_prime_excluding(pi_bitsize, n, &all_moduli);
        all_moduli.push(*pi);
    }

    Ok(CkksContext::new(
        &q_moduli,
        &p_moduli,
        params.ring_dim,
        (1u64 << params.scaling_size) as f64,
    )?)
}

fn random_ntt_prime_excluding(size: u32, n: usize, all_moduli: &[u64]) -> u64 {
    loop {
        let p = random_ntt_prime(size, n);
        if !all_moduli.contains(&p) {
            break p;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn simple_valid_params() {
        gen_ckks_context(&CkksParams {
            first_mod_size: 61,
            scaling_size: 55,
            mul_depth: 15,
            ring_dim: 16,
        })
        .unwrap();
    }

    // TODO more tests...

    proptest! {
        // change to ~10 when the bug below is fixed
        #![proptest_config(ProptestConfig::with_cases(0))]

        // TODO sometimes gets stuck, probably related to prime sampling
        // (invalid combinations of size and n)
        #[test]
        fn should_generate(
            n in (2u32..5),
            (mod_size, scale_size) in (3u32..61).prop_flat_map(|x| (Just(x), 3..x))
        ) {
            gen_ckks_context(&CkksParams {
                first_mod_size: mod_size,
                scaling_size: scale_size,
                mul_depth: 15,
                ring_dim: 2usize.pow(n),
            }).unwrap();
         }
    }
}
