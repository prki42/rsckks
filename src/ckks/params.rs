use crate::ckks::CkksContext;

struct CkksParams {
    first_mod_size: u32,
    scaling_size: u32,
    mul_depth: usize,
    ring_dim: usize,
}

enum InvalidParamsErr {}

pub fn gen_ckks_context(params: &CkksParams) -> Result<CkksContext, InvalidParamsErr> {
    Ok(todo!())
}
