use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use num_complex::Complex64;
use rand::RngExt;
use rsckks::ckks::{
    decryptor::Decryptor,
    encoder::Encoder,
    encryptor::Encryptor,
    evaluator::Evaluator,
    keygen::KeyGenerator,
    params::{CkksParams, gen_ckks_context},
};

fn main_bench(c: &mut Criterion) {
    let ctx = gen_ckks_context(&CkksParams {
        first_mod_size: 61,
        scaling_size: 55,
        mul_depth: 12,
        ring_dim: (1 << 15),
    })
    .unwrap();

    let keygen = KeyGenerator::new(&ctx);
    let sk = keygen.secret_key();
    let pk = keygen.public_key(&sk);
    let enc = Encryptor::new(&ctx, &pk);
    let dec = Decryptor::new(&ctx, &sk);
    let encoder = Encoder::new(&ctx);

    let mut rng = rand::rng();

    let slots: usize = 1 << 14;

    let val1: Vec<Complex64> = (0..slots)
        .map(|_| Complex64::new(rng.random(), rng.random()))
        .collect();
    let val2: Vec<Complex64> = (0..slots)
        .map(|_| Complex64::new(rng.random(), rng.random()))
        .collect();

    c.bench_function("encode", |b| {
        b.iter(|| {
            black_box(encoder.encode(&val1, &ctx));
        });
    });

    let pt1 = encoder.encode(&val1, &ctx);
    let pt2 = encoder.encode(&val2, &ctx);

    c.bench_function("encrypt", |b| {
        b.iter(|| {
            black_box(enc.encrypt(&pt1));
        });
    });

    let ct1 = enc.encrypt(&pt1);
    let ct2 = enc.encrypt(&pt2);

    let evaluator = Evaluator::new(&ctx);

    c.bench_function("add", |b| {
        b.iter(|| {
            black_box(evaluator.add(&ct1, &ct2).unwrap());
        });
    });

    c.bench_function("decrypt", |b| {
        b.iter(|| {
            black_box(dec.decrypt(&ct1));
        });
    });

    let pt1_dec = dec.decrypt(&ct1);

    c.bench_function("decode", |b| {
        b.iter(|| {
            black_box(encoder.decode(&pt1_dec, &ctx));
        });
    });
}

criterion_group!(benches, main_bench);
criterion_main!(benches);
