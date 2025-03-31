#define _POSIX_C_SOURCE 200112

#include "swk.c"
#include "params.c"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if defined(__amd64__)
#include <time.h>
#endif

#define BENCH_WARM (1 <<  2)
#define BENCH_ITER (1 << 11)
#define LEN(a)     (sizeof(a) / sizeof(*(a)))
#define S          "[+]      N      logq    W   W2         t\n"

FILE *f;

static u64
getiter(u64 N)
{
	switch (N) {
	case 1U << 14: return 1U << 11;
	case 1U << 15: return 1U <<  7;
	case 1U << 16: return 1U <<  5;
	case 1U << 17: return 1U <<  3;
	case 1U << 18: return 1U <<  2;
	default: exit(1);
	}
}

static void *
my_calloc(u64 len, u64 size)
{
	void *p = calloc(len, size);
	if (p == 0)
		perror("error: calloc"), exit(1);
	return p;
}

static u64
bench_time(void)
{
	#if defined(__aarch64__)
	u64 cntvct;
	__asm__ volatile("mrs %0, cntvct_el0" : "=r"(cntvct));
	return cntvct;
	#elif defined(__amd64__)
        u64 lo, hi;
	__asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
	return (hi << 32) | lo;
	#else
	#error bench_time
	#endif
}

static u64
bench_freq(void)
{
	#if defined(__aarch64__)
	u64 cntfrq;
	__asm__ volatile("mrs %0, cntfrq_el0" : "=r"(cntfrq));
	return cntfrq;
	#elif defined(__amd64__)
	u64 t0 = bench_time();
	struct timespec ts0 = { 0, 0 }, ts1 = { 0, 0 };
	clock_gettime(CLOCK_MONOTONIC_RAW, &ts0);
	do
		clock_gettime(CLOCK_MONOTONIC_RAW, &ts1);
	while (ts1.tv_sec - ts0.tv_sec < 1);
	u64 t = bench_time() - t0;
	u64 tv_sec = (u64)(ts1.tv_sec - ts0.tv_sec);
	u64 tv_nsec = (u64)(ts1.tv_nsec - ts0.tv_nsec);
	return t * 1000000000ULL / (tv_sec * 1000000000ULL + tv_nsec);
	#else
	#error bench_freq
	#endif
}

static void
prng_init(u64 state[4], u64 seed)
{
	u64 r = (seed += 0x9E3779B97f4A7C15);
	r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9;
	r = (r ^ (r >> 27)) * 0x94D049BB133111EB;
	state[0] = r ^ (r >> 31);

	r = (seed += 0x9E3779B97f4A7C15);
	r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9;
	r = (r ^ (r >> 27)) * 0x94D049BB133111EB;
	state[1] = r ^ (r >> 31);

	r = (seed += 0x9E3779B97f4A7C15);
	r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9;
	r = (r ^ (r >> 27)) * 0x94D049BB133111EB;
	state[2] = r ^ (r >> 31);

	r = (seed += 0x9E3779B97f4A7C15);
	r = (r ^ (r >> 30)) * 0xBF58476D1CE4E5B9;
	r = (r ^ (r >> 27)) * 0x94D049BB133111EB;
	state[3] = r ^ (r >> 31);
}

static u64
prng_squeeze(u64 state[4])
{
	#define ROL64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))
	u64 r = ROL64(state[0] + state[3], 23) + state[0];
	u64 t = state[1] << 17;
	state[2] ^= state[0];
	state[3] ^= state[1];
	state[1] ^= state[2];
	state[0] ^= state[3];
	state[2] ^= t;
	state[3] = ROL64(state[3], 45);
	return r;
}

static int
cmp(void const *op1, void const *op2) {
	u64 t1 = *(u64 *)op1;
	u64 t2 = *(u64 *)op2;
	return (t1 > t2) - (t1 < t2);
}

static void
bench_single(Context ctx, u64 *times, f64 freq)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	u64 len = 2 * W * LK * N, i;
	u64 *swk = my_calloc(len, 8), *p;
	p = swk, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);

	Ciphertext ct, in;
	len = 3 * L * N;
	ct.poly = my_calloc(len, sizeof *ct.poly);
	in.poly = my_calloc(len, sizeof *in.poly);
	in.drop = my_calloc(ctx.L, sizeof *in.drop);
	p = ct.poly, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);
	in.deg = 3, in.ntt = 1;

	i = BENCH_WARM;
	while (i--) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		swk_single(ctx, &in, swk, 2, 1);
	}

	u64 iter = getiter(N);
	for (i = 0; i < iter; ++i) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		u64 t0 = bench_time();
		swk_single(ctx, &in, swk, 2, 1);
		u64 t1 = bench_time();
		times[i] = t1 - t0;
	}

	qsort(times, iter, 8, cmp);
	f64 p50 = (f64)times[iter >> 1] / freq;
	f64 logq = 0.0;
	for (i = 0; i < L; ++i)
		logq += log2((f64)ctx.QPE.mods[i]);
	fprintf(stderr, "[>] %6"PRIu64"   %7.2f   %2"PRIu64"   %2"PRIu64"   %7.4f\n", N, logq, W, (u64)0, p50);
	fwrite(&p50, sizeof p50, 1, f);

	free(swk);
	free(ct.poly);
	free(in.poly);
	free(in.drop);
}

static void
bench_folded(Context ctx, u64 *times, f64 freq)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	u64 len = 2 * W * LK * N, i;
	u64 *fld = my_calloc(len, 8), *p;
	p = fld, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);

	Ciphertext ct, in;
	len = 3 * L * N;
	ct.poly = my_calloc(len, sizeof *ct.poly);
	in.poly = my_calloc(len, sizeof *in.poly);
	in.drop = my_calloc(ctx.L, sizeof *in.drop);
	p = ct.poly, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);
	in.deg = 3, in.ntt = 1;

	i = BENCH_WARM;
	while (i--) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		swk_folded(ctx, &in, fld, 2, 1);
	}

	u64 iter = getiter(N);
	for (i = 0; i < iter; ++i) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		u64 t0 = bench_time();
		swk_folded(ctx, &in, fld, 2, 1);
		u64 t1 = bench_time();
		times[i] = t1 - t0;
	}

	qsort(times, iter, 8, cmp);
	f64 p50 = (f64)times[iter >> 1] / freq;
	f64 logq = 0.0;
	for (i = 0; i < L; ++i)
		logq += log2((f64)ctx.QPE.mods[i]);
	fprintf(stderr, "[>] %6"PRIu64"   %7.2f   %2"PRIu64"   %2"PRIu64"   %7.4f\n", N, logq, W, (u64)0, p50);
	fwrite(&p50, sizeof p50, 1, f);

	free(fld);
	free(ct.poly);
	free(in.poly);
	free(in.drop);
}

static void
bench_double(Context ctx, u64 *times, f64 freq)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 R = ctx.R;
	u64 W = ctx.W;
	u64 W2 = ctx.W2;
	u64 len = 2 * W2 * W * R * N, i;
	u64 *dbl = my_calloc(len, 8), *p;
	p = dbl, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);

	Ciphertext ct, in;
	len = 3 * L * N;
	ct.poly = my_calloc(len, sizeof *ct.poly);
	in.poly = my_calloc(len, sizeof *in.poly);
	in.drop = my_calloc(ctx.L, sizeof *in.drop);
	p = ct.poly, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);
	in.deg = 3, in.ntt = 1;

	i = BENCH_WARM;
	while (i--) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		swk_double(ctx, &in, dbl, 2, 1);
	}

	u64 iter = getiter(N);
	for (i = 0; i < iter; ++i) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		u64 t0 = bench_time();
		swk_double(ctx, &in, dbl, 2, 1);
		u64 t1 = bench_time();
		times[i] = t1 - t0;
	}

	qsort(times, iter, 8, cmp);
	f64 p50 = (f64)times[iter >> 1] / freq;
	f64 logq = 0.0;
	for (i = 0; i < L; ++i)
		logq += log2((f64)ctx.QPE.mods[i]);
	fprintf(stderr, "[>] %6"PRIu64"   %7.2f   %2"PRIu64"   %2"PRIu64"   %7.4f\n", N, logq, W, W2, p50);
	fwrite(&p50, sizeof p50, 1, f);

	free(dbl);
	free(ct.poly);
	free(in.poly);
	free(in.drop);
}

static void
bench_switch(Context ctx, u64 *times, f64 freq, u64 mu, u64 kappa)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	u64 len = 2 * W * LK * N, i;
	u64 *fld = my_calloc(len, 8), *p;
	p = fld, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);

	Ciphertext ct, in;
	len = 3 * L * N;
	ct.poly = my_calloc(len, sizeof *ct.poly);
	in.poly = my_calloc(len, sizeof *in.poly);
	in.drop = my_calloc(ctx.L, sizeof *in.drop);
	p = ct.poly, i = len;
	while (i--)
		*p++ = prng_squeeze(ctx.prng);
	in.deg = 3, in.ntt = 1;

	i = BENCH_WARM;
	while (i--) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		swk_switch(ctx, &in, fld, 2, 1, mu, kappa);
	}

	u64 iter = getiter(N);
	for (i = 0; i < iter; ++i) {
		memcpy(in.poly, ct.poly, len * sizeof *in.poly);
		u64 t0 = bench_time();
		swk_switch(ctx, &in, fld, 2, 1, mu, kappa);
		u64 t1 = bench_time();
		times[i] = t1 - t0;
	}

	qsort(times, iter, 8, cmp);
	f64 p50 = (f64)times[iter >> 1] / freq;
	f64 logq = 0.0;
	for (i = 0; i < L; ++i)
		logq += log2((f64)ctx.QPE.mods[i]);
	fprintf(stderr, "[>] %6"PRIu64"   %7.2f   %2"PRIu64"   %2"PRIu64"   %7.4f\n", N, logq, W, (u64)0, p50);
	fwrite(&p50, sizeof p50, 1, f);

	free(fld);
	free(ct.poly);
	free(in.poly);
	free(in.drop);
}

Context
params2ctx(Params params)
{
	u64 N = params.N;
	u64 L = params.L;
	u64 K = params.K;
	u64 R = params.R;
	u64 W = params.W;
	u64 W2 = params.W2;
	u64 LK = L + K;

	u64 len;
	if (W2 == 0)
		len = (L + 2 * LK + W * LK) * N;
	else
		len = (W2 + 2 * LK + R * (LK + W + 2 * W2)) * N;

	Context ctx = {
		.N = N,
		.L = L,
		.R = R,
		.W = W,
		.W2 = W2,
		.LK = LK,
	};
	ctx.scratch_poly_bytes = len * 8;
	ctx.scratch_cnst = my_calloc(N, 8);
	ctx.scratch_poly = my_calloc(len, 8);
	ctx.QPE = ring_mods_init(params.mods, N, LK + R);
	ctx.p = ring_mods_init(params.p, N, 1);
	prng_init(ctx.prng, 0xcafebabecafebabe);
	return ctx;
}

int
main(int argc, char **argv)
{
	u64 times[BENCH_ITER];
	f64 freq = (f64)bench_freq();

	(void)argc;
	char fname[256];
	sprintf(fname, "%s.dat", argv[0]);
	f = fopen(fname, "w+");
	if (f == 0)
		perror("fopen"), exit(1);

	fputs(S, stderr);
	for (u64 i = 0; i < LEN(params_single); ++i) {
		Context ctx = params2ctx(params_single[i]);
		bench_single(ctx, times, freq);
		ring_mods_free(ctx.QPE);
		ring_mods_free(ctx.p);
		free(ctx.scratch_cnst);
		free(ctx.scratch_poly);
	}

	fputs(S, stderr);
	for (u64 i = 0; i < LEN(params_folded); ++i) {
		Context ctx = params2ctx(params_folded[i]);
		bench_folded(ctx, times, freq);
		ring_mods_free(ctx.QPE);
		ring_mods_free(ctx.p);
		free(ctx.scratch_cnst);
		free(ctx.scratch_poly);
	}

	fputs(S, stderr);
	for (u64 i = 0; i < LEN(params_double); ++i) {
		Context ctx = params2ctx(params_double[i]);
		bench_double(ctx, times, freq);
		ring_mods_free(ctx.QPE);
		ring_mods_free(ctx.p);
		free(ctx.scratch_cnst);
		free(ctx.scratch_poly);
	}

	fputs(S, stderr);
	for (u64 i = 0; i < LEN(params_switch); ++i) {
		u64 mu = params_switch[i].mu, kappa = params_switch[i].kappa;
		Context ctx = params2ctx(params_switch[i]);
		bench_switch(ctx, times, freq, mu, kappa);
		ring_mods_free(ctx.QPE);
		ring_mods_free(ctx.p);
		free(ctx.scratch_cnst);
		free(ctx.scratch_poly);
	}

	return 0;
}
