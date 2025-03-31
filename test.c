#include "swk.c"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static void *
my_calloc(u64 len, u64 size)
{
	void *p = calloc(len, size);
	if (p == 0)
		perror("error: calloc"), exit(1);
	return p;
}

/* PRNG
 * https://arxiv.org/abs/1805.01407
 */
#if defined(__amd64__)
#include <x86intrin.h>
#define BSR64(x)    __bsrq(x)
#define POPCNT32(x) __popcntd(x)
#define ROL64(x, r) __rolq(x, r)
#else
#define BSR64(x)    (63 - __builtin_clzll((u64)(x)))
#define POPCNT32(x) __builtin_popcount(x)
#define ROL64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))
#endif

typedef u64 prng_s[4];
struct prng {
	u64 data[4];
};

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

static void
prng_sample_secret(u64 state[4], s8 *rop, u64 len)
{
	s8 lut[] = { -1, 0, 0, 1 };

	len /= 32;
	while (len--) {
		u64 r = prng_squeeze(state);
		u64 i = 32;
		while (i--) {
			*rop++ = lut[r & 3];
			r >>= 2;
		}
	}
}

static void
prng_sample_error(u64 state[4], s8 *rop, u64 len)
{
	while (len--) {
		u64 r = prng_squeeze(state);
		u32 r0 = (r >> 32) & 0x001fffff;
		u32 r1 = r & 0x001fffff;
		*rop++ = (s8)(POPCNT32(r0) - POPCNT32(r1));
	}
}

static void
prng_sample_unimod(u64 state[4], u64 *rop, u64 mod, u64 len)
{
	// NOTE: mod is close to power-of-two,
	//       ignore bias for testing/PoC
	u64 mask = ((u64)2 << BSR64((s64)mod)) - 1;
	while (len--)
		*rop++ = prng_squeeze(state) & mask;
}


/* BGV
 *
 * Q:  ciphertext modulus
 * P:  key switching modulus
 * E:  double-decomposition modulus
 * p:  plaintext modulus
 *
 * N:  polynomial degree
 * L:  number of primes in Q
 * K:  number of primes in P
 * R:  number of primes in E
 * W:  decomposition number
 * W2: double-decomposition number
 */
static void
s8_to_u64(u64 *rop, s8 const *op, u64 mod, u64 len)
{
	while (len--) {
		s64 x = *op++;
		*rop++ = (u64)x + ((u64)(x >> 63) & mod);
	}
}

static void
s8_to_u64_neg(u64 *rop, s8 const *op, u64 mod, u64 len)
{
	while (len--) {
		s64 x = -*op++;
		*rop++ = (u64)x + ((u64)(x >> 63) & mod);
	}
}

static void
normalize(u64 *rop, u64 const *op, u64 mod1, u64 mod2, u64 len)
{
	u64 half = mod1 >> 1;
	u64 norm = mod2 - (mod1 % mod2);
	while (len--) {
		u64 x = *op++;
                *rop++ = (x + (x >= half) * norm) % mod2;
	}
}

void
moddown(Context ctx, Ciphertext *ct, u64 idx)
{
	if (ct->drop[idx] != 0)
		return;
	ct->drop[idx] = 1;

	u64 deg = ct->deg;
	u64 ntt = ct->ntt;

	u64 N = ctx.N;
	u64 L = ctx.L;
	ring_mods qidx = CTX_QPE(idx, 1);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods p = ctx.p;

	u64 (*norm)[N] = ctx.scratch_poly;
	u64 *qinv = ctx.scratch_cnst;
	for (u64 i = 0; i < L; ++i)
		qinv[i] = ring_coef_inv(*qidx.mods, Q.mods[i]);
	u64 tinv = *qidx.mods - ring_coef_inv(*p.mods, *qidx.mods);

	for (u64 d = 0; d < deg; ++d) {
		u64 *poly = ct->poly + d * L * N;
		u64 *pidx = poly + idx * N;
		ring_poly_mulcadd(pidx, pidx, &tinv, 0, qidx, RING_INC_1110);
		if (ntt)
			ring_poly_nttbwd(pidx, qidx);
		for (u64 i = 0; i < L; ++i)
			if (ct->drop[i] == 0)
				normalize(norm[i], pidx, *qidx.mods, Q.mods[i], N);

		for (u64 i = 0; i < L; ++i)
			Q.mods[i] |= RING_MOD_NOP(ct->drop[i]);
		if (ntt)
			ring_poly_nttfwd(*norm, Q);
		ring_poly_mulcadd(poly, *norm, p.mods, poly, Q, RING_INC_1101);
		ring_poly_mulcadd(poly, poly, qinv, 0, Q, RING_INC_1110);
		for (u64 i = 0; i < L; ++i)
			Q.mods[i] &= RING_MOD_USE;
	}
}

static void
encrypt(Context ctx, Ciphertext *ct, s8 const *sk, u64 const *msg)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	ring_mods Q = CTX_QPE(0, L);
	ring_mods p = ctx.p;

	u64 (*poly0)[N] = (u64 (*)[N])ct->poly;
	u64 (*poly1)[N] = poly0 + L;
	u64 (*poly2)[N] = poly1 + L;

	s8 *error = ctx.scratch_poly;
	u64 *m = (u64 *)(error + N);
	u64 (*e)[N] = (u64 (*)[N])(m + N);
	u64 (*s)[N] = e + L;

	memcpy(m, msg, N * sizeof *m);
	ring_poly_nttbwd(m, ctx.p);

	prng_sample_error(ctx.prng, error, N);
	for (u64 i = 0; i < L; ++i)
		s8_to_u64(e[i], error, Q.mods[i], N);
	ring_poly_mulcadd(*e, *e, p.mods, 0, Q, RING_INC_1100);
	ring_poly_add(*e, *e, m, Q, RING_INC_110);
	ring_poly_nttfwd(*e, Q);

	for (u64 i = 0; i < L; ++i)
		s8_to_u64_neg(s[i], sk, Q.mods[i], N);
	ring_poly_nttfwd(*s, Q);
	for (u64 i = 0; i < L; ++i)
		prng_sample_unimod(ctx.prng, poly1[i], Q.mods[i], N);
	ring_poly_muladd(*poly0, *poly1, *s, *e, Q, RING_INC_1111);
	memset(*poly2, 0, L * N * sizeof(u64));
	memset(ct->drop, 0, L * sizeof *ct->drop);
	ct->deg = 2;
	ct->ntt = 1;
}

static void
decrypt(Context ctx, u64 *msg, s8 const *sk, Ciphertext *ct)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	ring_mods Q0 = CTX_QPE(0, 1);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods p = ctx.p;

	// NOTE: switching down to one prime ensures that all
	//       ciphertext primes influence the decryption
	for (u64 i = 1; i < L; ++i)
		moddown(ctx, ct, i);
	u64 *poly0 = ct->poly;
	u64 *poly1 = poly0 + L * N;
	u64 ntt = ct->ntt;

	u64 *s = ctx.scratch_poly;
	s8_to_u64(s, sk, *Q0.mods, N);
	ring_poly_nttfwd(s, Q0);

	if (ntt == 0) {
		ring_poly_nttfwd(poly0, Q0);
		ring_poly_nttfwd(poly1, Q0);
	}
	ring_poly_muladd(msg, poly1, s, poly0, Q0, RING_INC_1111);
	ring_poly_nttbwd(msg, Q0);

	normalize(msg, msg, *Q0.mods, *p.mods, N);
	for (u64 i = 0; i < L; ++i)
		if (ct->drop[i])
			ring_poly_mulcadd(msg, msg, Q.mods + i, 0, p, RING_INC_1110);
	ring_poly_nttfwd(msg, p);
}

static void
genswk(Context ctx, u64 swk[2][ctx.W][ctx.LK][ctx.N], const s8 *sk)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods p = ctx.p;

	s8 (*error)[N] = ctx.scratch_poly;
	u64 (*e)[N] = (u64 (*)[N])(error + W);
	u64 (*s)[N] = e + LK;
	u64 (*s2)[N] = s + LK;

	u64 *P = ctx.scratch_cnst;
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		P[i] = prod(L, K, 1, QP.mods, (u64)-1, LK, mod);
	}

	prng_sample_error(ctx.prng, *error, W * N);
	for (u64 i = 0; i < LK; ++i)
		s8_to_u64_neg(s[i], sk, QP.mods[i], N);
	ring_poly_nttfwd(*s, QP);
	ring_poly_mul(*s2, *s, *s, QP, RING_INC_111);

	for (u64 k = 0; k < W; ++k) {
		for (u64 i = 0; i < LK; ++i)
			s8_to_u64(e[i], error[k], QP.mods[i], ctx.N);
		ring_poly_mulcadd(*e, *e, p.mods, 0, QP, RING_INC_1100);
		ring_poly_nttfwd(*e, QP);
		for (u64 i = 0; i < LK; ++i)
			prng_sample_unimod(ctx.prng, swk[1][k][i], QP.mods[i], N);
		ring_poly_muladd(*swk[0][k], *swk[1][k], *s, *e, QP, RING_INC_1111);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(k != i % W);
		ring_poly_mulcadd(*swk[0][k], *s2, P, *swk[0][k], QP, RING_INC_1111);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}
}

static void
fold_swk(Context ctx, u64 fld[2][ctx.W][ctx.LK][ctx.N], u64 swk[2][ctx.W][ctx.LK][ctx.N])
{
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods p = ctx.p;

	u64 *fold = ctx.scratch_cnst;
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		fold[i] = ring_coef_inv(prod(L, K, 1, QP.mods, i, LK, mod), mod);
		if (i < L)
			continue;
		u64 neginvp = mod - ring_coef_inv(*p.mods, mod);
		fold[i] = ring_coef_mul(fold[i], neginvp, mod);
	}

	for (u64 k = 0; k < W; ++k) {
		ring_poly_mulcadd(*fld[0][k], *swk[0][k], fold, 0, QP, RING_INC_1110);
		ring_poly_mulcadd(*fld[1][k], *swk[1][k], fold, 0, QP, RING_INC_1110);
	}
}

static void
decomp_swk(Context ctx, u64 dbl[2][ctx.W2][ctx.W][ctx.R][ctx.N], u64 swk[2][ctx.W][ctx.LK][ctx.N])
{
	u64 N = ctx.N;
	u64 R = ctx.R;
	u64 W = ctx.W;
	u64 W2 = ctx.W2;
	u64 LK = ctx.LK;
	u64 LKR = LK + ctx.R;
	ring_mods QPE = ctx.QPE;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods E = CTX_QPE(LK, R);

	u64 (*inv)[LKR][N] = ctx.scratch_poly;
	f64 (*approx)[N] = (f64 (*)[N])(inv + W);
	u64 *fbcQP2E_modQP = ctx.scratch_cnst;
	u64 (*fbcQP2E_modE)[R] = (u64 (*)[R])(fbcQP2E_modQP + LK);
	u64 (*rndQP_modE)[R] = fbcQP2E_modE + LK;
	u64 *fold = (u64 *)(rndQP_modE + W2);
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		fbcQP2E_modQP[i] = ring_coef_inv(prod(i % W2, LK, W2, QP.mods, i, LK, mod), mod);
	}
	for (u64 i = 0; i < R; ++i) {
		u64 mod = E.mods[i];
		for (u64 ii = 0; ii < LK; ++ii)
			fbcQP2E_modE[ii][i] = prod(ii % W2, LK, W2, QP.mods, ii, LK, mod);
		for (u64 k2 = 0; k2 < W2; ++k2)
			rndQP_modE[k2][i] = mod - prod(k2, LK, W2, QP.mods, (u64)-1, LK, mod);
		fold[i] = ring_coef_inv(prod(0, R, 1, E.mods, i, R, mod), mod);
	}

	for (u64 k = 0; k < W; ++k) {
		ring_poly_mulcadd(*inv[k], *swk[0][k], fbcQP2E_modQP, 0, QP, RING_INC_1110);
		ring_poly_nttbwd(*inv[k], QP);
	}
	for (u64 k = 0; k < W; ++k) {
		memset(approx, 0, W2 * sizeof *approx);
		for (u64 ii = 0; ii < LK; ++ii) {
			u64 k2 = ii % W2;
			ring_poly_mulcadd(*dbl[0][k2][k], inv[k][ii], fbcQP2E_modE[ii], *dbl[0][k2][k], E, RING_INC_1011);
		}
		for (u64 k2 = 0; k2 < W2; ++k2) {
			ring_poly_divcadd(approx[k2], *inv[k], QPE.mods, approx[k2], N, LK, RING_INC_0110);
			ring_poly_rndmulcadd(*dbl[0][k2][k], approx[k2], rndQP_modE[k2], *dbl[0][k2][k], E, RING_INC_1011);
			ring_poly_nttfwd(*dbl[0][k2][k], E);
		}
	}

	for (u64 k = 0; k < W; ++k) {
		ring_poly_mulcadd(*inv[k], *swk[1][k], fbcQP2E_modQP, 0, QP, RING_INC_1110);
		ring_poly_nttbwd(*inv[k], QP);
	}
	for (u64 k = 0; k < W; ++k) {
		memset(approx, 0, W2 * sizeof *approx);
		for (u64 ii = 0; ii < LK; ++ii) {
			u64 k2 = ii % W2;
			ring_poly_mulcadd(*dbl[1][k2][k], inv[k][ii], fbcQP2E_modE[ii], *dbl[1][k2][k], E, RING_INC_1011);
		}
		for (u64 k2 = 0; k2 < W2; ++k2) {
			ring_poly_divcadd(approx[k2], *inv[k], QPE.mods, approx[k2], N, LK, RING_INC_0110);
			ring_poly_rndmulcadd(*dbl[1][k2][k], approx[k2], rndQP_modE[k2], *dbl[1][k2][k], E, RING_INC_1011);
			ring_poly_nttfwd(*dbl[1][k2][k], E);
		}
	}

	for (u64 k = 0; k < ctx.W; ++k) {
		for (u64 k2 = 0; k2 < ctx.W2; ++k2) {
			ring_poly_mulcadd(*dbl[0][k2][k], *dbl[0][k2][k], fold, 0, E, RING_INC_1110);
			ring_poly_mulcadd(*dbl[1][k2][k], *dbl[1][k2][k], fold, 0, E, RING_INC_1110);
		}
	}
}


/*
 * TEST
 */
#define TEST_N        (1 << 15)
#define TEST_P        65537
#define TEST_MOD_LEN  28
#define TEST_MOD_LENQ 12
#define TEST_MOD_LENP 4
#define TEST_MOD_LENE 12
#define TEST_KSW_W    3
#define TEST_KSW_W2   8

u64 test_p = TEST_P;
u64 test_mod[TEST_MOD_LEN] = {
	288230376154267649,
	288230376155185153,
	288230376155250689,
	288230376156758017,
	288230376157413377,
	288230376158396417,
	288230376160755713,
	288230376161280001,
	288230376161673217,
	288230376161738753,
	288230376162459649,
	288230376164294657,
	288230376166850561,
	288230376167047169,
	288230376168030209,
	288230376168947713,
	1152921504606584833,
	1152921504598720513,
	1152921504597016577,
	1152921504595968001,
	1152921504595640321,
	1152921504593412097,
	1152921504592822273,
	1152921504592429057,
	1152921504589938689,
	1152921504586530817,
	1152921504585547777,
	1152921504583647233
};

static int
test_u64(u64 const *op1, u64 const *op2, u64 len)
{
	while (len--)
		if (*op1++ != *op2++)
			return 0;
	return 1;
}

static void
test(Context ctx)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 R = ctx.R;
	u64 W = ctx.W;
	u64 W2 = ctx.W2;
	u64 LK = L + K;
	ring_mods Q = CTX_QPE(0, L);
	ring_mods p = ctx.p;

	// NOTE: insecure PRNG, but we do not care for testing
	//       as long as its random enough (which it is)
	prng_init(ctx.prng, 0xdeadcafebabebeef);

	u64 *m = my_calloc(N, sizeof *m);
	for (u64 j = 0; j < N; ++j)
		m[j] = prng_squeeze(ctx.prng) % *p.mods;
	u64 *cmp = my_calloc(N, sizeof *cmp);
	for (u64 j = 0; j < N; ++j)
		cmp[j] = (m[j] * m[j]) % *p.mods;

	s8 *sk = my_calloc(N, sizeof *sk);
	prng_sample_secret(ctx.prng, sk, N);

	Ciphertext ct, fwd, bwd;
	ct.drop = my_calloc(L, sizeof *ct.drop);
	fwd.drop = my_calloc(L, sizeof *fwd.drop);
	bwd.drop = my_calloc(L, sizeof *bwd.drop);
	ct.poly = my_calloc(3 * L * N, sizeof *ct.poly);
	fwd.poly = my_calloc(3 * L * N, sizeof *fwd.poly);
	bwd.poly = my_calloc(3 * L * N, sizeof *bwd.poly);

	u64 (*swk)[W][LK][N] = my_calloc(2 * W * LK * N, 8);
	u64 (*fld)[W][LK][N] = my_calloc(2 * W * LK * N, 8);
	u64 (*dbl)[W2][W][R][N] = my_calloc(2 * W2 * W * R * N, 8);
	genswk(ctx, swk, sk);
	fold_swk(ctx, fld, swk);
	decomp_swk(ctx, dbl, swk);

	u64 poly0 = 0;
	u64 poly1 = L * N;
	u64 poly2 = 2 * L * N;
	encrypt(ctx, &fwd, sk, m);
	ring_poly_mul(fwd.poly + poly2, fwd.poly + poly1, fwd.poly + poly1, Q, RING_INC_111);
	ring_poly_mul(fwd.poly + poly1, fwd.poly + poly0, fwd.poly + poly1, Q, RING_INC_111);
	ring_poly_add(fwd.poly + poly1, fwd.poly + poly1, fwd.poly + poly1, Q, RING_INC_111);
	ring_poly_mul(fwd.poly + poly0, fwd.poly + poly0, fwd.poly + poly0, Q, RING_INC_111);
	fwd.deg = 3;

	#define TEST(s, cpy, testexpr) do {                             \
		fputs("[+] " s "... ", stderr);                         \
		memcpy(ct.poly, cpy.poly, 3 * L * N * sizeof *ct.poly); \
		memcpy(ct.drop, cpy.drop, L * sizeof *ct.drop);         \
		ct.deg = cpy.deg;                                       \
		ct.ntt = cpy.ntt;                                       \
		testexpr;                                               \
		decrypt(ctx, m, sk, &ct);                               \
		if (test_u64(m, cmp, N))                                \
			fputs("ok.\n", stderr);                         \
		else {                                                  \
			fputs("not ok.\n", stderr);                     \
			exit(1);                                        \
		}                                                       \
	} while (0)

	memcpy(bwd.poly, fwd.poly, 3 * L * N * sizeof *bwd.poly);
	memcpy(bwd.drop, fwd.drop, L * sizeof *bwd.drop);
	ring_poly_nttbwd(bwd.poly + poly0, Q);
	ring_poly_nttbwd(bwd.poly + poly1, Q);
	ring_poly_nttbwd(bwd.poly + poly2, Q);
	bwd.deg = fwd.deg;
	bwd.ntt = 0;

	TEST("FWD to FWD (single)", fwd, swk_single(ctx, &ct, (u64 *)swk, 2, 1));
	TEST("FWD to BWD (single)", fwd, swk_single(ctx, &ct, (u64 *)swk, 2, 0));
	TEST("BWD to FWD (single)", bwd, swk_single(ctx, &ct, (u64 *)swk, 2, 1));
	TEST("BWD to BWD (single)", bwd, swk_single(ctx, &ct, (u64 *)swk, 2, 0));

	TEST("FWD to FWD (folded)", fwd, swk_folded(ctx, &ct, (u64 *)fld, 2, 1));
	TEST("FWD to BWD (folded)", fwd, swk_folded(ctx, &ct, (u64 *)fld, 2, 0));
	TEST("BWD to FWD (folded)", bwd, swk_folded(ctx, &ct, (u64 *)fld, 2, 1));
	TEST("BWD to BWD (folded)", bwd, swk_folded(ctx, &ct, (u64 *)fld, 2, 0));

	TEST("FWD to FWD (double)", fwd, swk_double(ctx, &ct, (u64 *)dbl, 2, 1));
	TEST("FWD to BWD (double)", fwd, swk_double(ctx, &ct, (u64 *)dbl, 2, 0));
	TEST("BWD to FWD (double)", bwd, swk_double(ctx, &ct, (u64 *)dbl, 2, 1));
	TEST("BWD to BWD (double)", bwd, swk_double(ctx, &ct, (u64 *)dbl, 2, 0));

	moddown(ctx, &fwd, 0);
	moddown(ctx, &fwd, 1);
	memcpy(bwd.poly, fwd.poly, 3 * L * N * sizeof *bwd.poly);
	memcpy(bwd.drop, fwd.drop, L * sizeof *bwd.drop);
	ring_poly_nttbwd(bwd.poly + poly0, Q);
	ring_poly_nttbwd(bwd.poly + poly1, Q);
	ring_poly_nttbwd(bwd.poly + poly2, Q);
	bwd.deg = fwd.deg;
	bwd.ntt = 0;

	u64 mu = 1, kappa = 2;
	TEST("FWD to FWD (sw2to1)", fwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 1, mu, kappa));
	TEST("FWD to BWD (sw2to1)", fwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 0, mu, kappa));
	TEST("BWD to FWD (sw2to1)", bwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 1, mu, kappa));
	TEST("BWD to BWD (sw2to1)", bwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 0, mu, kappa));

	moddown(ctx, &fwd, 2);
	memcpy(bwd.poly, fwd.poly, 3 * L * N * sizeof *bwd.poly);
	memcpy(bwd.drop, fwd.drop, L * sizeof *bwd.drop);
	ring_poly_nttbwd(bwd.poly + poly0, Q);
	ring_poly_nttbwd(bwd.poly + poly1, Q);
	ring_poly_nttbwd(bwd.poly + poly2, Q);
	bwd.deg = fwd.deg;
	bwd.ntt = 0;

	mu = 2, kappa = 3;
	TEST("FWD to FWD (sw3to2)", fwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 1, mu, kappa));
	TEST("FWD to BWD (sw3to2)", fwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 0, mu, kappa));
	TEST("BWD to FWD (sw3to2)", bwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 1, mu, kappa));
	TEST("BWD to BWD (sw3to2)", bwd, swk_switch(ctx, &ct, (u64 *)fld, 2, 0, mu, kappa));
}

int
main(void)
{
	u64 N = TEST_N;
	u64 L = TEST_MOD_LENQ;
	u64 K = TEST_MOD_LENP;
	u64 R = TEST_MOD_LENE;
	u64 W = TEST_KSW_W;
	u64 W2 = TEST_KSW_W2;
	u64 LK = L + K;
	u64 LKR = LK + R;
	Context ctx = {
		.N  = N,
		.L  = L,
		.K  = K,
		.R  = R,
		.W  = W,
		.W2 = W2,
		.LK = LK,
	};
	u64 len = 3 * LKR * LKR * N;
	ctx.scratch_poly_bytes = len * 8;
	ctx.scratch_cnst = my_calloc(N, 8);
	ctx.scratch_poly = my_calloc(len, 8);
	ctx.QPE = ring_mods_init(test_mod, N, LKR);
	ctx.p = ring_mods_init(&test_p, N, 1);

	test(ctx);

	return 0;
}
