#ifdef __cplusplus
extern "C" {
	#include "ring.h"
}
#else
#include "ring.h"
#endif
#include <string.h>

#define CTX_QPE(idx, len) (ring_mods){ ctx.QPE.precomp + idx, ctx.QPE.mods + idx, ctx.N, len }

typedef struct context {
	u64 prng[4];
	void *scratch_poly;
	void *scratch_cnst;
	ring_mods QPE, p;
	u64 scratch_poly_bytes;
	u64 N, L, K, R, W, W2, LK;
} Context;

typedef struct ciphertext {
	u64 *poly;
	u64 *drop;
	u64 deg;
	u64 ntt;
} Ciphertext;

static u64
prod(u64 i0, u64 len, u64 inc, u64 *op, u64 skip, u64 max, u64 mod)
{
	len += i0;
	u64 rop = 1;
	for (u64 i = i0; i < len && i < max; i += inc)
		if (i != skip)
			rop = ring_coef_mul(rop, op[i], mod);
	return rop;
}

void
swk_single(Context ctx, Ciphertext *ct, u64 *key, u64 dsw, u64 domain)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods P = CTX_QPE(L, K);
	ring_mods p = ctx.p;

	u64 (*swk)[W][LK][N] = (u64 (*)[W][LK][N])key;
	u64 (*poly0)[N] = (u64 (*)[N])ct->poly;
	u64 (*poly1)[N] = poly0 + L;
	u64 (*poly2)[N] = poly1 + L;
	u64 (*polysw)[N] = (u64 (*)[N])(poly0 + dsw * L);
	u64 ntt = ct->ntt;

	memset(ctx.scratch_poly, 0, ctx.scratch_poly_bytes);
	u64 (*tmp)[N] = (u64 (*)[N])ctx.scratch_poly;
	u64 (*dot0)[N] = tmp + L;
	u64 (*dot1)[N] = dot0 + LK;
	u64 (*ext)[LK][N] = (u64 (*)[LK][N])(dot1 + LK);

	u64 *P_modQ = (u64 *)ctx.scratch_cnst;
	u64 *invP_modQ = P_modQ + L;
	u64 *neginvp_modP = invP_modQ + L;
	u64 *fbcQ2QP_modQ = neginvp_modP + K;
	u64 (*fbcQ2QP_modQP)[LK] = (u64 (*)[LK])(fbcQ2QP_modQ + L);
	u64 *fbcP2Q_modP = (u64 *)(fbcQ2QP_modQP + L);
	u64 (*fbcP2Q_modQ)[L] = (u64 (*)[L])(fbcP2Q_modP + K);
	for (u64 i = 0; i < L; ++i) {
		u64 mod = Q.mods[i];
		P_modQ[i] = prod(0, K, 1, P.mods, (u64)-1, K, mod);
		invP_modQ[i] = ring_coef_inv(P_modQ[i], mod);
		fbcQ2QP_modQ[i] = ring_coef_inv(prod(i % W, L, W, Q.mods, i, L, mod), mod);
		for (u64 ii = 0; ii < L; ++ii)
			fbcQ2QP_modQP[ii][i] = prod(ii % W, L, W, Q.mods, ii, L, mod);
		for (u64 ii = 0; ii < K; ++ii)
			fbcP2Q_modQ[ii][i] = prod(0, K, 1, P.mods, ii, K, mod);
	}
	for (u64 i = 0; i < K; ++i) {
		u64 mod = P.mods[i];
		neginvp_modP[i] = mod - ring_coef_inv(*p.mods, mod);
		fbcP2Q_modP[i] = ring_coef_inv(prod(0, K, 1, P.mods, i, K, mod), mod);
	}
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		for (u64 ii = 0; ii < L; ++ii)
			fbcQ2QP_modQP[ii][i] = prod(ii % W, L, W, Q.mods, ii, L, mod);
	}

	for (u64 i = 0; i < L; ++i)
		memcpy(ext[i % W][i], polysw[i], N * sizeof(u64));
	ring_poly_mulcadd(*polysw, *polysw, fbcQ2QP_modQ, 0, Q, RING_INC_1110);
	if (ntt == 1)
		ring_poly_nttbwd(*polysw, Q);
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W);
		for (u64 ii = k; ii < L; ii += W)
			ring_poly_mulcadd(*ext[k], polysw[ii], fbcQ2QP_modQP[ii], *ext[k], QP, RING_INC_1011);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W && ntt == 1);
		ring_poly_nttfwd(*ext[k], QP);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}

	for (u64 k = 0; k < W; ++k) {
		ring_poly_muladd(*dot0, *ext[k], *swk[0][k], *dot0, QP, RING_INC_1111);
		ring_poly_muladd(*dot1, *ext[k], *swk[1][k], *dot1, QP, RING_INC_1111);
	}
	if (ntt == 1)
		ring_poly_mulcadd(*dot0, *poly0, P_modQ, *dot0, Q, RING_INC_1111);
	if (dsw == 2 && ntt == 1)
		ring_poly_mulcadd(*dot1, *poly1, P_modQ, *dot1, Q, RING_INC_1111);

	ring_poly_nttbwd(dot0[L], P);
	ring_poly_nttbwd(dot1[L], P);
	ring_poly_mulcadd(dot0[L], dot0[L], neginvp_modP, 0, P, RING_INC_1110);
	ring_poly_mulcadd(dot1[L], dot1[L], neginvp_modP, 0, P, RING_INC_1110);

	ring_poly_mulcadd(dot0[L], dot0[L], fbcP2Q_modP, 0, P, RING_INC_1110);
	ring_poly_mulcadd(dot1[L], dot1[L], fbcP2Q_modP, 0, P, RING_INC_1110);

	memset(*tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*tmp, dot0[L + ii], fbcP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	ring_poly_mulcadd(*tmp, *tmp, p.mods, 0, Q, RING_INC_1100);
	if (ntt == 0)
		ring_poly_mulcadd(*tmp, *poly0, P_modQ, *tmp, Q, RING_INC_1111);
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	else
		ring_poly_nttbwd(*dot0, Q);
	ring_poly_add(*poly0, *dot0, *tmp, Q, RING_INC_111);
	ring_poly_mulcadd(*poly0, *poly0, invP_modQ, 0, Q, RING_INC_1110);
	memset(*tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*tmp, dot1[L + ii], fbcP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	ring_poly_mulcadd(*tmp, *tmp, p.mods, 0, Q, RING_INC_1100);
	if (dsw == 2 && ntt == 0)
		ring_poly_mulcadd(*tmp, *poly1, P_modQ, *tmp, Q, RING_INC_1111);
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	else
		ring_poly_nttbwd(*dot1, Q);
	ring_poly_add(*poly1, *dot1, *tmp, Q, RING_INC_111);
	ring_poly_mulcadd(*poly1, *poly1, invP_modQ, 0, Q, RING_INC_1110);
	memset(*poly2, 0, L * sizeof *poly2);

	ct->deg = 2;
	ct->ntt = domain;
}

void
swk_folded(Context ctx, Ciphertext *ct, u64 *key, u64 dsw, u64 domain)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods P = CTX_QPE(L, K);
	ring_mods p = ctx.p;

	u64 (*swk)[W][LK][N] = (u64 (*)[W][LK][N])key;
	u64 (*poly0)[N] = (u64 (*)[N])ct->poly;
	u64 (*poly1)[N] = poly0 + L;
	u64 (*poly2)[N] = poly1 + L;
	u64 (*polysw)[N] = (u64 (*)[N])(poly0 + dsw * L);
	u64 ntt = ct->ntt;

	memset(ctx.scratch_poly, 0, ctx.scratch_poly_bytes);
	u64 (*tmp)[N] = (u64 (*)[N])ctx.scratch_poly;
	u64 (*dot0)[N] = tmp + L;
	u64 (*dot1)[N] = dot0 + LK;
	u64 (*ext)[LK][N] = (u64 (*)[LK][N])(dot1 + LK);

	u64 *fbcQ2QP_modQ = (u64 *)ctx.scratch_cnst;
	u64 (*fbcQ2QP_modQP)[LK] = (u64 (*)[LK])(fbcQ2QP_modQ + L);
	u64 (*fbcP2Q_modQ)[L] = (u64 (*)[L])(fbcQ2QP_modQP + L);
	for (u64 i = 0; i < L; ++i) {
		u64 mod = Q.mods[i];
		fbcQ2QP_modQ[i] = ring_coef_inv(prod(i % W, L, W, Q.mods, i, L, mod), mod);
		for (u64 ii = 0; ii < K; ++ii)
			fbcP2Q_modQ[ii][i] = ring_coef_mul(*p.mods, ring_coef_inv(P.mods[ii], mod), mod);
	}
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		for (u64 ii = 0; ii < L; ++ii)
			fbcQ2QP_modQP[ii][i] = prod(ii % W, L, W, Q.mods, ii, L, mod);
	}

	for (u64 i = 0; i < L; ++i)
		memcpy(ext[i % W][i], polysw[i], N * sizeof(u64));
	ring_poly_mulcadd(*polysw, *polysw, fbcQ2QP_modQ, 0, Q, RING_INC_1110);
	if (ntt == 1)
		ring_poly_nttbwd(*polysw, Q);
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W);
		for (u64 ii = k; ii < L; ii += W)
			ring_poly_mulcadd(*ext[k], polysw[ii], fbcQ2QP_modQP[ii], *ext[k], QP, RING_INC_1011);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W && ntt == 1);
		ring_poly_nttfwd(*ext[k], QP);
		for (u64 i = 0; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}

	for (u64 k = 0; k < W; ++k) {
		ring_poly_muladd(*dot0, *ext[k], *swk[0][k], *dot0, QP, RING_INC_1111);
		ring_poly_muladd(*dot1, *ext[k], *swk[1][k], *dot1, QP, RING_INC_1111);
	}

	ring_poly_nttbwd(dot0[L], P);
	ring_poly_nttbwd(dot1[L], P);

	memset(*tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*tmp, dot0[L + ii], fbcP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	if (ntt == 0)
		ring_poly_add(*tmp, *tmp, *poly0, Q, RING_INC_111);
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	if (ntt == 1)
		ring_poly_add(*dot0, *dot0, *poly0, Q, RING_INC_111);
	if (domain == 0)
		ring_poly_nttbwd(*dot0, Q);
	ring_poly_add(*poly0, *dot0, *tmp, Q, RING_INC_111);
	memset(*tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*tmp, dot1[L + ii], fbcP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	if (dsw == 2 && ntt == 0)
		ring_poly_add(*tmp, *tmp, *poly1, Q, RING_INC_111);
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	if (dsw == 2 && ntt == 1)
		ring_poly_add(*dot1, *dot1, *poly1, Q, RING_INC_111);
	if (domain == 0)
		ring_poly_nttbwd(*dot1, Q);
	ring_poly_add(*poly1, *dot1, *tmp, Q, RING_INC_111);
	memset(*poly2, 0, L * sizeof *poly2);

	ct->deg = 2;
	ct->ntt = domain;
}

void
swk_double(Context ctx, Ciphertext *ct, u64 *key, u64 dsw, u64 domain)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 R = ctx.R;
	u64 W = ctx.W;
	u64 W2 = ctx.W2;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods P = CTX_QPE(L, K);
	ring_mods E = CTX_QPE(LK, R);
	ring_mods p = ctx.p;

	u64 (*swk)[W2][W][R][N] = (u64 (*)[W2][W][R][N])key;
	u64 (*poly0)[N] = (u64 (*)[N])ct->poly;
	u64 (*poly1)[N] = poly0 + L;
	u64 (*poly2)[N] = poly1 + L;
	u64 (*polysw)[N] = (u64 (*)[N])(poly0 + dsw * L);
	u64 ntt = ct->ntt;

	memset(ctx.scratch_poly, 0, ctx.scratch_poly_bytes);
	f64 (*tmp)[N] = (f64 (*)[N])ctx.scratch_poly;
	u64 (*dot0)[N] = (u64 (*)[N])(tmp + W2);
	u64 (*dot1)[N] = dot0 + LK;
	u64 (*ext)[R][N] = (u64 (*)[R][N])(dot1 + LK);
	u64 (*decomp0)[R][N] = ext + W;
	u64 (*decomp1)[R][N] = decomp0 + W2;

	u64 *foldE_modQP = (u64 *)ctx.scratch_cnst;
	u64 *fbcQ2E_modQ = foldE_modQP + LK;
	u64 (*fbcQ2E_modE)[R] = (u64 (*)[R])(fbcQ2E_modQ + L);
	u64 (*fbcE2QP_modQP)[LK] = (u64 (*)[LK])(fbcQ2E_modE + L);
	u64 (*fbcP2Q_modQ)[L] = (u64 (*)[L])(fbcE2QP_modQP + R);
	for (u64 i = 0; i < L; ++i) {
		u64 mod = Q.mods[i];
		fbcQ2E_modQ[i] = ring_coef_inv(prod(i % W, L, W, Q.mods, i, L, mod), mod);
		for (u64 ii = 0; ii < K; ++ii)
			fbcP2Q_modQ[ii][i] = ring_coef_mul(*p.mods, ring_coef_inv(P.mods[ii], mod), mod);
	}
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		u64 invP = ring_coef_inv(prod(L, K, 1, QP.mods, i, LK, mod), mod);
		foldE_modQP[i] = ring_coef_mul(invP, mod - prod(0, R, 1, E.mods, (u64)-1, R, mod), mod);
		for (u64 ii = 0; ii < R; ++ii)
			fbcE2QP_modQP[ii][i] = ring_coef_mul(invP, prod(0, R, 1, E.mods, ii, R, mod), mod);
	}
	for (u64 i = L; i < LK; ++i) {
		u64 mod = Q.mods[i];
		u64 neginvp = mod - ring_coef_inv(*p.mods, mod);
		foldE_modQP[i] = ring_coef_mul(foldE_modQP[i], neginvp, mod);
		for (u64 ii = 0; ii < R; ++ii)
			fbcE2QP_modQP[ii][i] = ring_coef_mul(fbcE2QP_modQP[ii][i], neginvp, mod);
	}
	for (u64 i = 0; i < R; ++i) {
		u64 mod = E.mods[i];
		for (u64 ii = 0; ii < L; ++ii)
			fbcQ2E_modE[ii][i] = prod(ii % W, L, W, Q.mods, ii, L, mod);
	}

	ring_poly_mulcadd(*polysw, *polysw, fbcQ2E_modQ, 0, Q, RING_INC_1110);
	if (ntt == 1)
		ring_poly_nttbwd(*polysw, Q);
	for (u64 k = 0; k < W; ++k) {
		for (u64 ii = k; ii < L; ii += W)
			ring_poly_mulcadd(*ext[k], polysw[ii], fbcQ2E_modE[ii], *ext[k], E, RING_INC_1011);
		ring_poly_nttfwd(*ext[k], E);
	}

	memset(decomp0, 0, W2 * sizeof *decomp0);
	memset(decomp1, 0, W2 * sizeof *decomp1);
	for (u64 k2 = 0; k2 < W2; ++k2) {
		for (u64 k = 0; k < W; ++k)
			ring_poly_muladd(*decomp0[k2], *ext[k], *swk[0][k2][k], *decomp0[k2], E, RING_INC_1111);
		ring_poly_nttbwd(*decomp0[k2], E);
		for (u64 k = 0; k < W; ++k)
			ring_poly_muladd(*decomp1[k2], *ext[k], *swk[1][k2][k], *decomp1[k2], E, RING_INC_1111);
		ring_poly_nttbwd(*decomp1[k2], E);
	}

	memset(tmp, 0, W2 * sizeof *tmp);
	for (u64 k2 = 0; k2 < W2; ++k2)
		ring_poly_divcadd(tmp[k2], *decomp0[k2], E.mods, tmp[k2], N, R, RING_INC_0110);
	memset(dot0, 0, sizeof *dot0);
	for (u64 i = 0; i < LK; ++i) {
		ring_mods q = CTX_QPE(i, 1);
		for (u64 ii = 0; ii < R; ++ii)
			ring_poly_mulcadd(dot0[i], decomp0[i % W2][ii], &fbcE2QP_modQP[ii][i], dot0[i], q, RING_INC_1011);
		ring_poly_rndmulcadd(dot0[i], tmp[i % W2], &foldE_modQP[i], dot0[i], q, RING_INC_1011);
	}
	memset(tmp, 0, W2 * sizeof *tmp);
	for (u64 k2 = 0; k2 < W2; ++k2)
		ring_poly_divcadd(tmp[k2], *decomp1[k2], E.mods, tmp[k2], N, R, RING_INC_0110);
	memset(dot1, 0, sizeof *dot1);
	for (u64 i = 0; i < LK; ++i) {
		ring_mods q = CTX_QPE(i, 1);
		for (u64 ii = 0; ii < R; ++ii)
			ring_poly_mulcadd(dot1[i], decomp1[i % W2][ii], &fbcE2QP_modQP[ii][i], dot1[i], q, RING_INC_1011);
		ring_poly_rndmulcadd(dot1[i], tmp[i % W2], &foldE_modQP[i], dot1[i], q, RING_INC_1011);
	}

	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*dot0, dot0[L + ii], fbcP2Q_modQ[ii], *dot0, Q, RING_INC_1011);
	if (ntt == 0 && domain == 0)
		ring_poly_add(*poly0, *dot0, *poly0, Q, RING_INC_111);
	else if (ntt == 0 && domain == 1) {
		ring_poly_add(*poly0, *dot0, *poly0, Q, RING_INC_111);
		ring_poly_nttfwd(*poly0, Q);
	} else if (ntt == 1 && domain == 0) {
		ring_poly_nttbwd(*poly0, Q);
		ring_poly_add(*poly0, *dot0, *poly0, Q, RING_INC_111);
	} else {
		ring_poly_nttfwd(*dot0, Q);
		ring_poly_add(*poly0, *dot0, *poly0, Q, RING_INC_111);
	}
	for (u64 ii = 0; ii < K; ++ii)
		ring_poly_mulcadd(*dot1, dot1[L + ii], fbcP2Q_modQ[ii], *dot1, Q, RING_INC_1011);
	if (dsw == 1) {
		memcpy(*poly1, *dot1, L * sizeof *poly1);
		if (domain == 1)
			ring_poly_nttfwd(*poly1, Q);
	} else if (ntt == 0 && domain == 0) {
		ring_poly_add(*poly1, *dot1, *poly1, Q, RING_INC_111);
	} else if (ntt == 0 && domain == 1) {
		ring_poly_add(*poly1, *dot1, *poly1, Q, RING_INC_111);
		ring_poly_nttfwd(*poly1, Q);
	} else if (ntt == 1 && domain == 0) {
		ring_poly_nttbwd(*poly1, Q);
		ring_poly_add(*poly1, *dot1, *poly1, Q, RING_INC_111);
	} else {
		ring_poly_nttfwd(*dot1, Q);
		ring_poly_add(*poly1, *dot1, *poly1, Q, RING_INC_111);
	}
	memset(*poly2, 0, L * sizeof *poly2);

	ct->deg = 2;
	ct->ntt = domain;
}

void
swk_switch(Context ctx, Ciphertext *ct, u64 *key, u64 dsw, u64 domain, u64 mu, u64 kappa)
{
	u64 N = ctx.N;
	u64 L = ctx.L;
	u64 K = ctx.K;
	u64 W = ctx.W;
	u64 LK = ctx.LK;
	ring_mods QP = CTX_QPE(0, LK);
	ring_mods Q = CTX_QPE(0, L);
	ring_mods Qsw = CTX_QPE(L - mu, mu);
	ring_mods QswP = CTX_QPE(L - mu, mu + K);
	ring_mods p = ctx.p;

	u64 (*swk)[W][LK][N] = (u64 (*)[W][LK][N])key;
	u64 (*poly0)[N] = (u64 (*)[N])ct->poly;
	u64 (*poly1)[N] = poly0 + L;
	u64 (*poly2)[N] = poly1 + L;
	u64 (*polysw)[N] = (u64 (*)[N])(poly0 + dsw * L);
	u64 ntt = ct->ntt;

	memset(ctx.scratch_poly, 0, ctx.scratch_poly_bytes);
	u64 (*tmp)[N] = (u64 (*)[N])ctx.scratch_poly;
	u64 (*dot0)[N] = tmp + L;
	u64 (*dot1)[N] = dot0 + LK;
	u64 (*ext)[LK][N] = (u64 (*)[LK][N])(dot1 + LK);

	u64 *qscale = (u64 *)ctx.scratch_cnst;
	u64 *fbcQ2QP_modQ = qscale + LK;
	u64 (*fbcQ2QP_modQP)[LK] = (u64 (*)[LK])(fbcQ2QP_modQ + L);
	u64 (*fbcQswP2Q_modQ)[L] = (u64 (*)[L])(fbcQ2QP_modQP + L);
	for (u64 i = 0; i < L; ++i) {
		u64 mod = Q.mods[i];
		fbcQ2QP_modQ[i] = ring_coef_inv(prod(i % W, L, W, Q.mods, i, L, mod), mod);
		for (u64 ii = 0; ii < mu + K; ++ii)
			fbcQswP2Q_modQ[ii][i] = ring_coef_mul(ring_coef_inv(QswP.mods[ii], mod), *p.mods, mod);
	}
	for (u64 i = 0; i < LK; ++i) {
		u64 mod = QP.mods[i];
		for (u64 ii = 0; ii < L; ++ii)
			fbcQ2QP_modQP[ii][i] = prod(ii % W, L, W, Q.mods, ii, L, mod);
		u64 qinv = ring_coef_inv(prod(L - mu, L, 1, Q.mods, i, L, mod), mod);
		qscale[i] = ring_coef_mul(qinv, prod(0, kappa, 1, Q.mods, (u64)-1, L, mod), mod);
	}
	for (u64 i = L - mu; i < L; ++i) {
		u64 mod = Q.mods[i];
		u64 neginvp = mod - ring_coef_inv(*p.mods, mod);
		qscale[i] = ring_coef_mul(qscale[i], neginvp, mod);
	}

	for (u64 i = 0; i < LK; ++i)
		QP.mods[i] |= RING_MOD_NOP(i < kappa);
	for (u64 i = kappa; i < L; ++i)
		memcpy(ext[i % W][i], polysw[i], N * sizeof(u64));
	ring_poly_mulcadd(*polysw, *polysw, fbcQ2QP_modQ, 0, Q, RING_INC_1110);
	if (ntt == 1) {
		memcpy(dot0, poly0, L * sizeof *dot0);
		ring_poly_nttbwd(*polysw, Q);
	}
	if (dsw == 2 && ntt == 1)
		memcpy(dot1, poly1, L * sizeof *dot1);
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = kappa; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W);
		for (u64 ii = k; ii < L; ii += W)
			ring_poly_mulcadd(*ext[k], polysw[ii], fbcQ2QP_modQP[ii], *ext[k], QP, RING_INC_1011);
		for (u64 i = kappa; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}
	for (u64 k = 0; k < W; ++k) {
		for (u64 i = kappa; i < LK; ++i)
			QP.mods[i] |= RING_MOD_NOP(i < L && k == i % W && ntt == 1);
		ring_poly_nttfwd(*ext[k], QP);
		for (u64 i = kappa; i < LK; ++i)
			QP.mods[i] &= RING_MOD_USE;
	}
	for (u64 k = 0; k < W; ++k) {
		ring_poly_muladd(*dot0, *ext[k], *swk[0][k], *dot0, QP, RING_INC_1111);
		ring_poly_muladd(*dot1, *ext[k], *swk[1][k], *dot1, QP, RING_INC_1111);
	}
	for (u64 i = 0; i < LK; ++i)
		QP.mods[i] &= RING_MOD_USE;

	ring_poly_nttbwd(dot0[L - mu], QswP);
	if (ntt == 0)
		ring_poly_add(dot0[L - mu], poly0[L - mu], dot0[L - mu], Qsw, RING_INC_111);
	ring_poly_mulcadd(dot0[L - mu], dot0[L - mu], &qscale[L - mu], 0, QswP, RING_INC_1110);
	ring_poly_nttbwd(dot1[L - mu], QswP);
	if (dsw == 2 && ntt == 0)
		ring_poly_add(dot1[L - mu], poly1[L - mu], dot1[L - mu], Qsw, RING_INC_111);
	ring_poly_mulcadd(dot1[L - mu], dot1[L - mu], &qscale[L - mu], 0, QswP, RING_INC_1110);

	memset(tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < mu + K; ++ii)
		ring_poly_mulcadd(*tmp, dot0[L - mu + ii], fbcQswP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	for (u64 i = 0; i < L; ++i)
		Q.mods[i] |= RING_MOD_NOP(i < kappa);
	if (domain == 0)
		ring_poly_nttbwd(*dot0, Q);
	if (ntt == 0 && domain == 0)
		ring_poly_add(*dot0, *dot0, *poly0, Q, RING_INC_111);
	if (ntt == 0 && domain == 1)
		ring_poly_mulcadd(*tmp, *poly0, qscale, *tmp, Q, RING_INC_1111);
	for (u64 i = 0; i < L; ++i)
		Q.mods[i] &= RING_MOD_USE;
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	ring_poly_mulcadd(*poly0, *dot0, qscale, *tmp, Q, RING_INC_1111);
	memset(tmp, 0, L * sizeof *tmp);
	for (u64 ii = 0; ii < mu + K; ++ii)
		ring_poly_mulcadd(*tmp, dot1[L - mu + ii], fbcQswP2Q_modQ[ii], *tmp, Q, RING_INC_1011);
	for (u64 i = 0; i < L; ++i)
		Q.mods[i] |= RING_MOD_NOP(i < kappa);
	if (domain == 0)
		ring_poly_nttbwd(*dot1, Q);
	if (dsw == 2 && ntt == 0 && domain == 0)
		ring_poly_add(*dot1, *dot1, *poly1, Q, RING_INC_111);
	if (dsw == 2 && ntt == 0 && domain == 1)
		ring_poly_mulcadd(*tmp, *poly1, qscale, *tmp, Q, RING_INC_1111);
	for (u64 i = 0; i < L; ++i)
		Q.mods[i] &= RING_MOD_USE;
	if (domain == 1)
		ring_poly_nttfwd(*tmp, Q);
	ring_poly_mulcadd(*poly1, *dot1, qscale, *tmp, Q, RING_INC_1111);
	memset(*poly2, 0, L * sizeof *poly2);

	for (u64 i = 0; i < L - mu; ++i)
		ct->drop[i] = 0;
	for (u64 i = L - mu; i < L; ++i)
		ct->drop[i] = 1;
	ct->deg = 2;
	ct->ntt = domain;
}
