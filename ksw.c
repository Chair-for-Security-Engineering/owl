#include "ksw.h"

__extension__ typedef unsigned __int128 uint128_t;

struct precomp {
	tiifhe_Digit negextQ[TIIFHE_QPELEN][TIIFHE_QPELEN];
	tiifhe_Digit neginvt[TIIFHE_QPELEN];
	tiifhe_Digit invP[TIIFHE_QPELEN];
	tiifhe_Digit foldE[TIIFHE_QPELEN];
} precomp;

struct precomp_ext {
	tiifhe_Digit mod[TIIFHE_QPELEN][TIIFHE_QPELEN];
	tiifhe_Digit div[TIIFHE_QPELEN][TIIFHE_QPELEN];
	tiifhe_Digit inv[TIIFHE_QPELEN];
} extQ, extP, extQP, extE;

tiifhe_Poly (*ksw_mul)[2][TIIFHE_OMEGA];
tiifhe_Poly (*ksw_rot)[2][TIIFHE_OMEGA];
tiifhe_Poly (*ksw_decomp_mul)[TIIFHE_OMEGA2][2][TIIFHE_OMEGA];
tiifhe_Poly (*ksw_decomp_rot)[TIIFHE_OMEGA2][2][TIIFHE_OMEGA];

static void
ksw_decomp(tiifhe_Poly (*decomp)[TIIFHE_OMEGA2][2][TIIFHE_OMEGA], tiifhe_Poly (*ksw)[2][TIIFHE_OMEGA], size_t poly)
{
	tiifhe_Poly (*inv)[TIIFHE_OMEGA], tmp;
	tiifhe_PolyApprox *approx;
	size_t i, ii, k, k2;

	inv = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *inv);
	approx = tiifhe_util_alloc(TIIFHE_OMEGA2, sizeof *approx);

	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			tiifhe_poly_mulc(i, &inv[i][k], &ksw[i][poly][k], extQP.inv[i]);
			tiifhe_poly_intt(i, &inv[i][k]);
		}
	}
	for (i = TIIFHE_QPLEN; i < TIIFHE_QPELEN; ++i) {
		size_t idx = i - TIIFHE_QPLEN;
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			memset(approx, 0, TIIFHE_OMEGA2 * sizeof *approx);
			for (ii = 0; ii < TIIFHE_QPLEN; ++ii) {
				k2 = ii % TIIFHE_OMEGA2;
				tiifhe_poly_mulcadd(i, &decomp[idx][k2][poly][k], &inv[ii][k], extQP.div[i][ii], &decomp[idx][k2][poly][k]);
				tiifhe_poly_approx_add(ii, &approx[k2], &inv[ii][k], &approx[k2]);
			}
			for (k2 = 0; k2 < TIIFHE_OMEGA2; ++k2) {
				tiifhe_poly_approx_mulc(i, &tmp, &approx[k2], tiifhe_mod_neg(i, extQP.mod[i][k2]));
				tiifhe_poly_add(i, &decomp[idx][k2][poly][k], &decomp[idx][k2][poly][k], &tmp);
				tiifhe_poly_ntt(i, &decomp[idx][k2][poly][k]);
			}
		}
	}

	tiifhe_util_dealloc(inv);
	tiifhe_util_dealloc(approx);
}

void
keyswitch_double(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0)
{
	tiifhe_Poly (*ksw)[TIIFHE_OMEGA2][2][TIIFHE_OMEGA];
	tiifhe_Poly (*decomp0)[TIIFHE_OMEGA2], (*decomp1)[TIIFHE_OMEGA2];
	tiifhe_Poly (*ext)[TIIFHE_OMEGA], *dot0, *dot1, tmp;
	tiifhe_PolyApprox *approx;
	size_t i, ii, k, k2;
	int ntt;

	ntt = ct[idx0].ntt;
	ksw = dsw == 2 ? ksw_decomp_mul : ksw_decomp_rot;

	approx = tiifhe_util_alloc(TIIFHE_OMEGA2, sizeof *approx);
	decomp0 = tiifhe_util_alloc(TIIFHE_QPELEN, sizeof *decomp0);
	decomp1 = tiifhe_util_alloc(TIIFHE_QPELEN, sizeof *decomp1);
	dot0 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot0);
	dot1 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot1);
	ext = tiifhe_util_alloc(TIIFHE_QPELEN, sizeof *ext);

	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		tiifhe_poly_mulc(i, &ct[i].poly[dsw], &ct[i].poly[dsw], extQ.inv[i]);
		if (ntt == 1)
			tiifhe_poly_intt(i, &ct[i].poly[dsw]);
	}
	for (i = TIIFHE_QPLEN; i < TIIFHE_QPELEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			for (ii = k; ii < TIIFHE_QLEN; ii += TIIFHE_OMEGA)
				tiifhe_poly_mulcadd(i, &ext[i][k], &ct[ii].poly[dsw], extQ.div[i][ii], &ext[i][k]);
			tiifhe_poly_ntt(i, &ext[i][k]);
		}
	}

	for (i = TIIFHE_QPLEN; i < TIIFHE_QPELEN; ++i) {
		for (k2 = 0; k2 < TIIFHE_OMEGA2; ++k2) {
			memset(&decomp0[i][k2], 0, sizeof decomp0[i][k2]);
			for (k = 0; k < TIIFHE_OMEGA; ++k)
				tiifhe_poly_muladd(i, &decomp0[i][k2], &ext[i][k], &ksw[i - TIIFHE_QPLEN][k2][0][k], &decomp0[i][k2]);
			tiifhe_poly_intt(i, &decomp0[i][k2]);

			memset(&decomp1[i][k2], 0, sizeof decomp1[i][k2]);
			for (k = 0; k < TIIFHE_OMEGA; ++k)
				tiifhe_poly_muladd(i, &decomp1[i][k2], &ext[i][k], &ksw[i - TIIFHE_QPLEN][k2][1][k], &decomp1[i][k2]);
			tiifhe_poly_intt(i, &decomp1[i][k2]);
		}
	}

	memset(approx, 0, TIIFHE_OMEGA2 * sizeof *approx);
	for (k2 = 0; k2 < TIIFHE_OMEGA2; ++k2) {
		for (ii = TIIFHE_QPLEN; ii < TIIFHE_QPELEN; ++ii)
			tiifhe_poly_approx_add(ii, &approx[k2], &decomp0[ii][k2], &approx[k2]);
	}
	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		memset(&dot0[i], 0, sizeof dot0[i]);
		for (ii = TIIFHE_QPLEN; ii < TIIFHE_QPELEN; ++ii)
			tiifhe_poly_mulcadd(i, &dot0[i], &decomp0[ii][i % TIIFHE_OMEGA2], extE.div[i][ii], &dot0[i]);
		tiifhe_poly_approx_mulc(i, &tmp, &approx[i % TIIFHE_OMEGA2], precomp.foldE[i]);
		tiifhe_poly_add(i, &dot0[i], &dot0[i], &tmp);
	}

	memset(approx, 0, TIIFHE_OMEGA2 * sizeof *approx);
	for (k2 = 0; k2 < TIIFHE_OMEGA2; ++k2) {
		for (ii = TIIFHE_QPLEN; ii < TIIFHE_QPELEN; ++ii)
			tiifhe_poly_approx_add(ii, &approx[k2], &decomp1[ii][k2], &approx[k2]);
	}
	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		memset(&dot1[i], 0, sizeof dot1[i]);
		for (ii = TIIFHE_QPLEN; ii < TIIFHE_QPELEN; ++ii)
			tiifhe_poly_mulcadd(i, &dot1[i], &decomp1[ii][i % TIIFHE_OMEGA2], extE.div[i][ii], &dot1[i]);
		tiifhe_poly_approx_mulc(i, &tmp, &approx[i % TIIFHE_OMEGA2], precomp.foldE[i]);
		tiifhe_poly_add(i, &dot1[i], &dot1[i], &tmp);
	}

	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &dot0[i], &dot0[ii], extP.div[i][ii], &dot0[i]);
		if (ntt == 0 && domain == 0)
			tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &ct[i].poly[0]);
		else if (ntt == 0 && domain == 1) {
			tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &ct[i].poly[0]);
			tiifhe_poly_ntt(i, &ct[i].poly[0]);
		} else if (ntt == 1 && domain == 0) {
			tiifhe_poly_intt(i, &ct[i].poly[0]);
			tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &ct[i].poly[0]);
		} else {
			tiifhe_poly_ntt(i, &dot0[i]);
			tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &ct[i].poly[0]);
		}

		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &dot1[i], &dot1[ii], extP.div[i][ii], &dot1[i]);
		if (dsw == 1) {
			memcpy(&ct[i].poly[1], &dot1[i], sizeof ct[i].poly[1]);
			if (domain == 1)
				tiifhe_poly_ntt(i, &ct[i].poly[1]);
		} else if (ntt == 0 && domain == 0) {
			tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &ct[i].poly[1]);
		} else if (ntt == 0 && domain == 1) {
			tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &ct[i].poly[1]);
			tiifhe_poly_ntt(i, &ct[i].poly[1]);
		} else if (ntt == 1 && domain == 0) {
			tiifhe_poly_intt(i, &ct[i].poly[1]);
			tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &ct[i].poly[1]);
		} else {
			tiifhe_poly_ntt(i, &dot1[i]);
			tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &ct[i].poly[1]);
		}

		memset(&ct[i].poly[2], 0, sizeof ct[i].poly[2]);
		ct[i].degree = 2;
		ct[i].ntt = domain;
	}

	tiifhe_util_dealloc(approx);
	tiifhe_util_dealloc(decomp0);
	tiifhe_util_dealloc(decomp1);
	tiifhe_util_dealloc(dot0);
	tiifhe_util_dealloc(dot1);
	tiifhe_util_dealloc(ext);
}

void
keyswitch_single(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0)
{
	tiifhe_Poly (*ksw)[2][TIIFHE_OMEGA], (*ext)[TIIFHE_OMEGA];
	tiifhe_Poly *dot0, *dot1, tmp;
	size_t i, ii, k;
	int ntt;

	ntt = ct[idx0].ntt;
	ksw = dsw == 2 ? ksw_mul : ksw_rot;

	dot0 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot0);
	dot1 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot1);
	ext = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ext);

	/* base extension */
	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		k = i % TIIFHE_OMEGA;
		memcpy(&ext[i][k], &ct[i].poly[dsw], sizeof ext[i][k]);
		tiifhe_poly_mulc(i, &ct[i].poly[dsw], &ct[i].poly[dsw], extQ.inv[i]);
		if (ntt == 1)
			tiifhe_poly_intt(i, &ct[i].poly[dsw]);
	}
	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			for (ii = k; ii < idx0; ii += TIIFHE_OMEGA)
				;
			for (; ii < TIIFHE_QLEN; ii += TIIFHE_OMEGA)
				tiifhe_poly_mulcadd(i, &ext[i][k], &ct[ii].poly[dsw], extQ.div[i][ii], &ext[i][k]);
		}
	}

	/* dot product */
	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (ntt == 1 && i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			tiifhe_poly_ntt(i, &ext[i][k]);
		}

		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot0[i], &ext[i][k], &ksw[i][0][k], &dot0[i]);
		if (ntt == 1)
			tiifhe_poly_mulcadd(i, &dot0[i], &ct[i].poly[0], tiifhe_q[i].P, &dot0[i]);

		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot1[i], &ext[i][k], &ksw[i][1][k], &dot1[i]);
		if (dsw == 2 && ntt == 1)
			tiifhe_poly_mulcadd(i, &dot1[i], &ct[i].poly[1], tiifhe_q[i].P, &dot1[i]);
	}

	/* delta computation */
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_intt(i, &dot0[i]);
		tiifhe_poly_mulc(i, &dot0[i], &dot0[i], precomp.neginvt[i]);

		tiifhe_poly_intt(i, &dot1[i]);
		tiifhe_poly_mulc(i, &dot1[i], &dot1[i], precomp.neginvt[i]);
	}

	/* modulus switching */
	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_mulc(i, &dot0[i], &dot0[i], extP.inv[i]);
		tiifhe_poly_mulc(i, &dot1[i], &dot1[i], extP.inv[i]);
	}
	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot0[ii], extP.div[i][ii], &tmp);
		tiifhe_poly_mulc(i, &tmp, &tmp, tiifhe_q[i].t);
		if (ntt == 0)
			tiifhe_poly_mulcadd(i, &tmp, &ct[i].poly[0], tiifhe_q[i].P, &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		else
			tiifhe_poly_intt(i, &dot0[i]);
		tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &tmp);
		tiifhe_poly_mulc(i, &ct[i].poly[0], &ct[i].poly[0], precomp.invP[i]);

		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot1[ii], extP.div[i][ii], &tmp);
		tiifhe_poly_mulc(i, &tmp, &tmp, tiifhe_q[i].t);
		if (dsw == 2 && ntt == 0)
			tiifhe_poly_mulcadd(i, &tmp, &ct[i].poly[1], tiifhe_q[i].P, &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		else
			tiifhe_poly_intt(i, &dot1[i]);
		tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &tmp);
		tiifhe_poly_mulc(i, &ct[i].poly[1], &ct[i].poly[1], precomp.invP[i]);

		memset(&ct[i].poly[2], 0, sizeof ct[i].poly[2]);
		ct[i].degree = 2;
		ct[i].ntt = domain;
	}

	tiifhe_util_dealloc(dot0);
	tiifhe_util_dealloc(dot1);
	tiifhe_util_dealloc(ext);
}

void
keyswitch_single_folded(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0)
{
	tiifhe_Poly (*ksw)[2][TIIFHE_OMEGA], (*ext)[TIIFHE_OMEGA];
	tiifhe_Poly *dot0, *dot1, tmp;
	size_t i, ii, k;
	int ntt;

	ntt = ct[idx0].ntt;
	ksw = dsw == 2 ? ksw_mul : ksw_rot;

	dot0 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot0);
	dot1 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot1);
	ext = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ext);

	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		k = i % TIIFHE_OMEGA;
		memcpy(&ext[i][k], &ct[i].poly[dsw], sizeof ext[i][k]);
		tiifhe_poly_mulc(i, &ct[i].poly[dsw], &ct[i].poly[dsw], extQ.inv[i]);
		if (ntt == 1)
			tiifhe_poly_intt(i, &ct[i].poly[dsw]);
	}
	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			for (ii = k; ii < idx0; ii += TIIFHE_OMEGA)
				;
			for (; ii < TIIFHE_QLEN; ii += TIIFHE_OMEGA)
				tiifhe_poly_mulcadd(i, &ext[i][k], &ct[ii].poly[dsw], extQ.div[i][ii], &ext[i][k]);
		}
	}

	for (i = idx0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (ntt == 1 && i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			tiifhe_poly_ntt(i, &ext[i][k]);
		}

		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot0[i], &ext[i][k], &ksw[i][0][k], &dot0[i]);
		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot1[i], &ext[i][k], &ksw[i][1][k], &dot1[i]);
	}

	for (i = TIIFHE_QLEN; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_intt(i, &dot0[i]);
		tiifhe_poly_intt(i, &dot1[i]);
	}

	for (i = idx0; i < TIIFHE_QLEN; ++i) {
		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot0[ii], extP.div[i][ii], &tmp);
		if (ntt == 0)
			tiifhe_poly_add(i, &tmp, &tmp, &ct[i].poly[0]);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		if (ntt == 1)
			tiifhe_poly_add(i, &dot0[i], &dot0[i], &ct[i].poly[0]);
		if (domain == 0)
			tiifhe_poly_intt(i, &dot0[i]);
		tiifhe_poly_add(i, &ct[i].poly[0], &dot0[i], &tmp);

		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot1[ii], extP.div[i][ii], &tmp);
		if (dsw == 2 && ntt == 0)
			tiifhe_poly_add(i, &tmp, &tmp, &ct[i].poly[1]);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		if (dsw == 2 && ntt == 1)
			tiifhe_poly_add(i, &dot1[i], &dot1[i], &ct[i].poly[1]);
		if (domain == 0)
			tiifhe_poly_intt(i, &dot1[i]);
		tiifhe_poly_add(i, &ct[i].poly[1], &dot1[i], &tmp);

		memset(&ct[i].poly[2], 0, sizeof ct[i].poly[2]);
		ct[i].degree = 2;
		ct[i].ntt = domain;
	}

	tiifhe_util_dealloc(dot0);
	tiifhe_util_dealloc(dot1);
	tiifhe_util_dealloc(ext);
}

void
keyswitch_single_switch(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t mu, size_t kappa)
{

	tiifhe_Poly (*ksw)[2][TIIFHE_OMEGA], (*ext)[TIIFHE_OMEGA];
	tiifhe_Poly *dot0, *dot1, tmp;
	tiifhe_Digit qswitch[TIIFHE_QPLEN];
	struct precomp_ext extmuP;
	size_t i, ii, k;
	int ntt;

	ntt = ct[kappa].ntt;
	ksw = dsw == 2 ? ksw_mul : ksw_rot;

	dot0 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot0);
	dot1 = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *dot1);
	ext = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ext);

	/* precomp */
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		extmuP.inv[i] = i >= TIIFHE_QLEN - mu && i < TIIFHE_QLEN ? precomp.neginvt[i] : 1;
		for (ii = TIIFHE_QLEN - mu; ii < TIIFHE_QLEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			extmuP.inv[i] = tiifhe_mod_mul(i, extmuP.inv[i], inv);
			extmuP.div[i][ii] = tiifhe_mod_mul(i, tiifhe_q[i].t, inv);
		}
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii)
			extmuP.div[i][ii] = extP.div[i][ii];

		qswitch[i] = 1;
		for (ii = 0; ii < kappa; ++ii)
			qswitch[i] = tiifhe_mod_mul(i, qswitch[i], tiifhe_q[ii].value);
		qswitch[i] = tiifhe_mod_mul(i, qswitch[i], extmuP.inv[i]);
	}

	for (i = kappa; i < TIIFHE_QLEN; ++i) {
		k = i % TIIFHE_OMEGA;
		memcpy(&ext[i][k], &ct[i].poly[dsw], sizeof ext[i][k]);
		tiifhe_poly_mulc(i, &ct[i].poly[dsw], &ct[i].poly[dsw], extQ.inv[i]);
		if (ntt == 1) {
			memcpy(&dot0[i], &ct[i].poly[0], sizeof dot0[i]);
			if (dsw == 2)
				memcpy(&dot1[i], &ct[i].poly[1], sizeof dot1[i]);
			tiifhe_poly_intt(i, &ct[i].poly[dsw]);
		}
	}
	for (i = kappa; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			for (ii = k; ii < kappa; ii += TIIFHE_OMEGA)
				;
			for (; ii < TIIFHE_QLEN; ii += TIIFHE_OMEGA)
				tiifhe_poly_mulcadd(i, &ext[i][k], &ct[ii].poly[dsw], extQ.div[i][ii], &ext[i][k]);
		}
	}

	for (i = kappa; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			if (ntt == 1 && i < TIIFHE_QLEN && k == i % TIIFHE_OMEGA)
				continue;

			tiifhe_poly_ntt(i, &ext[i][k]);
		}

		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot0[i], &ext[i][k], &ksw[i][0][k], &dot0[i]);
		for (k = 0; k < TIIFHE_OMEGA; ++k)
			tiifhe_poly_muladd(i, &dot1[i], &ext[i][k], &ksw[i][1][k], &dot1[i]);
	}

	/* delta computation */
	for (i = TIIFHE_QLEN - mu; i < TIIFHE_QPLEN; ++i) {
		tiifhe_poly_intt(i, &dot0[i]);
		if (i < TIIFHE_QLEN && ntt == 0)
			tiifhe_poly_add(i, &dot0[i], &ct[i].poly[0], &dot0[i]);
		tiifhe_poly_mulc(i, &dot0[i], &dot0[i], qswitch[i]);

		tiifhe_poly_intt(i, &dot1[i]);
		if (dsw == 2 && i < TIIFHE_QLEN && ntt == 0)
			tiifhe_poly_add(i, &dot1[i], &ct[i].poly[1], &dot1[i]);
		tiifhe_poly_mulc(i, &dot1[i], &dot1[i], qswitch[i]);
	}

	/* modulus switching */
	for (i = 0; i < kappa; ++i) {
		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN - mu; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot0[ii], extmuP.div[i][ii], &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		tiifhe_poly_mulcadd(i, &ct[i].poly[0], &dot0[i], qswitch[i], &tmp);

		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN - mu; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot1[ii], extmuP.div[i][ii], &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		tiifhe_poly_mulcadd(i, &ct[i].poly[1], &dot1[i], qswitch[i], &tmp);

		memset(&ct[i].poly[2], 0, sizeof ct[i].poly[2]);
		ct[i].degree = 2;
		ct[i].drop = 0;
		ct[i].ntt = domain;
	}
	for (i = kappa; i < TIIFHE_QLEN - mu; ++i) {
		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN - mu; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot0[ii], extmuP.div[i][ii], &tmp);
		if (domain == 0)
			tiifhe_poly_intt(i, &dot0[i]);
		if (ntt == 0 && domain == 0)
			tiifhe_poly_add(i, &dot0[i], &dot0[i], &ct[i].poly[0]);
		if (ntt == 0 && domain == 1)
			tiifhe_poly_mulcadd(i, &tmp, &ct[i].poly[0], qswitch[i], &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		tiifhe_poly_mulcadd(i, &ct[i].poly[0], &dot0[i], qswitch[i], &tmp);

		memset(&tmp, 0, sizeof tmp);
		for (ii = TIIFHE_QLEN - mu; ii < TIIFHE_QPLEN; ++ii)
			tiifhe_poly_mulcadd(i, &tmp, &dot1[ii], extmuP.div[i][ii], &tmp);
		if (domain == 0)
			tiifhe_poly_intt(i, &dot1[i]);
		if (dsw == 2 && ntt == 0 && domain == 0)
			tiifhe_poly_add(i, &dot1[i], &dot1[i], &ct[i].poly[1]);
		if (dsw == 2 && ntt == 0 && domain == 1)
			tiifhe_poly_mulcadd(i, &tmp, &ct[i].poly[1], qswitch[i], &tmp);
		if (domain == 1)
			tiifhe_poly_ntt(i, &tmp);
		tiifhe_poly_mulcadd(i, &ct[i].poly[1], &dot1[i], qswitch[i], &tmp);

		memset(&ct[i].poly[2], 0, sizeof ct[i].poly[2]);
		ct[i].degree = 2;
		ct[i].drop = 0;
		ct[i].ntt = domain;
	}
	for (i = TIIFHE_QLEN - mu; i < TIIFHE_QLEN; ++i) {
		ct[i].degree = 2;
		ct[i].drop = 1;
	}

	tiifhe_util_dealloc(dot0);
	tiifhe_util_dealloc(dot1);
	tiifhe_util_dealloc(ext);
}

void
keyswitch_init(void)
{

	ksw_mul = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ksw_mul);
	ksw_rot = tiifhe_util_alloc(TIIFHE_QPLEN, sizeof *ksw_rot);
	ksw_decomp_mul = tiifhe_util_alloc(TIIFHE_ELEN, sizeof *ksw_decomp_mul);
	ksw_decomp_rot = tiifhe_util_alloc(TIIFHE_ELEN, sizeof *ksw_decomp_rot);
}

void
keyswitch_precomp(const tiifhe_KeySecret *sk)
{
	tiifhe_KeySwitch *kmul, *krot;
	size_t i, ii, k;

	kmul = tiifhe_bgv_alloc_ksw();
	krot = tiifhe_bgv_alloc_ksw();

	/* precomp */
	for (i = 0; i < TIIFHE_QPELEN; ++i) {
		precomp.neginvt[i] = tiifhe_mod_neg(i, tiifhe_mod_inv(i, tiifhe_q[i].t));
		precomp.invP[i] = tiifhe_mod_inv(i, tiifhe_q[i].P);

		for (k = 0; k < TIIFHE_OMEGA; ++k)
			extQ.mod[i][k] = 1;
		for (ii = 0; ii < TIIFHE_QLEN; ++ii) {
			k = ii % TIIFHE_OMEGA;
			extQ.mod[i][k] = tiifhe_mod_mul(i, extQ.mod[i][k], tiifhe_q[ii].value);
		}
		extQ.inv[i] = 1;
		for (ii = 0; ii < TIIFHE_QLEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			if (i % TIIFHE_OMEGA == ii % TIIFHE_OMEGA)
				extQ.inv[i] = tiifhe_mod_mul(i, extQ.inv[i], inv);
			extQ.div[i][ii] = tiifhe_mod_mul(i, extQ.mod[i][ii % TIIFHE_OMEGA], inv);
		}

		extP.inv[i] = 1;
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			extP.inv[i] = tiifhe_mod_mul(i, extP.inv[i], inv);
			extP.div[i][ii] = tiifhe_mod_mul(i, tiifhe_q[i].P, inv);
		}
		extP.div[i][i] = tiifhe_mod_inv(i, extP.inv[i]);
	}

	tiifhe_bgv_keygen_switch_mul(kmul, sk);
	tiifhe_bgv_keygen_switch_rot(krot, sk, 1);
	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			memcpy(&ksw_mul[i][0][k], &kmul[i].poly[k], sizeof ksw_mul[i][0][k]);
			tiifhe_poly_sample_uniform(i, &ksw_mul[i][1][k], &kmul[i].seed);

			memcpy(&ksw_rot[i][0][k], &krot[i].poly[k], sizeof ksw_rot[i][0][k]);
			tiifhe_poly_sample_uniform(i, &ksw_rot[i][1][k], &krot[i].seed);
		}
		kmul[i].seed.init = 1;
		krot[i].seed.init = 1;
	}

	tiifhe_util_dealloc(kmul);
	tiifhe_util_dealloc(krot);
}

void
keyswitch_precomp_fold(void)
{
	size_t i, ii, k, k2;

	for (i = 0; i < TIIFHE_QPELEN; ++i) {
		extP.inv[i] = i >= TIIFHE_QLEN ? precomp.neginvt[i] : 1;
		for (ii = TIIFHE_QLEN; ii < TIIFHE_QPLEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			extP.inv[i] = tiifhe_mod_mul(i, extP.inv[i], inv);
			extP.div[i][ii] = tiifhe_mod_mul(i, tiifhe_q[i].t, inv);
		}

		for (k = 0; k < TIIFHE_OMEGA2; ++k)
			extQP.mod[i][k] = 1;
		for (ii = 0; ii < TIIFHE_QPLEN; ++ii) {
			k = ii % TIIFHE_OMEGA2;
			extQP.mod[i][k] = tiifhe_mod_mul(i, extQP.mod[i][k], tiifhe_q[ii].value);
		}
		extQP.inv[i] = 1;
		for (ii = 0; ii < TIIFHE_QPLEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			if (i % TIIFHE_OMEGA2 == ii % TIIFHE_OMEGA2)
				extQP.inv[i] = tiifhe_mod_mul(i, extQP.inv[i], inv);
			extQP.div[i][ii] = tiifhe_mod_mul(i, extQP.mod[i][ii % TIIFHE_OMEGA2], inv);
		}

		extE.inv[i] = 1;
		for (ii = TIIFHE_QPLEN; ii < TIIFHE_QPELEN; ++ii) {
			tiifhe_Digit inv = tiifhe_mod_inv(i, tiifhe_q[ii].value);
			extE.inv[i] = tiifhe_mod_mul(i, extE.inv[i], inv);
			extE.div[i][ii] = tiifhe_mod_mul(i, tiifhe_q[i].E, inv);
			extE.div[i][ii] = tiifhe_mod_mul(i, extE.div[i][ii], extP.inv[i]);
		}

		precomp.foldE[i] = tiifhe_mod_neg(i, tiifhe_mod_mul(i, tiifhe_q[i].E, extP.inv[i]));
		for (k = 0; k < TIIFHE_OMEGA; ++k)
			precomp.negextQ[i][k] = tiifhe_mod_neg(i, extQ.mod[i][k]);
	}

	memset(ksw_decomp_mul, 0, TIIFHE_ELEN * sizeof *ksw_decomp_mul);
	memset(ksw_decomp_rot, 0, TIIFHE_ELEN * sizeof *ksw_decomp_rot);
	ksw_decomp(ksw_decomp_mul, ksw_mul, 0);
	ksw_decomp(ksw_decomp_mul, ksw_mul, 1);
	ksw_decomp(ksw_decomp_rot, ksw_rot, 0);
	ksw_decomp(ksw_decomp_rot, ksw_rot, 1);

	for (i = 0; i < TIIFHE_QPLEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			tiifhe_poly_mulc(i, &ksw_mul[i][0][k], &ksw_mul[i][0][k], extP.inv[i]);
			tiifhe_poly_mulc(i, &ksw_mul[i][1][k], &ksw_mul[i][1][k], extP.inv[i]);
			tiifhe_poly_mulc(i, &ksw_rot[i][0][k], &ksw_rot[i][0][k], extP.inv[i]);
			tiifhe_poly_mulc(i, &ksw_rot[i][1][k], &ksw_rot[i][1][k], extP.inv[i]);
		}
	}
	for (; i < TIIFHE_QPELEN; ++i) {
		for (k = 0; k < TIIFHE_OMEGA; ++k) {
			for (k2 = 0; k2 < TIIFHE_OMEGA2; ++k2) {
				size_t idx = i - TIIFHE_QPLEN;
				tiifhe_poly_mulc(i, &ksw_decomp_mul[idx][k2][0][k], &ksw_decomp_mul[idx][k2][0][k], extE.inv[i]);
				tiifhe_poly_mulc(i, &ksw_decomp_mul[idx][k2][1][k], &ksw_decomp_mul[idx][k2][1][k], extE.inv[i]);
				tiifhe_poly_mulc(i, &ksw_decomp_rot[idx][k2][0][k], &ksw_decomp_rot[idx][k2][0][k], extE.inv[i]);
				tiifhe_poly_mulc(i, &ksw_decomp_rot[idx][k2][1][k], &ksw_decomp_rot[idx][k2][1][k], extE.inv[i]);
			}
		}
	}
}
