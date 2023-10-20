#include "ksw.h"


tiifhe_Message *
msg_alloc_rand(void)
{
	gmp_randstate_t state;
	tiifhe_Message *m;
	size_t j;

	gmp_randinit_default(state);

	m = tiifhe_msg_alloc(1);
	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_urandomm(m->value[j], state, tiifhe_t.value);
	}

	gmp_randclear(state);

	return m;
}

int
msg_equal(const tiifhe_Message *m1, const tiifhe_Message *m2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		if (mpz_cmp(m1->value[j], m2->value[j]) != 0)
			return 0;
	}

	return 1;
}

void
msg_mul(tiifhe_Message *rop, const tiifhe_Message *op1, const tiifhe_Message *op2)
{
	size_t j;

	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_mul(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], tiifhe_t.value);
	}
}

void
msg_rot1(tiifhe_Message *rop, const tiifhe_Message *m)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	mpz_set(tmp, m->value[0]);
	for (j = 0; j < TIIFHE_N / 2 - 1; ++j) {
		mpz_set(rop->value[j], m->value[j + 1]);
	}
	mpz_set(rop->value[j], tmp);

	mpz_set(tmp, m->value[TIIFHE_N / 2]);
	for (j = TIIFHE_N / 2; j < TIIFHE_N - 1; ++j) {
		mpz_set(rop->value[j], m->value[j + 1]);
	}
	mpz_set(rop->value[j], tmp);

	mpz_clear(tmp);
}

int
main(void)
{
	tiifhe_KeySecret sk;
	tiifhe_KeyPublic *pk;
	tiifhe_Ciphertext *ct, *cpy, *tmp;
	tiifhe_Message *m, *cmp;
	size_t mu, kappa, i;

	tiifhe_bgv_init();
	keyswitch_init();

	pk = tiifhe_bgv_alloc_pk();
	ct = tiifhe_bgv_alloc_ct();
	cpy = tiifhe_bgv_alloc_ct();
	tmp = tiifhe_bgv_alloc_ct();
	m = tiifhe_msg_alloc(1);
	cmp = msg_alloc_rand();

	/* keys and ciphertext */
	printf("[+] seed size: %ld\n", sizeof(tiifhe_Seed));
	puts("[+] computing keys");
	tiifhe_bgv_keygen_secret(&sk);
	tiifhe_bgv_keygen_public(pk, &sk);
	keyswitch_precomp(&sk);

	tiifhe_msg_pack(m, cmp);
	tiifhe_bgv_encrypt(ct, pk, m);

	/* NTT to NTT */
	puts("[+] running NTT to NTT");
	tiifhe_bgv_mul(ct, ct, ct);
	keyswitch_single(ct, 2, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_mul(cmp, cmp, cmp);
	assert(msg_equal(m, cmp));

	/* NTT to coefficient */
	puts("[+] running NTT to coefficient");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to coefficient */
	puts("[+] running coefficient to coefficient");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to NTT */
	puts("[+] running coefficient to NTT");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single(ct, 1, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	puts("[+] folding precomputation");
	keyswitch_precomp_fold();

	/* NTT to NTT (folded) */
	puts("[+] running NTT to NTT (folded)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single_folded(ct, 1, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* NTT to coefficient (folded) */
	puts("[+] running NTT to coefficient (folded)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single_folded(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to coefficient (folded) */
	puts("[+] running coefficient to coefficient (folded)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single_folded(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to NTT (folded) */
	puts("[+] running coefficient to NTT (folded)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_single_folded(ct, 1, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	msg_rot1(cmp, cmp);
	mu = 1, kappa = 2;

	/* NTT to NTT (switch 1/2) */
	puts("[+] running NTT to NTT (switch 1/2)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 1, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	/* NTT to coefficient (switch 1/2) */
	puts("[+] running NTT to coefficient (switch 1/2)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 0, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	for (i = 0; i < TIIFHE_QLEN; ++i) {
		assert(ct[i].ntt == 1);
		tiifhe_poly_intt(i, &ct[i].poly[0]);
		tiifhe_poly_intt(i, &ct[i].poly[1]);
		ct[i].ntt = 0;
	}

	/* coefficient to coefficient (switch 1/2) */
	puts("[+] running coefficient to coefficient (switch 1/2)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 0, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	/* coefficient to NTT (switch 1/2) */
	puts("[+] running coefficient to NTT (switch 1/2)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 1, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	mu = 2, kappa = 3;
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		assert(ct[i].ntt == 0);
		tiifhe_poly_ntt(i, &ct[i].poly[0]);
		tiifhe_poly_ntt(i, &ct[i].poly[1]);
		ct[i].ntt = 1;
	}

	/* NTT to NTT (switch 2/3) */
	puts("[+] running NTT to NTT (switch 2/3)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 1, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	/* NTT to coefficient (switch 2/3) */
	puts("[+] running NTT to coefficient (switch 2/3)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 0, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	for (i = kappa; i < TIIFHE_QLEN; ++i) {
		assert(ct[i].ntt == 1);
		tiifhe_poly_intt(i, &ct[i].poly[0]);
		tiifhe_poly_intt(i, &ct[i].poly[1]);
		ct[i].ntt = 0;
	}

	/* coefficient to coefficient (switch 2/3) */
	puts("[+] running coefficient to coefficient (switch 2/3)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 0, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);
	assert(msg_equal(m, cmp));

	/* coefficient to NTT (switch 2/3) */
	puts("[+] running coefficient to NTT (switch 2/3)");
	tiifhe_bgv_copy(tmp, ct);
	for (i = 0; i < kappa; ++i)
		tiifhe_bgv_modswitch(tmp);
	tiifhe_bgv_rot(tmp, tmp, 1);
	keyswitch_single_switch(tmp, 1, 1, mu, kappa);

	tiifhe_bgv_copy(cpy, tmp);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);

	tiifhe_msg_pack(m, cmp);
	tiifhe_bgv_encrypt(ct, pk, m);

	/* NTT to NTT (double) */
	puts("[+] running NTT to NTT (double)");
	tiifhe_bgv_mul(ct, ct, ct);
	keyswitch_double(ct, 2, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, cpy);
	tiifhe_msg_unpack(m, m);

	msg_mul(cmp, cmp, cmp);
	assert(msg_equal(m, cmp));

	/* NTT to coefficient (double) */
	puts("[+] running NTT to coefficient (double)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_double(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to coefficient (double) */
	puts("[+] running coefficient to coefficient (double)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_double(ct, 1, 0, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	/* coefficient to NTT (double) */
	puts("[+] running coefficient to NTT (double)");
	tiifhe_bgv_rot(ct, ct, 1);
	keyswitch_double(ct, 1, 1, 0);

	tiifhe_bgv_copy(cpy, ct);
	tiifhe_bgv_decrypt(m, &sk, ct);
	tiifhe_msg_unpack(m, m);

	msg_rot1(cmp, cmp);
	assert(msg_equal(m, cmp));

	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct);
	tiifhe_bgv_dealloc(cpy);
	tiifhe_bgv_dealloc(tmp);
	tiifhe_msg_dealloc(m, 1);
	tiifhe_msg_dealloc(cmp, 1);

	tiifhe_bgv_deinit();

	return 0;
}
