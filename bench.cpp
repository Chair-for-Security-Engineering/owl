#include <benchmark/benchmark.h>
#include <gmp.h>

#if OWL_FOLD
	#define keyswitch   OWL_KEYSWITCH
	#define precomp(sk) keyswitch_precomp(sk), keyswitch_precomp_fold();
#else
	#define keyswitch   OWL_KEYSWITCH
	#define precomp(sk) keyswitch_precomp(sk);
#endif

extern "C" {
	#include "ksw.h"
}

void
msg_rand(tiifhe_Message *m)
{
	gmp_randstate_t state;
	size_t j;

	gmp_randinit_default(state);
	for (j = 0; j < TIIFHE_N; ++j) {
		mpz_urandomm(m->value[j], state, tiifhe_t.value);
	}
	gmp_randclear(state);
}

void
ct_rand(tiifhe_Ciphertext *ct, tiifhe_KeyPublic *pk, int domain)
{
	tiifhe_Message *m;
	size_t i;

	m = tiifhe_msg_alloc(1);
	msg_rand(m);

	tiifhe_bgv_encrypt(ct, pk, m);
	for (i = 0; i < TIIFHE_QLEN; ++i) {
		assert(ct[i].ntt == 1);
		if (domain)
			continue;

		tiifhe_poly_intt(i, &ct[i].poly[0]);
		tiifhe_poly_intt(i, &ct[i].poly[1]);
		ct[i].ntt = 0;
	}

	tiifhe_msg_dealloc(m, 1);
}

void
BM_intt_to_intt(benchmark::State& state)
{
	tiifhe_KeySecret sk;
	tiifhe_KeyPublic *pk;
	tiifhe_Ciphertext *ct;

	pk = tiifhe_bgv_alloc_pk();
	ct = tiifhe_bgv_alloc_ct();

	tiifhe_bgv_keygen_secret(&sk);
	tiifhe_bgv_keygen_public(pk, &sk);
	keyswitch_precomp(&sk);
	ct_rand(ct, pk, 0);

	for (auto _ : state) {
		state.PauseTiming();
		tiifhe_bgv_rot(ct, ct, 1);

		state.ResumeTiming();
		keyswitch(ct, 1, 0, 0);

		benchmark::DoNotOptimize(ct);
		benchmark::ClobberMemory();
	}

	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct);
}

void
BM_intt_to_ntt(benchmark::State& state)
{
	tiifhe_KeySecret sk;
	tiifhe_KeyPublic *pk;
	tiifhe_Ciphertext *ct;

	pk = tiifhe_bgv_alloc_pk();
	ct = tiifhe_bgv_alloc_ct();

	tiifhe_bgv_keygen_secret(&sk);
	tiifhe_bgv_keygen_public(pk, &sk);
	keyswitch_precomp(&sk);
	ct_rand(ct, pk, 0);

	for (auto _ : state) {
		state.PauseTiming();
		tiifhe_bgv_rot(ct, ct, 1);

		state.ResumeTiming();
		keyswitch(ct, 1, 1, 0);

		benchmark::DoNotOptimize(ct);
		benchmark::ClobberMemory();
	}

	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct);
}

void
BM_ntt_to_intt(benchmark::State& state)
{
	tiifhe_KeySecret sk;
	tiifhe_KeyPublic *pk;
	tiifhe_Ciphertext *ct;

	pk = tiifhe_bgv_alloc_pk();
	ct = tiifhe_bgv_alloc_ct();

	tiifhe_bgv_keygen_secret(&sk);
	tiifhe_bgv_keygen_public(pk, &sk);
	keyswitch_precomp(&sk);
	ct_rand(ct, pk, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiifhe_bgv_rot(ct, ct, 1);

		state.ResumeTiming();
		keyswitch(ct, 1, 0, 0);

		benchmark::DoNotOptimize(ct);
		benchmark::ClobberMemory();
	}

	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct);
}

void
BM_ntt_to_ntt(benchmark::State& state)
{
	tiifhe_KeySecret sk;
	tiifhe_KeyPublic *pk;
	tiifhe_Ciphertext *ct;

	pk = tiifhe_bgv_alloc_pk();
	ct = tiifhe_bgv_alloc_ct();

	tiifhe_bgv_keygen_secret(&sk);
	tiifhe_bgv_keygen_public(pk, &sk);
	keyswitch_precomp(&sk);
	ct_rand(ct, pk, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiifhe_bgv_rot(ct, ct, 1);

		state.ResumeTiming();
		keyswitch(ct, 1, 1, 0);

		benchmark::DoNotOptimize(ct);
		benchmark::ClobberMemory();
	}

	tiifhe_bgv_dealloc(pk);
	tiifhe_bgv_dealloc(ct);
}

BENCHMARK(BM_intt_to_intt);
BENCHMARK(BM_intt_to_ntt);
BENCHMARK(BM_ntt_to_intt);
BENCHMARK(BM_ntt_to_ntt);

int
main(int argc, char** argv)
{

	printf("[+] N:  %d\n", TIIFHE_N);
	printf("[+] t:  %s\n", TIIFHE_T);
	printf("[+] l:  %d\n", TIIFHE_QLEN);
	printf("[+] k:  %d\n", TIIFHE_PLEN);
	printf("[+] r:  %d\n", TIIFHE_ELEN);
	printf("[+] w:  %d\n", TIIFHE_OMEGA);
	printf("[+] w2: %d\n", TIIFHE_OMEGA2);

	tiifhe_bgv_init();
	keyswitch_init();
	::benchmark::Initialize(&argc, argv);
	::benchmark::RunSpecifiedBenchmarks();
}
