extern "C" {
	#include "ring.h"
}
#include "hexl/hexl.hpp"
using namespace intel;

ring_mods
ring_mods_init(u64 *mods, i32 deg, i32 len)
{
	void **p = new void *[len];
	ring_mods rop = { p, mods, deg, len };
	while (len--)
		*p++ = (void *)(new hexl::NTT(deg, *mods++));
	return rop;
}

u64
ring_coef_inv(u64 op, u64 mod)
{
	return hexl::InverseMod(op, mod);
}

u64
ring_coef_mul(u64 op1, u64 op2, u64 mod)
{
	return hexl::MultiplyMod(op1, op2, mod);
}

void
ring_poly_add(u64 *rop, u64 const *op1, u64 const *op2, ring_mods mods, u8 inc)
{
	u64 *p = mods.mods;
	i32 deg = mods.deg, len = mods.len;
	i32 inc0 = ((inc >> 0) & 1) * deg;
	i32 inc1 = ((inc >> 1) & 1) * deg;
	i32 inc2 = ((inc >> 2) & 1) * deg;
	while (len--) {
		u64 mod = *p++;
		if ((s64)mod > 0)
			hexl::EltwiseAddMod(rop, op1, op2, deg, mod);
		rop += inc0, op1 += inc1, op2 += inc2;
	}
}

void
ring_poly_mul(u64 *rop, u64 const *op1, u64 const *op2, ring_mods mods, u8 inc)
{
	u64 *p = mods.mods;
	i32 deg = mods.deg, len = mods.len;
	i32 inc0 = ((inc >> 0) & 1) * deg;
	i32 inc1 = ((inc >> 1) & 1) * deg;
	i32 inc2 = ((inc >> 2) & 1) * deg;
	while (len--) {
		u64 mod = *p++;
		if ((s64)mod > 0)
			hexl::EltwiseMultMod(rop, op1, op2, deg, mod, 1);
		rop += inc0, op1 += inc1, op2 += inc2;
	}
}

void
ring_poly_muladd(u64 *rop, u64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc)
{
	u64 buf[128];
	u64 *p = mods.mods;
	i32 deg = mods.deg, len = mods.len;
	i32 dec0 = ((~inc >> 0) & 1) * deg;
	i32 dec1 = ((~inc >> 1) & 1) * deg;
	i32 dec2 = ((~inc >> 2) & 1) * deg;
	i32 dec3 = ((~inc >> 3) & 1) * deg;
	while (len--) {
		u64 mod = *p++;
		i32 i = deg / 128;
		while (i--) {
			if ((s64)mod > 0) {
				hexl::EltwiseMultMod(buf, op1, op2, 128, mod, 1);
				hexl::EltwiseAddMod(rop, buf, op3, 128, mod);
			}
			rop += 128, op1 += 128, op2 += 128, op3 += 128;
		}
		rop -= dec0, op1 -= dec1, op2 -= dec2, op3 -= dec3;
	}
}

void
ring_poly_mulcadd(u64 *rop, u64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc)
{
	u64 *p = mods.mods;
	i32 deg = mods.deg, len = mods.len;
	i32 inc0 = ((inc >> 0) & 1) * deg;
	i32 inc1 = ((inc >> 1) & 1) * deg;
	i32 inc2 = ((inc >> 2) & 1);
	i32 inc3 = ((inc >> 3) & 1) * deg;
	while (len--) {
		u64 mod = *p++;
		if ((s64)mod > 0)
			hexl::EltwiseFMAMod(rop, op1, *op2 % mod, op3, deg, mod, 1);
		rop += inc0, op1 += inc1, op2 += inc2, op3 += inc3;
	}
}

void
ring_poly_rndmulcadd(u64 *rop, f64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc)
{
	u64 *p = mods.mods;
	i32 deg = mods.deg, len = mods.len;
	i32 dec0 = ((~inc >> 0) & 1) * deg;
	i32 dec1 = ((~inc >> 1) & 1) * deg;
	i32 inc2 = (( inc >> 2) & 1);
	i32 dec3 = ((~inc >> 3) & 1) * deg;
	while (len--) {
		u64 mod = *p++;
		u64 c = *op2 % mod;
		i32 i = deg / 8;
		while (i--) {
			if ((s64)mod > 0) {
				rop[0] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[0] + 0.5), c, mod), op3[0], mod);
				rop[1] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[1] + 0.5), c, mod), op3[1], mod);
				rop[2] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[2] + 0.5), c, mod), op3[2], mod);
				rop[3] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[3] + 0.5), c, mod), op3[3], mod);
				rop[4] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[4] + 0.5), c, mod), op3[4], mod);
				rop[5] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[5] + 0.5), c, mod), op3[5], mod);
				rop[6] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[6] + 0.5), c, mod), op3[6], mod);
				rop[7] = hexl::AddUIntMod(hexl::MultiplyMod((u64)(op1[7] + 0.5), c, mod), op3[7], mod);
			}
			rop += 8, op1 += 8, op3 += 8;
		}
		rop -= dec0, op1 -= dec1, op2 += inc2, op3 -= dec3;
	}
}

void
ring_poly_divcadd(f64 *rop, u64 const *op1, u64 const *op2, f64 const *op3, i32 deg, i32 len, u8 inc)
{
	i32 dec0 = ((~inc >> 0) & 1) * deg;
	i32 dec1 = ((~inc >> 1) & 1) * deg;
	i32 inc2 = (( inc >> 2) & 1);
	i32 dec3 = ((~inc >> 3) & 1) * deg;
	while (len--) {
		f64 c = (f64)*op2;
		i32 i = deg / 8;
		while (i--) {
			rop[0] = (f64)op1[0] / c + op3[0];
			rop[1] = (f64)op1[1] / c + op3[1];
			rop[2] = (f64)op1[2] / c + op3[2];
			rop[3] = (f64)op1[3] / c + op3[3];
			rop[4] = (f64)op1[4] / c + op3[4];
			rop[5] = (f64)op1[5] / c + op3[5];
			rop[6] = (f64)op1[6] / c + op3[6];
			rop[7] = (f64)op1[7] / c + op3[7];
			rop += 8, op1 += 8, op3 += 8;
		}
		rop -= dec0, op1 -= dec1, op2 += inc2, op3 -= dec3;
	}
}

void
ring_poly_nttbwd(u64 *op, ring_mods mods)
{
	u64 *p = mods.mods;
	void **pp = mods.precomp;
	i32 deg = mods.deg, len = mods.len;
	while (len--) {
		auto *ntt = (hexl::NTT *)*pp++;
		u64 mod = *p++;
		if ((s64)mod > 0)
			// NOTE: in-place NTT internally
			ntt->ComputeInverse(op, op, 1, 1);
		op += deg;
	}
}

void
ring_poly_nttfwd(u64 *op, ring_mods mods)
{
	u64 *p = mods.mods;
	void **pp = mods.precomp;
	i32 deg = mods.deg, len = mods.len;
	while (len--) {
		auto *ntt = (hexl::NTT *)*pp++;
		u64 mod = *p++;
		if ((s64)mod > 0)
			// NOTE: in-place NTT internally
			ntt->ComputeForward(op, op, 1, 1);
		op += deg;
	}
}
