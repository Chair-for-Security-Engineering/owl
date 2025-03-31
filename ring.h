#include <stdint.h>

typedef  int8_t  s8;
typedef uint8_t  u8;
typedef  int32_t s32;
typedef uint32_t u32;
typedef  int64_t s64;
typedef uint64_t u64;
typedef double   f64;

__extension__ typedef __int128 s128;
__extension__ typedef unsigned __int128 u128;

typedef struct ring_mods {
	void **precomp;
	u64 *mods;
	u64 deg;
	u64 len;
} ring_mods;

#define RING_MOD_NOP(b) ((u64)(b) << 63)
#define RING_MOD_USE    0x7fffffffffffffffULL

#define RING_INC_110    0x03
#define RING_INC_111    0x07
#define RING_INC_0110   0x06
#define RING_INC_1011   0x0d
#define RING_INC_1100   0x03
#define RING_INC_1101   0x0b
#define RING_INC_1110   0x07
#define RING_INC_1111   0x0f

ring_mods ring_mods_init(u64 *mods, u64 deg, u64 len);
void      ring_mods_free(ring_mods mods);
u64       ring_coef_inv(u64 op, u64 mod);
u64       ring_coef_mul(u64 op1, u64 op2, u64 mod);

void ring_poly_add(u64 *rop, u64 const *op1, u64 const *op2, ring_mods mods, u8 inc);
void ring_poly_mul(u64 *rop, u64 const *op1, u64 const *op2, ring_mods mods, u8 inc);
void ring_poly_muladd(u64 *rop, u64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc);
void ring_poly_mulcadd(u64 *rop, u64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc);
void ring_poly_rndmulcadd(u64 *rop, f64 const *op1, u64 const *op2, u64 const *op3, ring_mods mods, u8 inc);
void ring_poly_divcadd(f64 *rop, u64 const *op1, u64 const *op2, f64 const *op3, u64 deg, u64 len, u8 inc);
void ring_poly_nttbwd(u64 *op, ring_mods mod);
void ring_poly_nttfwd(u64 *op, ring_mods mod);
