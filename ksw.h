#ifndef KSW_H
#define KSW_H

#include "tiifhe.h"

void keyswitch_double(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0);
void keyswitch_single(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0);
void keyswitch_single_folded(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t idx0);
void keyswitch_single_switch(tiifhe_Ciphertext *ct, size_t dsw, int domain, size_t mu, size_t kappa);

void keyswitch_init(void);
void keyswitch_precomp(const tiifhe_KeySecret *sk);
void keyswitch_precomp_decomp(void);
void keyswitch_precomp_fold(void);

#endif /* KSW_H */
