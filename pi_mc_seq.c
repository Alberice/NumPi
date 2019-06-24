/* Pi via Monte Carlo, sequencial */
#include <stdio.h>
#include <gmp.h>

#define N 1000000000

int main() {
	gmp_randstate_t randstate;
	gmp_randinit_default(randstate);
	mp_bitcnt_t nbits = 100;
	mpf_t x, y, p;
	mpf_init(x);
	mpf_init(y);
	mpf_init(p);
	for (int i=0; i<N; i++) {
		mpf_urandomb(x, randstate, nbits);
		mpf_urandomb(y, randstate, nbits);
		mpf_mul(x, x, x); // x=x^2
		mpf_mul(y, y, y); // y=y^2
		mpf_add(x, x, y); // x=norma^2
		if (mpf_cmp_ui(x, 1) <= 0) {
			mpf_add_ui(p, p, 1);
		}
	}
	mpf_div_ui(p, p, N/4);
	gmp_printf("%Ff\n", p);
	return 0;
}
