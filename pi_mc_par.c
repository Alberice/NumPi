/* Pi via Monte Carlo, paralelo */
#include <stdio.h>
#include <gmp.h>
#include <pthread.h>

#define N 1000000000
#define NJOBS 10

void* work(void *out) {
	gmp_randstate_t randstate;
	gmp_randinit_default(randstate);
	mp_bitcnt_t nbits = 100;
	mpf_t x, y, p;
	mpf_init(x);
	mpf_init(y);
	mpf_init(p);
	for (int i=0; i<N/NJOBS; i++) {
		mpf_urandomb(x, randstate, nbits);
		mpf_urandomb(y, randstate, nbits);
		mpf_mul(x, x, x); // x=x^2
		mpf_mul(y, y, y); // y=y^2
		mpf_add(x, x, y); // x=norma^2
		if (mpf_cmp_ui(x, 1) <= 0) {
			mpf_add_ui(p, p, 1);
		}
	}
	mpf_set(*(mpf_t*)out, p);
	return 0;
}

int main() {
	pthread_t threads[NJOBS];
	mpf_t ret[NJOBS];
	for (int i=0; i<NJOBS; i++) {
		mpf_init(ret[i]);
		pthread_create(&threads[i], NULL, work, &ret[i]);
	}
	for (int i=0; i<NJOBS; i++) {
		pthread_join(threads[i], NULL);
	}
	mpf_t p;
	mpf_init(p);
	for (int i=0; i<NJOBS; i++) {
		mpf_add(p, p, ret[i]);
	}
	mpf_div_ui(p, p, N/4);
	gmp_printf("%Ff\n", p);
	return 0;
}
