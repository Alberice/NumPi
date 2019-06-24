/* Black-Scholes via Monte Carlo, paralelo */
#include <stdio.h>
#include <pthread.h>
#include <gmp.h>
#include <mpfr.h>

#define NJOBS 10

long M;
mpf_t x[5];
mpfr_t S, E, r, s, T, a, b, c, t, s1, s2;

typedef struct {
	mpfr_t s1, s2;
} args_t;

void* work(void *args) {
	gmp_randstate_t randstate;
	gmp_randinit_default(randstate);
	mpfr_t a, b, c, t, s1, s2;
	mpfr_init(a);
	mpfr_init(b);
	mpfr_init(c);
	mpfr_init(t);
	mpfr_init_set_ui(s1, 0, MPFR_RNDN);
	mpfr_init_set_ui(s2, 0, MPFR_RNDN);
	for (int i=0; i<M/NJOBS; i++) {
		// a = (r - 1/2 sigma^2) * T
		mpfr_mul(a, s, s, MPFR_RNDN); // sigma^2
		mpfr_div_ui(a, a, 2, MPFR_RNDN); // *1/2
		mpfr_sub(a, r, a, MPFR_RNDN); // r - a
		mpfr_mul(a, a, T, MPFR_RNDN); // *T
		// b = sigma * sqrt(T) * randomNumber
		mpfr_urandomb(c, randstate);
		mpfr_sqrt(b, T, MPFR_RNDN);
		mpfr_mul(b, b, c, MPFR_RNDN);
		mpfr_mul(b, b, s, MPFR_RNDN);
		// t
		mpfr_add(t, a, b, MPFR_RNDN);
		mpfr_exp(t, t, MPFR_RNDN);
		mpfr_mul(t, S, t, MPFR_RNDN);
		// trials
		mpfr_neg(a, r, MPFR_RNDN);
		mpfr_mul(a, a, T, MPFR_RNDN);
		mpfr_exp(a, a, MPFR_RNDN);
		mpfr_sub(t, t, E, MPFR_RNDN);
		if (mpfr_cmp_ui(t, 0) > 0) {
			mpfr_mul(t, t, a, MPFR_RNDN);
			mpfr_add(s1, s1, t, MPFR_RNDN);
			mpfr_mul(t, t, t, MPFR_RNDN);
			mpfr_add(s2, s2, t, MPFR_RNDN);
		}
	}
	args_t *ret = (args_t *) args;
	mpfr_swap(ret->s1, s1);
	mpfr_swap(ret->s2, s2);
	return 0;
}

int main() {
	for (int i=0; i<5; i++) {
		mpf_init(x[i]);
		gmp_scanf(" %Ff", x[i]);
	}
	scanf(" %ld", &M);
	mpfr_init_set_f(S, x[0], MPFR_RNDN);
	mpfr_init_set_f(E, x[1], MPFR_RNDN);
	mpfr_init_set_f(r, x[2], MPFR_RNDN);
	mpfr_init_set_f(s, x[3], MPFR_RNDN);
	mpfr_init_set_f(T, x[4], MPFR_RNDN);
	mpfr_init(a);
	mpfr_init(b);
	mpfr_init(c);
	mpfr_init(t);
	mpfr_init_set_ui(s1, 0, MPFR_RNDN);
	mpfr_init_set_ui(s2, 0, MPFR_RNDN);

	pthread_t threads[NJOBS];
	args_t ret[NJOBS];
	for (int i=0; i<NJOBS; i++) {
		mpfr_init(ret[i].s1);
		mpfr_init(ret[i].s2);
		pthread_create(&threads[i], NULL, work, &ret[i]);
	}
	for (int i=0; i<NJOBS; i++) {
		pthread_join(threads[i], NULL);
	}
	for (int i=0; i<NJOBS; i++) {
		mpfr_add(s1, s1, ret[i].s1, MPFR_RNDN);
		mpfr_add(s2, s2, ret[i].s2, MPFR_RNDN);
	}
	mpfr_div_ui(s1, s1, M, MPFR_RNDN); // media
	mpfr_div_ui(s2, s2, M, MPFR_RNDN);
	mpfr_mul(a, s1, s1, MPFR_RNDN);
	mpfr_sub(b, s2, a, MPFR_RNDN);
	mpfr_sqrt(a, b, MPFR_RNDN); // stddev
	mpfr_sqrt_ui(b, M, MPFR_RNDN); // sqrt(M)
	mpfr_set_str(c, "1.96", 10, MPFR_RNDN);
	mpfr_mul(c, c, a, MPFR_RNDN);
	mpfr_div(c, c, b, MPFR_RNDN); // confwidth
	mpfr_sub(a, s1, c, MPFR_RNDN);
	mpfr_add(b, s1, c, MPFR_RNDN);

	mpfr_printf("S       %.0Rf\n", S);
	mpfr_printf("E       %.0Rf\n", E);
	mpfr_printf("r       %.0Rf\n", r);
	mpfr_printf("sigma   %.0Rf\n", s);
	mpfr_printf("T       %.0Rf\n", T);
	printf("M       %d\n", M);
	mpfr_printf("Confidence interval: (%Rf, %Rf)\n", a, b);
	return 0;
}
