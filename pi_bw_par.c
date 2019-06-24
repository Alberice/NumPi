/* Pi via Borwein quadrático, paralelo */
#include <stdio.h>
#include <gmp.h>
#include <pthread.h>

#define N 15
#define Q 10

pthread_mutex_t lock_a = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_b = PTHREAD_MUTEX_INITIALIZER;

mpf_t a[Q];
mpf_t b[Q];
int ai, an, bi, bn;

void* work_a(void *arg) {
	mpf_t x, a1;
	mpf_init(x);
	mpf_init(a1);
	for (int i=1; i<N; i++) {
		while (an < Q) {
			if (!pthread_mutex_trylock(&lock_a) && an < Q) {
				int m = (ai + an - 1) % Q;
				int n = (ai + an) % Q;
				mpf_sqrt(x, a[m]);
				mpf_ui_div(a1, 1, x);
				mpf_add(a1, a1, x);
				mpf_div_ui(a1, a1, 2);
				mpf_swap(a[n], a1);
				an++;
				printf("a: %d\n", an);
				pthread_mutex_unlock(&lock_a);
			}
			while (an == Q);
		}
	}
	return 0;
}

void* work_b(void *arg) {
	mpf_t x, b1;
	mpf_init(x);
	mpf_init(b1);
	for (int i=1; i<N; i++) {
		while (bn < Q) {
			if (!pthread_mutex_trylock(&lock_b) && bn < Q) {
				int m = (bi + bn - 1) % Q;
				int n = (bi + bn) % Q;
				mpf_add_ui(b1, b[m], 1);
				mpf_mul(b1, b1, x);
				while (pthread_mutex_trylock(&lock_a)) {
					int ma = (ai + an - 1) % Q;
					mpf_add(x, a[ma], b[m]);
					pthread_mutex_unlock(&lock_a);
					break;
				}
				mpf_div(b1, b1, x);
				mpf_swap(b[n], b1);
				bn++;
				printf("b: %d %d\n", an, bn);
				pthread_mutex_unlock(&lock_b);
			}
			while (bn == Q);
		}
	}
	return 0;
}

void* work_p(void *arg) {
	mpf_t p0, p1, x;
	mpf_init_set_ui(p0, 2);
	mpf_sqrt(p0, p0);
	mpf_add_ui(p0, p0, 2);
	mpf_init(x);
	for (int i=1; i<N; i++) {
		while (bn > 1) {
			if (!pthread_mutex_trylock(&lock_b) && bn > 1) {
				while (an > 1) {
					if (!pthread_mutex_trylock(&lock_a) && an > 1) {
						int ma = (ai + an - 1) % Q;
						int mb = (bi + bn - 1) % Q;
						mpf_add_ui(p1, a[ma], 1);
						mpf_mul(p1, p1, p0);
						mpf_mul(p1, p1, b[mb]);
						mpf_add_ui(x, b[mb], 1);
						mpf_div(p1, p1, x);
						mpf_swap(p0, p1);
						ai = (ai + 1) % Q;
						an--;
						bi = (bi + 1) % Q;
						bn--;
						printf("p: %d %d %d\n", an, bn, i);
						pthread_mutex_unlock(&lock_a);
					}
				}
				pthread_mutex_unlock(&lock_b);
			}
		}
	}
	gmp_printf("%Ff\n", p0);
	printf("%d %d\n", an, bn);
}

int main() {
	ai = bi = 0;
	an = bn = 1;
	mpf_init_set_ui(a[0], 2);
	mpf_sqrt(a[0], a[0]);
	mpf_init_set_ui(b[0], 0);
	for (int i=1; i<Q; i++) {
		mpf_init(a[i]);
		mpf_init(b[i]);
	}
	pthread_t ta, tb, tp;
	pthread_create(&ta, NULL, work_a, NULL);
	pthread_create(&tb, NULL, work_b, NULL);
	pthread_create(&tp, NULL, work_p, NULL);
	pthread_join(ta, NULL);
	pthread_join(tb, NULL);
	pthread_join(tp, NULL);
	return 0;
}
