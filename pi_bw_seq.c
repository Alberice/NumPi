/* Pi via Borwein quadrático, sequencial */
#include <stdio.h>
#include <gmp.h>

#define N 1000000000

int main() {
	mpf_t a, b, c, p, a1, b1, p1;
	mpf_init_set_ui(a, 2);
	mpf_sqrt(a, a);
	mpf_init_set_ui(b, 0);
	mpf_init(c);
	mpf_init_set_ui(p, 2);
	mpf_add(p, p, a);
	mpf_init(a1);
	mpf_init(b1);
	mpf_init(p1);
	for (int i=1; i<N; i++) {
		mpf_sqrt(c, a);
		mpf_ui_div(a1, 1, c);
		mpf_add(a1, a1, c);
		mpf_div_ui(a1, a1, 2);

		mpf_add_ui(b1, b, 1);
		mpf_mul(b1, b1, c);
		mpf_add(c, a, b);
		mpf_div(b1, b1, c);

		mpf_add_ui(p1, a1, 1);
		mpf_mul(p1, p1, p);
		mpf_mul(p1, p1, b1);
		mpf_add_ui(c, b1, 1);
		mpf_div(p1, p1, c);

		mpf_swap(a, a1);
		mpf_swap(b, b1);
		mpf_swap(p, p1);
	}
	gmp_printf("%Ff\n", p);
	return 0;
}
