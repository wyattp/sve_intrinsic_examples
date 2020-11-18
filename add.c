#include <stdio.h>
#include <arm_sve.h>

double a[4096];
double b[4096];
double c[4096];

volatile int * d;

int main(void) {

    int64_t i = 0;
    int64_t n = 4096;
    svbool_t pg = svwhilelt_b64 (i, n);

    for (; i < n;) {
        //a[i] = b[i] + c[i] + d[i];
        svfloat64_t vb = svld1_f64 (pg, b);
        svfloat64_t vc = svld1_f64 (pg, c);
        svfloat64_t res = svadd_f64_z (pg, vb, vc);
        svst1_f64 (pg, a, res);
        i += svcntd();
        pg = svwhilelt_b64 (i, n);
        printf ("%d\n", i);
    }

    printf ("%f\n", a[0]);

    return 0;
}