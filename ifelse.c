#include <stdio.h>
#include <arm_sve.h>

double a[4096];
double b[4096];
double c[4096];
double d[4096];
int64_t bools[4096];

int main(void) {

    int64_t i = 0;
    int64_t n = 4096;
    svbool_t pg = svwhilelt_b64 (i, n);


    /*

    for ...
        if bools[i]
            ai = bi + ci;
        else
            ai = bi + di;

    */

    for (; i < n;) {
        //a[i] = b[i] + c[i] + d[i];

        svint64_t b = svld1_s64 (pg, bools+i);
        svbool_t cond1 = svcmpeq_n_s64 (pg, b, 0);
        svbool_t cond2 = svcmpne_n_s64 (pg, b, 0);

        svfloat64_t vb = svld1_f64 (pg, b+i);
        svfloat64_t vc = svld1_f64 (pg, c+i);
        svfloat64_t vd = svld1_f64 (pg, d+i);

        svfloat64_t res = svadd_f64_z (cond1, vb, vc);
        res = svadd_f64_m (cond2, vb, vd);

        svst1_f64 (pg, a+i, res);

        i += svcntd();
        pg = svwhilelt_b64 (i, n);
        printf ("%d\n", i);
    }

    printf ("%f\n", a[0]);

    return 0;
}
