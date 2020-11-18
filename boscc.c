#include <stdio.h>
#include <arm_sve.h>
#include <stdlib.h>

double a[4096];
double b[4096];
double c[4096];
double d[4096];
int64_t bools[4096];

void init_arrays (int n) {
    for (int i = 0; i < 4096; i++) {
        a[i] = (double) i;
        b[i] = (double) i / 10;
        c[i] = (double) i / 5;
        d[i] = (double) i / 2;
        bools[i] = i % n;
    }
}

void reg_kern (void) {

    for (int i = 0; i < 4096; i++) {

        a[i:i+8]*3
        b[i:i+8]

        if ( bools[i] ) {
            a[i] = b[i] + c[i];
        }
        else {
            a[i] = b[i] + d[i];
        }
    }

    for (int i = 0; i < 100; i++) {
        printf ("%f ", a[i]);
    }
    printf ("\n");

}

void boscc_kern (void) {
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

   int boscc_count = 0;

    for (; i < n;) {
        //a[i] = b[i] + c[i] + d[i];

        svint64_t bools_v = svld1_s64 (pg, bools+i);

        // else cond
        svbool_t cond1 = svcmpeq_n_s64 (pg, bools_v, 0);

        // if cond
        svbool_t cond2 = svcmpne_n_s64 (pg, bools_v, 0);

        // if side all true
        int boscc_if = svcntp_b64(pg, cond2) == svcntd();
        if (boscc_if) {
            svfloat64_t vb = svld1_f64 (pg, b+i);
            svfloat64_t vc = svld1_f64 (pg, c+i);
            svfloat64_t res = svadd_f64_z (pg, vb, vc);

            svst1_f64 (pg, a+i, res);

            boscc_count++;
        } else { // execute if cnvrt
            // if converted
            // if side
            svfloat64_t vb = svld1_f64 (pg, b+i);
            svfloat64_t vc = svld1_f64 (pg, c+i);

            // else side
            svfloat64_t vd = svld1_f64 (pg, d+i);

            svfloat64_t res1 = svadd_f64_m (cond2, vb, vc);
            svfloat64_t res2 = svadd_f64_m (cond1, vb, vd);
            svfloat64_t res = svsel (cond2, res1, res2);

            svst1_f64 (pg, a+i, res);
        }

        i += svcntd();
        pg = svwhilelt_b64 (i, n);
    }

    for (int i = 0; i < 100; i++) {
        printf ("%f ", a[i]);
    }
    printf ("\n");

    printf ("boscc count: %d\n", boscc_count);
}

void boscc_kern_unrolled (void) {
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

   int boscc_count = 0;

   svbool_t all_true = svptrue_b64();

    for (; i < n;) {
        //a[i] = b[i] + c[i] + d[i];

        svint64_t bools_v = svld1_s64 (pg, bools+i);

        // else cond
        svbool_t cond1 = svcmpeq_n_s64 (pg, bools_v, 0);

        // if cond
        svbool_t cond2 = svcmpne_n_s64 (pg, bools_v, 0);

        // if side all true
        int boscc_if = svcntp_b64(pg, cond2) == svcntd();
        if (boscc_if) {
            svfloat64_t vb = svld1_f64 (pg, b+i);
            svfloat64_t vc = svld1_f64 (pg, c+i);
            svfloat64_t res = svadd_f64_z (pg, vb, vc);

            svst1_f64 (pg, a+i, res);

            boscc_count++;
        } else { // execute if cnvrt
            // if converted
            // if side
            svfloat64_t vb = svld1_f64 (pg, b+i);
            svfloat64_t vc = svld1_f64 (pg, c+i);

            // else side
            svfloat64_t vd = svld1_f64 (pg, d+i);

            svfloat64_t res1 = svadd_f64_m (cond2, vb, vc);
            svfloat64_t res2 = svadd_f64_m (cond1, vb, vd);
            svfloat64_t res = svsel (cond2, res1, res2);

            svst1_f64 (pg, a+i, res);
        }

        i += svcntd();
        pg = svwhilelt_b64 (i, n);
    }

    for (int i = 0; i < 100; i++) {
        printf ("%f ", a[i]);
    }
    printf ("\n");

    printf ("boscc count: %d\n", boscc_count);
}


void test_kern (void) {
    int64_t i = 0;
    int64_t n = 4096;
    svbool_t pg = svwhilelt_b64 (i, n);


    /*

    for ...
        if bools[i]
            ai = bi + ci;

    */

   int boscc_count = 0;

    for (; i < n;) {
        // sum += ai

        // index vector
        svint64_t idx = svindex_s64 (i, 1);

        svint64_t bools_v = svld1_s64 (pg, bools+i);
        svint64_t bools_v2 = svld1_s64 (pg, bools+i);

        // if cond
        svbool_t cond = svcmpne_n_s64 (pg, bools_v, 0);
        svbool_t cond2 = svcmpne_n_s64 (pg, bools_v2, 0);

        svint64_t bools_v2comp = svcompact_s64 (cond2, bools_v2);

        svint64_t collect = svsplice_s64 (cond, bools_v, bools_v2comp);

        svbool_t cond_collect = svcmpne_n_s64 (pg, collect, 0);

        int active = svcntp_b64 (cond_collect);

        svfloat64_t a1,a2;

        if (active) {
            a1 = svld1_f64 (pg, a+i);
            a2 = svld1_f64 (pg, a+i+svlen_f64());

            a2_comp = svcompact_f64 (cond2, bools_v2);

            a12 = svsplice_s64 (cond, a1, a2_comp);
        }

        i += svcntd();
        pg = svwhilelt_b64 (i, n);
    }

    for (int i = 0; i < 100; i++) {
        printf ("%f ", a[i]);
    }
    printf ("\n");

    printf ("boscc count: %d\n", boscc_count);
}

int main(int argc, char **argv) {

    int n = atoi (argv[1]);

    init_arrays (n);

    reg_kern();
    boscc_kern();

    return 0;
}
