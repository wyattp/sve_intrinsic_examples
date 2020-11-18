#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <arm_sve.h>
#include <stdlib.h>

struct struct_prm {
    int64_t     Ntypes;
    double *    Cn1;
    double *    Cn2;
    int64_t *   Cno;
};

int64_t i;
int64_t iaci;
double * atype;
double sigma_iw6;
double epsilon_iw;
double sigmaw3;

#define EPSILONW 3.3
#define PI 3.1415926
#define RHOW 7.54

struct struct_prm * prm;

void unknown_func (struct struct_prm *p) {
    prm = malloc(sizeof(struct struct_prm));

    prm->Ntypes = 4096;

    prm->Cno = malloc(sizeof(int64_t)*prm->Ntypes*prm->Ntypes);
    prm->Cn1 = malloc(sizeof(double)*prm->Ntypes*prm->Ntypes);
    prm->Cn2 = malloc(sizeof(double)*prm->Ntypes*prm->Ntypes);

    atype = malloc(sizeof(double) * prm->Ntypes*prm->Ntypes);

    for (int i = 0; i < prm->Ntypes; i++) {
        prm->Cno[prm->Ntypes*i+i] = i+1;
        //prm->Cn1[i] = (double) 1./(double)i;
        //prm->Cn2[i] = (double) 1./(double)i;
        prm->Cn1[i+1] = i%3;
        prm->Cn2[i+1] = 1;
    }
}



void ifconversion_kern() {

    svbool_t loop_pred = svwhilelt_b64_s64 (i, prm->Ntypes);
    svbool_t all_true = svptrue_b64 ();
    svint64_t i_v = svindex_s64(0,1);

    while (svptest_any (all_true, loop_pred)) {
        // iaci = prm->Cno[prm->Ntypes * i + i] - 1; 
        svint64_t offsets = svmul_n_s64_z (loop_pred, i_v, prm->Ntypes);
        offsets = svadd_s64_z (loop_pred, offsets, i_v);
        svint64_t iaci_v = svld1_gather_s64offset_s64 (loop_pred, prm->Cno, offsets);
        iaci_v = svadd_n_s64_z (loop_pred, iaci_v, -1);

        // build predicate for
        // if (prm->Cn1[iaci] == 0.0 || prm->Cn2[iaci] == 0.0)
        svfloat64_t cn1ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn1, iaci_v);
        svfloat64_t cn2ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn2, iaci_v);
        svbool_t cond1_pred = svcmpeq_n_f64 (loop_pred, cn1ld, 0);
        svbool_t cond2_pred = svcmpeq_n_f64 (loop_pred, cn2ld, 0);
        svbool_t cond_pred = svorr_b_z (loop_pred, cond1_pred, cond2_pred);

        // if
        svst1_f64 (cond_pred, atype+i, svdup_f64(0) );

        svbool_t not_cond_pred = svnot_b_z (loop_pred, cond_pred);

        // else
        svfloat64_t sigmaw3_v = svdup_n_f64(sigmaw3);
        svfloat64_t div = svdiv_f64_z (not_cond_pred, cn1ld, cn2ld);
        svfloat64_t sigma_iw6_v = svmul_f64_z (not_cond_pred, sigmaw3_v, svsqrt_f64_z (not_cond_pred, div));
        svfloat64_t epsilon_iw_v = svmul_f64_z(not_cond_pred, svmul_f64_z (not_cond_pred, svdup_n_f64 (0.5),svsqrt_f64_z(not_cond_pred, svdiv_f64_z (not_cond_pred, svdup_n_f64(EPSILONW),cn1ld))),cn2ld);
        svfloat64_t st = svdiv_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svdup_n_f64(-16),svdup_n_f64(PI)),svdup_n_f64(RHOW)),epsilon_iw_v),sigma_iw6_v) , svdup_n_f64(3.0) );
        svst1_f64 (not_cond_pred, atype+i, st);

        i_v = svadd_n_s64_z (loop_pred, i_v, svcntd());
        i += svcntd();
        loop_pred = svwhilelt_b64_s64 (i, prm->Ntypes);
    }

    printf ("%f\n", atype[0]);
}

void boscc_unrolled_kern() {

    i = 0;
    svbool_t loop_pred1 = svwhilelt_b64_s64 (i, prm->Ntypes);
    svbool_t all_true = svptrue_b64 ();
    svint64_t i_v = svindex_s64(0,1);

    i += svcntd();
    svbool_t loop_pred2 = svwhilelt_b64_s64 (i, prm->Ntypes);
    svint64_t i_v2 = svindex_s64(i,1);

    while (svptest_any (all_true, loop_pred1)) {
        // iaci = prm->Cno[prm->Ntypes * i + i] - 1; 
        svint64_t offsets = svmul_n_s64_z (loop_pred1, i_v, prm->Ntypes);
        offsets = svadd_s64_z (loop_pred1, offsets, i_v);
        svint64_t iaci_v = svld1_gather_s64offset_s64 (loop_pred1, prm->Cno, offsets);
        iaci_v = svadd_n_s64_z (loop_pred1, iaci_v, -1);

        svint64_t offsets2 = svmul_n_s64_z (loop_pred2, i_v2, prm->Ntypes);
        offsets2 = svadd_s64_z (loop_pred2, offsets2, i_v2);
        svint64_t iaci_v_2 = svld1_gather_s64offset_s64 (loop_pred2, prm->Cno, offsets2);
        iaci_v_2 = svadd_n_s64_z (loop_pred2, iaci_v_2, -1);

        // build predicate for
        // if (prm->Cn1[iaci] == 0.0 || prm->Cn2[iaci] == 0.0)
        svfloat64_t cn1ld = svld1_gather_s64offset_f64 (loop_pred1, prm->Cn1, iaci_v);
        svbool_t cond1_pred = svcmpeq_n_f64 (loop_pred1, cn1ld, 0);

        svfloat64_t cn2ld = svld1_gather_s64offset_f64 (loop_pred1, prm->Cn2, iaci_v);
        svbool_t cond2_pred = svcmpeq_n_f64 (loop_pred1, cn2ld, 0);
        svbool_t cond_pred = svorr_b_z (loop_pred1, cond1_pred, cond2_pred);
        svbool_t not_cond_pred = svnot_b_z (loop_pred1, cond_pred);

        svfloat64_t cn1ld2 = svld1_gather_s64offset_f64 (loop_pred2, prm->Cn1, iaci_v);
        svbool_t cond1_pred2 = svcmpeq_n_f64 (loop_pred2, cn1ld2, 0);
        svfloat64_t cn2ld2 = svld1_gather_s64offset_f64 (loop_pred2, prm->Cn2, iaci_v);
        svbool_t cond2_pred2 = svcmpeq_n_f64 (loop_pred2, cn2ld2, 0);
        svbool_t cond_pred2 = svorr_b_z (loop_pred2, cond1_pred2, cond2_pred2);
        svbool_t not_cond_pred2 = svnot_b_z (loop_pred2, cond_pred2);

        // make uniform vector
        svfloat64_t splice = svsplice_f64 (cond_pred, cn1ld, svcompact_f64 (cond_pred2, cn1ld2));
        svfloat64_t splice2 = svsplice_f64 (not_cond_pred2, cn1ld2, svcompact_f64 (not_cond_pred, cn1ld));

        // merged vector
        svbool_t ucmp = svcmpeq_n_f64 (all_true, splice, 0);
        svbool_t not_ucmp = svnot_b_z (all_true, ucmp);

        // remainder vector
        svbool_t ucmp2 = svcmpeq_n_f64 (all_true, splice2, 0);
        svbool_t not_ucmp2 = svnot_b_z (all_true, ucmp2);

        if (!svptest_any(all_true, not_ucmp)) {
            svint64_t idxmerge = svsplice_s64 (cond_pred, i_v, svcompact(cond_pred2, i_v2));
            svst1_scatter_s64index_f64 (ucmp, atype, idxmerge, svdup_n_f64(0.0));
        }
        else {
            // IF CONVERTED CODE

            // if
            svint64_t idxmerge = svsplice_s64 (cond_pred, i_v, svcompact(cond_pred2, i_v2));
            svst1_f64 (ucmp, atype, svdup_f64(0) );

            svfloat64_t cn1merge = svsplice_f64 (cond_pred, cn1ld, svcompact_f64 (cond_pred,cn1ld2));
            svfloat64_t cn2merge = svsplice_f64 (cond_pred, cn2ld, svcompact_f64 (cond_pred, cn2ld2));

            // else
            svfloat64_t sigmaw3_v = svdup_n_f64(sigmaw3);
            svfloat64_t div = svdiv_f64_z (not_ucmp, cn1merge, cn2merge);
            svfloat64_t sigma_iw6_v = svmul_f64_z (not_ucmp, sigmaw3_v, svsqrt_f64_z (not_ucmp, div));
            svfloat64_t epsilon_iw_v = svmul_f64_z(not_ucmp, svmul_f64_z (not_ucmp, svdup_n_f64 (0.5),svsqrt_f64_z(not_ucmp, svdiv_f64_z (not_ucmp, svdup_n_f64(EPSILONW),cn1ld))),cn2ld);
            svfloat64_t st = svdiv_f64_z(not_ucmp, svmul_f64_z(not_ucmp, svmul_f64_z(not_ucmp, svmul_f64_z(not_ucmp, svmul_f64_z(not_ucmp, svdup_n_f64(-16),svdup_n_f64(PI)),svdup_n_f64(RHOW)),epsilon_iw_v),sigma_iw6_v) , svdup_n_f64(3.0) );
            svst1_scatter_s64index_f64(not_ucmp, atype, idxmerge, st);
        }

        if (!svptest_any(all_true, not_ucmp2)) {
            //svst1_f64 (cond_pred, atype+i, svdup_n_f64(0));
            svint64_t idxremainder = svsplice_s64 (not_cond_pred, i_v, svcompact(not_cond_pred2, i_v2));
            svst1_scatter_s64index_f64 (ucmp2, atype, idxremainder, svdup_n_f64(0.0));
        }
        else {
            // IF CONVERTED CODE

            // if
            svint64_t idxremainder = svsplice_s64 (not_cond_pred, i_v, svcompact(not_cond_pred2, i_v2));
            svst1_scatter_s64index_f64(ucmp2, atype, idxremainder, svdup_f64(0.0));

            svfloat64_t cn1merge = svsplice_f64 (not_cond_pred, cn1ld, svcompact_f64 (not_cond_pred, cn1ld2));
            svfloat64_t cn2merge = svsplice_f64 (not_cond_pred, cn2ld, svcompact_f64 (not_cond_pred, cn2ld2));

            // else
            svfloat64_t sigmaw3_v = svdup_n_f64(sigmaw3);
            svfloat64_t div = svdiv_f64_z (not_ucmp2, cn1merge, cn2merge);
            svfloat64_t sigma_iw6_v = svmul_f64_z (not_ucmp2, sigmaw3_v, svsqrt_f64_z (not_ucmp2, div));
            svfloat64_t epsilon_iw_v = svmul_f64_z(not_ucmp2, svmul_f64_z (not_ucmp2, svdup_n_f64 (0.5),svsqrt_f64_z(not_ucmp2, svdiv_f64_z (not_ucmp2, svdup_n_f64(EPSILONW),cn1ld))),cn2ld);
            svfloat64_t st = svdiv_f64_z(not_ucmp2, svmul_f64_z(not_ucmp2, svmul_f64_z(not_ucmp2, svmul_f64_z(not_ucmp2, svmul_f64_z(not_ucmp2, svdup_n_f64(-16),svdup_n_f64(PI)),svdup_n_f64(RHOW)),epsilon_iw_v),sigma_iw6_v) , svdup_n_f64(3.0) );
            svst1_scatter_s64index_f64(not_ucmp2, atype, idxremainder, st);
        }

        i += svcntd();
        i_v = svadd_n_s64_z (loop_pred1, i_v, svcntd());
        loop_pred1 = svwhilelt_b64_s64 (i, prm->Ntypes);

        i += svcntd();
        i_v2 = svadd_n_s64_z (loop_pred2, i_v2, svcntd());
        loop_pred2 = svwhilelt_b64_s64 (i, prm->Ntypes);
    }

    printf ("%f\n", atype[0]);
}

void boscc_kern() {
    svbool_t loop_pred = svwhilelt_b64_s64 (i, prm->Ntypes);
    svbool_t all_true = svptrue_b64 ();
    svint64_t i_v = svindex_s64(0,1);

    while (svptest_any (all_true, loop_pred)) {
        // iaci = prm->Cno[prm->Ntypes * i + i] - 1; 
        svint64_t offsets = svmul_n_s64_z (loop_pred, i_v, prm->Ntypes);
        offsets = svadd_s64_z (loop_pred, offsets, i_v);
        svint64_t iaci_v = svld1_gather_s64offset_s64 (loop_pred, prm->Cno, offsets);
        iaci_v = svadd_n_s64_z (loop_pred, iaci_v, -1);

        // build predicate for
        // if (prm->Cn1[iaci] == 0.0 || prm->Cn2[iaci] == 0.0)
        svfloat64_t cn1ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn1, iaci_v);
        svfloat64_t cn2ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn2, iaci_v);
        svbool_t cond1_pred = svcmpeq_n_f64 (loop_pred, cn1ld, 0);
        svbool_t cond2_pred = svcmpeq_n_f64 (loop_pred, cn2ld, 0);
        svbool_t cond_pred = svorr_b_z (loop_pred, cond1_pred, cond2_pred);
        svbool_t not_cond_pred = svnot_b_z (loop_pred, cond_pred);
        if (!svptest_any(all_true, not_cond_pred)) {
            // boscc guarded branch, will skip multiple vector instructions if this is taken
            // works well if branch probability is high -> pgo or estimate statically

            // if
            svst1_f64 (cond_pred, atype+i, svdup_f64(0) );
            printf ("boscc reg kern took boscc branch1\n");
        }
        else {
            // IF CONVERTED CODE
            // Only one more statement is added
            // Can take this into account when building cost model
            //      i.e. boscc guard only one CFG path and IF convert other
            //          or boscc guard multiple CFG paths
            //          Moll et al do a 'trickle down' boscc -> skip as the predicate becomes all false

            // if
            svst1_f64 (cond_pred, atype+i, svdup_f64(0) );

            svfloat64_t cn1ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn1, iaci_v);
            svfloat64_t cn2ld = svld1_gather_s64offset_f64 (loop_pred, prm->Cn2, iaci_v);

            // else
            svfloat64_t sigmaw3_v = svdup_n_f64(sigmaw3);
            svfloat64_t div = svdiv_f64_z (not_cond_pred, cn1ld, cn2ld);
            svfloat64_t sigma_iw6_v = svmul_f64_z (not_cond_pred, sigmaw3_v, svsqrt_f64_z (not_cond_pred, div));
            svfloat64_t epsilon_iw_v = svmul_f64_z(not_cond_pred, svmul_f64_z (not_cond_pred, svdup_n_f64 (0.5),svsqrt_f64_z(not_cond_pred, svdiv_f64_z (not_cond_pred, svdup_n_f64(EPSILONW),cn1ld))),cn2ld);
            svfloat64_t st = svdiv_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svmul_f64_z(not_cond_pred, svdup_n_f64(-16),svdup_n_f64(PI)),svdup_n_f64(RHOW)),epsilon_iw_v),sigma_iw6_v) , svdup_n_f64(3.0) );
            svst1_f64 (not_cond_pred, atype+i, st);
        }

        i_v = svadd_n_s64_z (loop_pred, i_v, svcntd());
        i += svcntd();
        loop_pred = svwhilelt_b64_s64 (i, prm->Ntypes);
    }

    printf ("%f\n", atype[0]);
}

void scalar_kern() {
    for (i = 0; i < prm->Ntypes; i++) {
        iaci = prm->Cno[prm->Ntypes * i + i] - 1; 
        if (prm->Cn1[iaci] == 0.0 || prm->Cn2[iaci] == 0.0) {
            atype[i] = 0.0; 
        } else {
            sigma_iw6 = sigmaw3 * sqrt(prm->Cn1[iaci] / prm->Cn2[iaci]);
            epsilon_iw =
                0.5 * sqrt(EPSILONW / prm->Cn1[iaci]) * prm->Cn2[iaci];
            atype[i] = -16. * PI * RHOW * epsilon_iw * sigma_iw6 / 3.;
        }    
    }
    printf ("%f\n", atype[0]);
}

int main(void) {

    unknown_func(prm);

    //scalar_kern();
    ifconversion_kern();
    //boscc_kern();
    //boscc_unrolled_kern();

    printf ("%f\n", atype[0]);

    return 0;
}