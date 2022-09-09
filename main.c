#include "main.h"
#include "param.h"

/* A <- A + lambda B */
void mygsl_matrix_addscale(gsl_matrix *A, double lambda , gsl_matrix *B) {
    for (size_t i = 0; i < SIZE; i++)
        for (size_t j = 0; i < SIZE; i++)
            gsl_matrix_set(A, i, j, lambda * gsl_matrix_get(B, i, j));
}

void m2(double t, double h, space *S) {
    double t_bar = t + h/2;
    double v_bar = v(t_bar);
    gsl_matrix_memcpy(S->Omega2, S->H0);    /* H0 */
    mygsl_matrix_addscale(S->Omega2, v_bar, S->W);   /* H0 + vbar W */
}

void get_qptr3(gsl_matrix *A, double *q, double *p, double *tr3) {
    /* q = det(A) */
    *q = MGET(A,0,0) * ( MGET(A,1,1)*MGET(A,2,2) - MGET(A,1,2)*MGET(A,2,1) )
       - MGET(A,0,1) * ( MGET(A,1,0)*MGET(A,2,2) - MGET(A,1,2)*MGET(A,2,0) )
       + MGET(A,0,2) * ( MGET(A,1,0)*MGET(A,2,1) - MGET(A,1,1)*MGET(A,2,0) );

    /* p = tr(A^2)/2 */
    *p = (MGET(A,0,0) * MGET(A,0,0) + MGET(A,1,1) * MGET(A,1,1) + MGET(A,2,2) * MGET(A,2,2))/2.0
       +  MGET(A,0,1) * MGET(A,1,0) + MGET(A,1,2) * MGET(A,2,1) + MGET(A,2,0) * MGET(A,0,2);

    *tr3 = (MGET(A,0,0) + MGET(A,1,1) + MGET(A,2,2))/3.0;
}

void order3(double *l0, double *l1, double *l2) {
    double a;
    if (*l0 > *l1)
        a = *l0, *l0 = *l1, *l1 = a;
    if (*l0 > *l2)
        a = *l0, *l0 = *l2, *l2 = a;
    if (*l1 > *l2)
        a = *l1, *l1 = *l2, *l2 = a;
}

void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    size_t i, j;
    for (i = 0; i < SIZE; i++) {
        VSET(y, i, gsl_complex_mul_real(VGET(x, 0), MGET(A, i, 0)));
        for (size_t j = 1; j < 3; j++)
            VSET(y, i, gsl_complex_add(VGET(y, i), gsl_complex_mul_real(VGET(x, j), MGET(A, i, j))));
    }
}

/* expA <- exp(tA) */
void matrix_exp_vec(gsl_matrix *A, double t, space *S) {
    double q, p, tr3;
    get_qptr3(A, &q, &p, &tr3);
    int sign = (q > 0) ? (-1) : (1);
    double l0 = sign*2*sqrt(p/3.0)*cos((acos(q*sqrt(0.75/p)))/3.0);
    double l1 = sign*2*sqrt(p/3.0)*cos((acos(q*sqrt(0.75/p)))/3.0 - 2.0*M_PI/3.0);
    double l2 = sign*2*sqrt(p/3.0)*cos((acos(q*sqrt(0.75/p)))/3.0 - 4.0*M_PI/3.0);
    order3(&l0, &l1, &l2);

    double a = l1 - l0;
    double b = l2 - l0;

    realmatrix_complexvec(A, S->psi, S->Apsi);
    realmatrix_complexvec(A, S->Apsi, S->AApsi);
}
