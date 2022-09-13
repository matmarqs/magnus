#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_vector.h>         /* vector definitions */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
//#include <gsl/gsl_linalg.h>         /* linear algebra */
#include <gsl/gsl_spline.h>         /* interpolation of real data */
#include <gsl/gsl_poly.h>           /* find cubic roots */
#include <gsl/gsl_eigen.h>          /* solve eigensystems */

#define DIM 3
#define VGET(V, i)          gsl_vector_complex_get((V), (i))
#define VSET(V, i, x)       gsl_vector_complex_set((V), (i), (x))
#define MGET(M, i, j)       gsl_matrix_get((M), (i), (j))
#define MSET(M, i, j, x)    gsl_matrix_set((M), (i), (j), (x))
#define DEG_TO_RAD(ang)     ((M_PI / 180.0) * ang)

#define NUM_NU      (3)     /* numero de neutrinos */
#define T_0         (0.0)
#define T_FINAL     (1.0)
#define MASS_ORDER  ('N')   /* 'N' for normal, 'I' for inverted */

/* mixing */                    /* http://www.nu-fit.org/?q=node/238 */
#define DM2_21      (7.42e-5)
#define DM2_32      (2.515e-3)
#define THETA12     (33.44)
#define THETA23     (49.2)      /* todos os angulos em graus */
#define THETA13     (8.57)
#define DELTACP     (195.0)

/* precision */
#define PASSO   (1e-3)

/* constants */
#define G_F     (1.1663787e-23) /* constante de Fermi em eV^{-2} */
#define N_A     (6.02214076e+23)    /* numero de Avogadro */
#define G_FxN_A (7.0240967108658126)    /* Fermi x Avogadro em eV^{-2} */
#define EV_CM   (50677.30716548338) /* eV.cm = 1 / [ 100 * (c in m/s) * (hbar in eV.s) ] */
#define R_SUN   (6.957e+10) /* R_SUN = 6.957e+10 centimeters */


typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
} interpol;

typedef struct {
    gsl_matrix *vec;
    gsl_vector *val;
    gsl_eigen_symmv_workspace *work;
} eigen_sys;

typedef struct {
    gsl_matrix *H0, *W;
    gsl_matrix *Omega2, *Omega4;
    gsl_vector_complex *psi, *psi_aux;
    interpol *interp;
    double *a, *b;
    eigen_sys *eig;
} Space;


/* DEFINITIONS */
void realmatrix_trans_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y);
void free_space(Space *space);
void sqrt2_GF_NA(double *Ne, int N);
Space *init_space(double *x, double *Ne, int N);
void expi_matrix_vec(gsl_matrix *A, double t, Space *S);
void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y);
void print_matrix(gsl_matrix *M, size_t dim);
void print_vec(gsl_vector_complex *psi, size_t dim);


int main() {
    double *x = malloc(5 * sizeof(double)), *Ne = malloc(5 * sizeof(double));
    for (size_t i = 0; i < 5; i++) {
        Ne[i] = x[i] = (double) i;
    }
    Space *space = init_space(x, Ne, 5);
    gsl_vector_complex *psi = space->psi;
    gsl_vector_complex_set(psi, 0, gsl_complex_rect(1/sqrt(2), 0.0));
    gsl_vector_complex_set(psi, 1, gsl_complex_rect(1/sqrt(3), 0.0));
    gsl_vector_complex_set(psi, 2, gsl_complex_rect(1/sqrt(6), 0.0));
    printf("norm(psi) = %e\n\n", gsl_blas_dznrm2(psi));
    gsl_matrix *A = space->Omega2;
    gsl_matrix_set_zero(A);
    MSET(A, 0, 0, 1.0);    MSET(A, 0, 1, -2.0); MSET(A, 0, 2, 1.0);
    MSET(A, 1, 0, -2.0); MSET(A, 1, 1, 2.0);    MSET(A, 1, 2, -2.0);
    MSET(A, 2, 0, 1.0);   MSET(A, 2, 1, -2.0);   MSET(A, 2, 2, 4.0);
    printf("A =\n");
    print_matrix(A, DIM);
    //printf("\n");
    //print_vec(psi, DIM);
    expi_matrix_vec(A, 2.0, space);
    //print_matrix(A, DIM);
    printf("\n");
    printf("norm[exp(i A t) psi] = %e\n\n", gsl_blas_dznrm2(psi));
    //printf("exp(i A t) psi =\n\n");
    //print_vec(psi, DIM);
    free_space(space);
    free(x); free(Ne);
    return 0;
}


/* print matrix for debuggig */
void print_vec(gsl_vector_complex *psi, size_t dim) {
    for (size_t j = 0; j < dim; j++)
        printf("%e + %e i  %c", GSL_REAL(gsl_vector_complex_get(psi, j)),
                                GSL_IMAG(gsl_vector_complex_get(psi, j)),
                                (j == dim-1) ? '\n' : ' ');
}


/* print matrix for debuggig */
void print_matrix(gsl_matrix *M, size_t dim) {
    for (size_t i = 0; i < dim; i++)  /* OUT OF RANGE ERROR */
        for (size_t j = 0; j < dim; j++)
            printf("%e%c", gsl_matrix_get(M, i, j),
                  (j == dim-1) ? '\n' : ' ');
}


/* apply real matrix A to complex vector x and stores in y = A.x */
void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    size_t i, j; gsl_complex z;
    for (i = 0; i < DIM; i++) {
        z = gsl_complex_mul_real(VGET(x, 0), MGET(A, i, 0));
        for (j = 1; j < DIM; j++)
            z = gsl_complex_add(z, gsl_complex_mul_real(VGET(x, j), MGET(A, i, j)));
        VSET(y, i, z);
    }
}


/* apply real matrix A^T (transpose) to complex vector x and stores in y = A^T.x */
void realmatrix_trans_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    size_t i, j; gsl_complex z;
    for (i = 0; i < DIM; i++) {
        z = gsl_complex_mul_real(VGET(x, 0), MGET(A, 0, i));
        for (j = 1; j < DIM; j++)
            z = gsl_complex_add(z, gsl_complex_mul_real(VGET(x, j), MGET(A, j, i)));
        VSET(y, i, z);
    }
}


/* psi <- exp(i t A) . psi */
void expi_matrix_vec(gsl_matrix *A, double t, Space *S) {
    gsl_eigen_symmv(A, S->eig->val, S->eig->vec, S->eig->work);
    /* U = eigvec, A = U D U^T */
    realmatrix_trans_complexvec(S->eig->vec, S->psi, S->psi_aux);   /* psi_aux = U^T psi */
    for (size_t j = 0; j < DIM; j++) {  /* psi_aux = e^(iDt) U^T . psi */
        VSET(S->psi_aux, j,
            gsl_complex_mul(
                VGET(S->psi_aux, j),
                gsl_complex_polar(1.0, t * gsl_vector_get(S->eig->val, j))
            )
        );
    }
    realmatrix_complexvec(S->eig->vec, S->psi_aux, S->psi); /* psi = U e^(iDt) U^T . psi */
}


/* gera o workspace e retorna o pointer dele */
Space *init_space(double *x, double *Ne, int N) {
    /* initializing lower-level things */
    Space *space = malloc(sizeof(Space));
    space->interp = malloc(sizeof(interpol));
    space->interp->acc = gsl_interp_accel_alloc();
    /* allocating memory for eigensystem */
    space->eig = malloc(sizeof(eigen_sys));
    space->eig->val = gsl_vector_alloc(DIM);
    space->eig->vec = gsl_matrix_alloc(DIM, DIM);
    space->eig->work = gsl_eigen_symmv_alloc(DIM);
    /* scaling and interpolating data */
    space->interp->spline = gsl_spline_alloc(gsl_interp_steffen, N);  /* STEFFEN */
    sqrt2_GF_NA(Ne, N);
    gsl_spline_init(space->interp->spline, x, Ne, N);
    /* defining the parameters */
    double th12 = DEG_TO_RAD(THETA12),
           //th23,
           //d_CP,
           th13;
    /* checking number of neutrinos and setting parameters accordingly */
    if (NUM_NU == 2) {
        th13 /*= th23 = d_CP*/ = 0.0;
    }
    else {
        th13 = DEG_TO_RAD(THETA13);// th23 = DEG_TO_RAD(THETA23); d_CP = DEG_TO_RAD(DELTACP);
    }
    /* calculating mixing matrix */
    double s12 = sin(th12), c12 = cos(th12),
           //s23 = sin(th23), c23 = cos(th23),
           s13 = sin(th13), c13 = cos(th13);
    double
    W11=c13*c13*c12*c12,  W12=c12*s12*c13*c13,  W13=c12*c13*s13,
    W21=c12*s12*c13*c13,  W22=s12*s12*c13*c13,  W23=s12*c13*s13,
    W31=c12*s13*c13    ,  W32=s12*c13*s13    ,  W33=s13*s13;
    gsl_matrix *W = gsl_matrix_alloc(DIM, DIM);
    MSET(W, 0, 0, W11);  MSET(W, 0, 1, W12);  MSET(W, 0, 2, W13);
    MSET(W, 1, 0, W21);  MSET(W, 1, 1, W22);  MSET(W, 1, 2, W23);
    MSET(W, 2, 0, W31);  MSET(W, 2, 1, W32);  MSET(W, 2, 2, W33);
    /* we have H0 = (1 / (E in MeV) ) . diag(-b, 0, a) */
    gsl_matrix *H0 = gsl_matrix_alloc(DIM, DIM); gsl_matrix_set_zero(H0);
    space->a = malloc(sizeof(double)); space->b = malloc(sizeof(double));
    *space->a = DM2_32 * (1e-6 / 2.0) * EV_CM * R_SUN;
    *space->b = DM2_21 * (1e-6 / 2.0) * EV_CM * R_SUN;
    /* allocating memory for workspace */
    space->W = W;
    space->H0 = H0;
    space->Omega2 = gsl_matrix_alloc(DIM, DIM);
    space->Omega4 = gsl_matrix_alloc(DIM, DIM);
    /* setting electron neutrino initial state */
    space->psi = gsl_vector_complex_alloc(DIM);
    VSET(space->psi, 0, gsl_complex_rect(c12*c13, 0.0));
    VSET(space->psi, 1, gsl_complex_rect(s12*c13, 0.0));
    VSET(space->psi, 2, gsl_complex_rect(s13, 0.0));
    space->psi_aux = gsl_vector_complex_alloc(DIM);
    return space;
}


/* multiply our N_e data by sqrt(2) . G_F */
void sqrt2_GF_NA(double *Ne, int N) {
    for (int i = 0; i < N; i++)
        Ne[i] *= M_SQRT2 * G_FxN_A;
}


/* frees workspace */
void free_space(Space *space) {
    /* vectors and matrices */
    gsl_matrix_free(space->H0); gsl_matrix_free(space->W);
    gsl_matrix_free(space->Omega2); gsl_matrix_free(space->Omega4);
    gsl_vector_complex_free(space->psi);
    gsl_vector_complex_free(space->psi_aux);
    /* parameters */
    free(space->a); free(space->b);
    /* interpol */
    gsl_interp_accel_free(space->interp->acc);
    gsl_spline_free(space->interp->spline);
    free(space->interp);
    /* eigensystem */
    gsl_vector_free(space->eig->val);
    gsl_matrix_free(space->eig->vec);
    gsl_eigen_symmv_free(space->eig->work);
    free(space->eig);
    /* last */
    free(space);
}
