/* standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* GSL */
#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_vector.h>         /* vector definitions */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
//#include <gsl/gsl_linalg.h>         /* linear algebra */
#include <gsl/gsl_spline.h>         /* interpolation of real data */
#include <gsl/gsl_poly.h>           /* find cubic roots */
#include <gsl/gsl_eigen.h>          /* solve eigensystems */

/* macros */
#define DIM 3
#define VGET(V, i)          gsl_vector_complex_get((V), (i))
#define VSET(V, i, x)       gsl_vector_complex_set((V), (i), (x))
#define MGET(M, i, j)       gsl_matrix_get((M), (i), (j))
#define MSET(M, i, j, x)    gsl_matrix_set((M), (i), (j), (x))
#define DEG_TO_RAD(ang)     ((M_PI / 180.0) * ang)
#define SQRT3               1.73205080756887729352  /* sqrt(3) */


typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
} interpol;

typedef struct {
    /* real symmetric */
    gsl_matrix *vec;
    gsl_vector *val;
    gsl_eigen_symmv_workspace *work;
    /* complex hermitian */
    gsl_matrix_complex *cvec;
    gsl_vector *cval;
    gsl_eigen_hermv_workspace *cwork;
} eigen_sys;

typedef struct {
    gsl_matrix *H0, *W, *commH0_W;
    gsl_matrix *Omega2;
    gsl_matrix_complex *Omega4;
    gsl_vector_complex *psi, *psi_aux;
    gsl_vector_complex *elec;
    interpol *interp;
    double *a, *b;
    eigen_sys *eig;
} Space;


/* DEFINITIONS */
Space *init_space(double *x, double *Ne, int N);
void free_space(Space *space);
void setH0(Space *S, double E);
void step(double t, double h, Space *space, int method);
void m2(double t, double h, Space *space);
void m4(double t, double h, Space *space);
void expi_matrix_vec(gsl_matrix *A, double t, Space *S);
void expi_cmatrix_vec(gsl_matrix_complex *A, double t, Space *S);
void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y);
void realmatrix_trans_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y);
void comm(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C);
double v(double t, void *interpo);
void sqrt2_GF_NA(double *Ne, int N);
double amp2(gsl_vector_complex *psi, size_t i);
double surv(gsl_vector_complex *psi, gsl_vector_complex *neutrino);
int readalloc(FILE *stream, double **x_ptr, double **y_ptr, int chunk);
void print_matrix(gsl_matrix *M);
void print_vec(gsl_vector_complex *psi);
