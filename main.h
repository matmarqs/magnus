#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_vector.h>         /* vector definitions */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
#include <gsl/gsl_linalg.h>         /* linear algebra */
#include <gsl/gsl_spline.h>         /* interpolation of real data */

/* macros */
#define DIM 3
#define VGET(V, i)          gsl_vector_complex_get((V), (i))
#define VSET(V, i, x)       gsl_vector_complex_set((V), (i), (x))
#define MGET(M, i, j)       gsl_matrix_get((M), (i), (j))
#define MSET(M, i, j, x)    gsl_matrix_set((M), (i), (j), (x))
#define DEG_TO_RAD(ang)     ((M_PI / 180.0) * ang)


typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
} interpol;

typedef struct {
    gsl_matrix *H0, *W;
    gsl_matrix *Omega2, *Omega4;
    gsl_vector_complex *psi, *Apsi, *AApsi;
    interpol *interp;
    double *a, *b;
} Space;

/* DEFINITIONS */
void setH0(Space *S, double E);
void free_space(Space *space);
Space *init_space(double *x, double *Ne, int N);
int readalloc(FILE *stream, double **x_ptr, double **y_ptr, int chunk);
double surv(gsl_vector_complex *psi);
double norm(gsl_vector_complex *psi);
void matrix_exp_vec(gsl_matrix *A, double t, Space *S);
gsl_complex exp1(double x, double t);
void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y);
void order3(double *l0, double *l1, double *l2);
void get_qptr3(gsl_matrix *A, double *q, double *p, double *tr3);
void m2(double t, double h, Space *space);
double v(double t, void *interpo);
void sqrt2_GF_NA(double *Ne, int N);