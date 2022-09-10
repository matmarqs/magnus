#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_linalg.h>         /* linear algebra */

/* DEFINITIONS */
void magnus2(double t, double h, gsl_matrix_complex *H);

typedef struct {
    gsl_matrix *H0, *W;
    gsl_matrix *Omega2;
    gsl_vector_complex *psi, *Apsi, *AApsi;
} space;

#define VGET(V, i)      gsl_vector_complex_get((V), (i))
#define VSET(V, i, x)   gsl_vector_complex_set((V), (i), (x))
#define MGET(M, i, j)   gsl_matrix_get((M), (i), (j))
#define SIZE 3
#define COMPLEX_1       gsl_complex_rect(1.0, 0.0)
