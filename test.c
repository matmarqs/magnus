#include "main.h"
#include "param.h"

#define CMGET(M, i, j)      gsl_matrix_complex_get((M), (i), (j))
#define CMSET(M, i, j, x)   gsl_matrix_complex_set((M), (i), (j), (x))

int main() {

    /**************/
    /* ALLOCATING */
    /**************/
    double *x = malloc(5 * sizeof(double)), *Ne = malloc(5 * sizeof(double));
    for (size_t i = 0; i < 5; i++) {
        Ne[i] = x[i] = (double) i;
    }
    Space *space = init_space(x, Ne, 5);
    gsl_vector_complex *psi = space->psi;
    gsl_vector_complex_set(psi, 0, gsl_complex_rect(1/sqrt(2), 0.0));
    gsl_vector_complex_set(psi, 1, gsl_complex_rect(1/sqrt(3), 0.0));
    gsl_vector_complex_set(psi, 2, gsl_complex_rect(1/sqrt(6), 0.0));
    gsl_matrix *A = gsl_matrix_alloc(DIM, DIM),
               *B = gsl_matrix_alloc(DIM, DIM),
               *C = gsl_matrix_alloc(DIM, DIM);

    //gsl_matrix *A = space->Omega4;
    //gsl_matrix_set_zero(A);

    /*************/
    /*  TESTING  */
    /*************/
    MSET(A, 0, 0,  1.0); MSET(A, 0, 1,  1.0); MSET(A, 0, 2,  1.0);
    MSET(A, 1, 0,  1.0); MSET(A, 1, 1,  1.0); MSET(A, 1, 2,  1.0);
    MSET(A, 2, 0,  1.0); MSET(A, 2, 1,  1.0); MSET(A, 2, 2,  1.0);
//
    MSET(B, 0, 0,  1.0); MSET(B, 0, 1,  2.0); MSET(B, 0, 2,  1.0);
    MSET(B, 1, 0,  2.0); MSET(B, 1, 1,  2.0); MSET(B, 1, 2,  2.0);
    MSET(B, 2, 0,  1.0); MSET(B, 2, 1,  2.0); MSET(B, 2, 2,  4.0);
//
    comm(A, B, C);
//
    printf("A =\n");
    print_matrix(A);
    printf("\n");
    printf("B =\n");
    print_matrix(B);
    printf("\n");
    printf("C = [A, B] =\n");
    print_matrix(C);
    printf("\n");

    ////print_vec(psi, DIM);
    //expi_matrix_vec(A, 2.0, space);
    ////print_matrix(A, DIM);
    //printf("\n");
    //printf("norm[exp(i A t) psi] = %e\n\n", gsl_blas_dznrm2(psi));
    ////printf("exp(i A t) psi =\n\n");
    ////print_vec(psi, DIM);

    /***********/
    /* FREEING */
    /***********/
    gsl_matrix_free(A); gsl_matrix_free(B); gsl_matrix_free(C);
    free_space(space);
    free(x); free(Ne);
    return 0;

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
    space->eig->cval = gsl_vector_alloc(DIM);
    space->eig->cvec = gsl_matrix_complex_alloc(DIM, DIM);
    space->eig->cwork = gsl_eigen_hermv_alloc(DIM);
    /* scaling and interpolating data */
    space->interp->spline = gsl_spline_alloc(gsl_interp_steffen, N);  /* STEFFEN */
    sqrt2_GF_NA(Ne, N);
    gsl_spline_init(space->interp->spline, x, Ne, N);
    /* defining the parameters */
    //double th12 = DEG_TO_RAD(THETA12),
    //       //th23,
    //       //d_CP,
    //       th13;
    ///* checking number of neutrinos and setting parameters accordingly */
    //if (NUM_NU == 2) {
    //    th13 /*= th23 = d_CP*/ = 0.0;
    //}
    //else {
    //    th13 = DEG_TO_RAD(THETA13);// th23 = DEG_TO_RAD(THETA23); d_CP = DEG_TO_RAD(DELTACP);
    //}
    ///* calculating mixing matrix */
    //double s12 = sin(th12), c12 = cos(th12),
    //       //s23 = sin(th23), c23 = cos(th23),
    //       s13 = sin(th13), c13 = cos(th13);
    /* calculating mixing matrix */
    double s12 = sqrt(0.308) , c12 = sqrt(1 - s12*s12),
           s13 = sqrt(0.0234), c13 = sqrt(1 - s13*s13);
    double
    W11=c13*c13*c12*c12,  W12=c12*s12*c13*c13,  W13=c12*c13*s13,
    W21=c12*s12*c13*c13,  W22=s12*s12*c13*c13,  W23=s12*c13*s13,
    W31=c12*s13*c13    ,  W32=s12*c13*s13    ,  W33=s13*s13;
    gsl_matrix *W = gsl_matrix_alloc(DIM, DIM); /* DIM = 3 */
    MSET(W, 0, 0, W11);  MSET(W, 0, 1, W12);  MSET(W, 0, 2, W13);
    MSET(W, 1, 0, W21);  MSET(W, 1, 1, W22);  MSET(W, 1, 2, W23);
    MSET(W, 2, 0, W31);  MSET(W, 2, 1, W32);  MSET(W, 2, 2, W33);
    space->W = W;
    /* we have H0 = (1 / (E in MeV) ) . diag(-b, 0, a) */
    space->H0 = gsl_matrix_alloc(DIM, DIM); gsl_matrix_set_zero(space->H0);
    space->a = malloc(sizeof(double)); space->b = malloc(sizeof(double));
    //*space->a = DM2_32 * (1e-6 / 2.0) * EV_CM * R_SUN;
    //*space->b = DM2_21 * (1e-6 / 2.0) * EV_CM * R_SUN;
    *space->a = 4.35196e+06;
    *space->b = 0.030554;
    space->commH0_W = gsl_matrix_alloc(DIM, DIM);
    space->Omega2 = gsl_matrix_alloc(DIM, DIM);
    space->Omega4 = gsl_matrix_complex_alloc(DIM, DIM);
    /* setting electron neutrino initial state */
    space->elec = gsl_vector_complex_alloc(DIM);
    VSET(space->elec, 0, gsl_complex_rect(c12*c13, 0.0));
    VSET(space->elec, 1, gsl_complex_rect(s12*c13, 0.0));
    VSET(space->elec, 2, gsl_complex_rect(s13, 0.0));
    space->psi = gsl_vector_complex_alloc(DIM);
    space->psi_aux = gsl_vector_complex_alloc(DIM);
    return space;
}


/* frees workspace */
void free_space(Space *space) {
    /* vectors and matrices */
    gsl_matrix_free(space->H0); gsl_matrix_free(space->W);
    gsl_matrix_free(space->commH0_W);
    gsl_matrix_free(space->Omega2); gsl_matrix_complex_free(space->Omega4);
    gsl_vector_complex_free(space->psi);
    gsl_vector_complex_free(space->psi_aux);
    gsl_vector_complex_free(space->elec);
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
    gsl_vector_free(space->eig->cval);
    gsl_matrix_complex_free(space->eig->cvec);
    gsl_eigen_hermv_free(space->eig->cwork);
    free(space->eig);
    /* last */
    free(space);
}


/* H0 = (1 / (E in MeV) ) . diag(-b, 0, a) */
void setH0(Space *S, double E) {
    MSET(S->H0, 0, 0, -*S->b / E);
    MSET(S->H0, 2, 2, MASS_ORDER == 'N' ? (*S->a / E) : (-*S->a / E));
}


/* apply order M2 or M4 */
void step(double t, double h, int order, Space *space, function *f) {
    if (order == 2) {
        m2(t, h, space, f);    /* Magnus 2 */
        expi_matrix_vec(space->Omega2, -h, space);  /* psi <- exp(-i Omega2 h) psi */
    }
    else {
        m4(t, h, space, f);   /* Magnus 4 */
        expi_cmatrix_vec(space->Omega4, -h, space); /* psi <- exp(-i Omega4 h) psi */
    }
}

/* Omega2 matrix for Magnus Expansion of order 2 */
void m2(double t, double h, Space *space, function *f) {
    double f_bar = f(t + h*0.5, space->interp);   /* exponential midpoint rule */
    gsl_matrix_memcpy(space->Omega2, space->W); /* W */
    gsl_matrix_scale(space->Omega2, f_bar);     /* f_bar W */
    gsl_matrix_add(space->Omega2, space->H0);   /* H0 + f_bar W */
}


/* Omega4 matrix for Magnus Expansion of order 4 */
void m4(double t, double h, Space *space, function *f) {
    double f_pls = f(t + (1.0 + 1.0/SQRT3)*h*0.5, space->interp),
           f_min = f(t + (1.0 - 1.0/SQRT3)*h*0.5, space->interp);
    for (size_t i = 0; i < space->H0->size1; i++)
        for (size_t j = 0; j < space->H0->size2; j++)
            gsl_matrix_complex_set(space->Omega4, i, j,
                gsl_complex_rect(
                    (MGET(space->H0, i, j) + 0.5*(f_pls + f_min)*MGET(space->W, i, j)) * h,
                    (SQRT3 / 12.0) * (f_pls - f_min) * MGET(space->commH0_W, i, j) * h * h
                )
            );
    printf("\n");
    print_cmatrix(space->Omega4);
    printf("\n");
}


/* psi <- exp(i t A) . psi */
void expi_matrix_vec(gsl_matrix *A, double t, Space *S) {
    gsl_eigen_symmv(A, S->eig->val, S->eig->vec, S->eig->work);
    /* U = eigvec, A = U D U^T */
    realmatrix_trans_complexvec(S->eig->vec, S->psi, S->psi_aux);   /* psi_aux = U^T psi */
    for (size_t i = 0; i < A->size1; i++) {  /* psi_aux = e^(iDt) U^T . psi */
        VSET(S->psi_aux, i,
            gsl_complex_mul(
                VGET(S->psi_aux, i),
                gsl_complex_polar(1.0, t * gsl_vector_get(S->eig->val, i))
            )
        );
    }
    realmatrix_complexvec(S->eig->vec, S->psi_aux, S->psi); /* psi = U . psi_aux */
}


/* psi <- exp(i t A) . psi */
void expi_cmatrix_vec(gsl_matrix_complex *A, double t, Space *S) {
    gsl_eigen_hermv(A, S->eig->cval, S->eig->cvec, S->eig->cwork);
    /* U = eigvec, A = U D U^dagger */
    gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, S->eig->cvec,   /* psi_aux = U^dagger psi */
                           S->psi, GSL_COMPLEX_ZERO, S->psi_aux);
    for (size_t i = 0; i < A->size1; i++) {  /* psi_aux = e^(iDt) U^dagger . psi */
        VSET(S->psi_aux, i,
            gsl_complex_mul(
                VGET(S->psi_aux, i),
                gsl_complex_polar(1.0, t * gsl_vector_get(S->eig->cval, i))
            )
        );
    }
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, S->eig->cvec,   /* psi = U . psi_aux */
                     S->psi_aux, GSL_COMPLEX_ZERO, S->psi);
}


/* apply real matrix A to complex vector x and stores in y = A.x */
void realmatrix_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    size_t i, j; gsl_complex z;
    for (i = 0; i < A->size1; i++) {
        z = gsl_complex_mul_real(VGET(x, 0), MGET(A, i, 0));
        for (j = 1; j < A->size2; j++)
            z = gsl_complex_add(z, gsl_complex_mul_real(VGET(x, j), MGET(A, i, j)));
        VSET(y, i, z);
    }
}


/* apply real matrix A^T (transpose) to complex vector x and stores in y = A^T.x */
void realmatrix_trans_complexvec(gsl_matrix *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    size_t i, j; gsl_complex z;
    for (i = 0; i < A->size2; i++) {
        z = gsl_complex_mul_real(VGET(x, 0), MGET(A, 0, i));
        for (j = 1; j < A->size1; j++)
            z = gsl_complex_add(z, gsl_complex_mul_real(VGET(x, j), MGET(A, j, i)));
        VSET(y, i, z);
    }
}


/* stores C = [A, B] = AB - BA */
void comm(gsl_matrix *A, gsl_matrix *B, gsl_matrix *C) {
    size_t i, k, j; double x;
    for (i = 0; i < A->size1; i++)
        for (k = 0; k < B->size2; k++) {
            x = MGET(A, i, 0)*MGET(B, 0, k) - MGET(B, i, 0)*MGET(A, 0, k);
            for (j = 1; j < A->size2; j++)
                x += MGET(A, i, j)*MGET(B, j, k) - MGET(B, i, j)*MGET(A, j, k);
            MSET(C, i, k, x);
        }
}


/* v = sqrt(2) . G_F . N_e */
double v(double t, void *interpo) {  /* void pointer because this defines a gsl_function */
    interpol *interp = (interpol *) interpo;
    return gsl_spline_eval(interp->spline, t, interp->acc);
}


double ne_exp(double t, void *params) {
    (void) params;  /* avoid unused parameter warning */
    return 6.5956e+04 * exp(-10.54 * t);
}


/* multiply our N_e data by sqrt(2) . G_F */
void sqrt2_GF_NA(double *Ne, int N) {
    for (int i = 0; i < N; i++)
        Ne[i] *= M_SQRT2 * G_FxN_A;
}


double amp2(gsl_vector_complex *psi, size_t i) {
    return gsl_complex_abs2(VGET(psi, i));
}


/* probability of an interaction neutrino */
double surv(gsl_vector_complex *psi, gsl_vector_complex *neutrino) {
    double s = gsl_complex_abs2(gsl_complex_mul(VGET(neutrino, 0), VGET(psi, 0)));
    for (size_t j = 1; j < psi->size; j++)
        s += gsl_complex_abs2(gsl_complex_mul(VGET(neutrino, j), VGET(psi, j)));
    return s;
}


/* reads data from file, storages them in (x_ptr, y_ptr) and returns N */
int readalloc(FILE *stream, double **x_ptr, double **y_ptr, int chunk) {
    char *line = NULL;
    size_t len = 0;
    int i = 0, N = 0;
    double *x = (double *) malloc(chunk * sizeof(double));
    double *y = (double *) malloc(chunk * sizeof(double));
    /* getting input */
    while (getline(&line, &len, stream) != -1) {
        if (line[0] == '#')
            continue;
        else {
            sscanf(line, "%lf%lf", &x[N], &y[N]);
            i++; N++;
            if (i > chunk-1) {
                x = (double *) realloc(x, (N+chunk)*sizeof(double));
                y = (double *) realloc(y, (N+chunk)*sizeof(double));
                i = 0;
            }
        }
    }
    /* resizing the arrays to correct size */
    *x_ptr = (double *) realloc(x, N * sizeof(double));
    *y_ptr = (double *) realloc(y, N * sizeof(double));
    /* freeing resources */
    free(line);
    return N;
}


/* print matrix for debuggig */
void print_matrix(gsl_matrix *M) {
    for (size_t i = 0; i < M->size1; i++)  /* OUT OF RANGE ERROR */
        for (size_t j = 0; j < M->size2; j++)
            printf("%+.3e %c", gsl_matrix_get(M, i, j),
                  (j == M->size2-1) ? '\n' : ' ');
}


/* print matrix for debuggig */
void print_cmatrix(gsl_matrix_complex *M) {
    for (size_t i = 0; i < M->size1; i++)  /* OUT OF RANGE ERROR */
        for (size_t j = 0; j < M->size2; j++)
            printf("%+.3e + %+.3e i  %c", GSL_REAL(gsl_matrix_complex_get(M, i, j)),
                                          GSL_IMAG(gsl_matrix_complex_get(M, i, j)),
                                          (j == M->size2-1) ? '\n' : ' ');
}


/* print matrix for debuggig */
void print_vec(gsl_vector_complex *psi) {
    for (size_t j = 0; j < psi->size; j++)
        printf("%+.3e + %+.3e i  %c", GSL_REAL(gsl_vector_complex_get(psi, j)),
                                      GSL_IMAG(gsl_vector_complex_get(psi, j)),
                                      (j == psi->size-1) ? '\n' : ' ');
}

/***************/
/* unused code */
/***************/

////gsl_blas_zdotc(space->psi, space->elec, &z);
//printf("%*ld%15.5e%15.5e\n", it_width,
//        i, ti, //gsl_blas_dznrm2(space->psi),
//        ///*amp2(space->psi, 0),*/ GSL_REAL(VGET(space->psi, 0)), GSL_IMAG(VGET(space->psi, 0)),   /* psi_1 */
//        ///*amp2(space->psi, 1),*/ GSL_REAL(VGET(space->psi, 1)), GSL_IMAG(VGET(space->psi, 1)),   /* psi_2 */
//        ///*amp2(space->psi, 2) */ GSL_REAL(VGET(space->psi, 2)), GSL_IMAG(VGET(space->psi, 2)),   /* psi_3 */
//        gsl_blas_dznrm2(space->psi)   /* norm of psi */
//        //gsl_complex_abs2(z)
//);
