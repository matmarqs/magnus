#include "main.h"
#include "param.h"

/* TODO: expi_matrix_vec nao conserva norma do vetor */

/* DEFINITIONS */

int main(int argc, char *argv[]) {

                            /*************************/
                            /***   READING DATA    ***/
                            /*************************/
    FILE *elecdens_file = fopen("elecdens.txt", "r");
    FILE *energy_file = fopen("8b-energy.txt", "r");
    FILE *distr_file = fopen("8b-distr.txt", "r");
    double *x, *Ne;
    int N = readalloc(elecdens_file, &x, &Ne, 2500); /* we have 2458 lines */
    double *E, *p_E;
    double *r0, *p_r0;
    /*int num_E = */readalloc(energy_file, &E, &p_E, 900);   /* 844 lines */
    /*int num_r = */readalloc(distr_file, &r0, &p_r0, 1300); /* 1226 lines */
    fclose(elecdens_file); fclose(energy_file); fclose(distr_file);


                            /*************************/
                            /***   INITIALIZING    ***/
                            /*************************/
    /* generating all parameters and allocating memory */
    Space *space = init_space(x, Ne, N);

    /* declaring */
    double t0 = T_0, ti;
    //char *format = "%15.5e%15.5e%15.5e%15.5e\n";


                            /*************************/
                            /***        ODE        ***/
                            /*************************/
    long num_it = lround((T_FINAL - t0) / PASSO);   /* number of iterations */
    int it_width = (int) log10(num_it) + 1;   /* variable to format and print things */
    double energy = 6.44;    /* energia mais provavel em MeV */
    setH0(space, energy);
    for (long i = 0; i <= 2; i++) {
        ti = i * PASSO + t0;    /* constant step size (variable will be implemented) */
        printf("%*ld%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n", it_width,
                i, ti,
                GSL_REAL(VGET(space->psi, 0)), GSL_IMAG(VGET(space->psi, 0)),   /* psi_1 */
                GSL_REAL(VGET(space->psi, 1)), GSL_IMAG(VGET(space->psi, 1)),   /* psi_2 */
                GSL_REAL(VGET(space->psi, 2)), GSL_IMAG(VGET(space->psi, 2)),   /* psi_3 */
                gsl_blas_dznrm2(space->psi));   /* norm of psi */
        m2(ti, PASSO, space);   /* Magnus 2 */
        expi_matrix_vec(space->Omega2, -PASSO, space);   /* psi <- exp(-i Omega2 h) psi */
    }
    printf("energy = %.5e\n", energy);
    /* this printf below is for the survival probability */
    //printf(format, t0, energy, surv(space->psi), gsl_blas_dznrm2(space->psi));


                            /*************************/
                            /*   FREEING RESOURCES   */
                            /*************************/
    /* space */
    free_space(space);
    /* data */
    free(x); free(Ne);
    free(E); free(p_E);
    free(r0); free(p_r0);
    return 0;

}


/* gera o workspace e retorna o pointer dele */
Space *init_space(double *x, double *Ne, int N) {
    /* initializing lower-level things */
    Space *space = malloc(sizeof(Space));
    space->interp = malloc(sizeof(interpol));
    space->interp->acc = gsl_interp_accel_alloc();
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
    space->Apsi = gsl_vector_complex_alloc(DIM);
    space->AApsi = gsl_vector_complex_alloc(DIM);
    return space;
}


/* frees workspace */
void free_space(Space *space) {
    /* vectors and matrices */
    gsl_matrix_free(space->H0); gsl_matrix_free(space->W);
    gsl_matrix_free(space->Omega2); gsl_matrix_free(space->Omega4);
    gsl_vector_complex_free(space->psi);
    gsl_vector_complex_free(space->Apsi);
    gsl_vector_complex_free(space->AApsi);
    /* parameters */
    free(space->a); free(space->b);
    /* interpol */
    gsl_interp_accel_free(space->interp->acc);
    gsl_spline_free(space->interp->spline);
    free(space->interp);
    /* last */
    free(space);
}


/* H0 = (1 / (E in MeV) ) . diag(-b, 0, a) */
void setH0(Space *S, double E) {
    MSET(S->H0, 0, 0, -*S->b / E);
    MSET(S->H0, 2, 2, MASS_ORDER == 'N' ? (*S->a / E) : (-*S->a / E));
}


/* Omega2 matrix for Magnus Expansion of order 2 */
void m2(double t, double h, Space *space) {
    double v_bar = v(t + h/2, space->interp);   /* exponential midpoint rule */
    gsl_matrix_memcpy(space->Omega2, space->W); /* W */
    gsl_matrix_scale(space->Omega2, v_bar);     /* v_bar W */
    gsl_matrix_add(space->Omega2, space->H0);   /* H0 + vbar W */
}


/* expA <- exp(i t A) */
void expi_matrix_vec(gsl_matrix *A, double t, Space *S) {
    /* tr3 = tr(A)/3 */
    double tr3 = (MGET(A,0,0) + MGET(A,1,1) + MGET(A,2,2))/3.0;
    /* making A traceless */
    MSET(A, 0, 0, MGET(A, 0, 0) - tr3);
    MSET(A, 1, 1, MGET(A, 1, 1) - tr3);
    MSET(A, 2, 2, MGET(A, 2, 2) - tr3);
    /* calculating base parameters q, p */
    double q, p;
    get_qp(A, &q, &p);
    /* eigenvalues of matrix A */
    double l0, l1, l2;
    gsl_poly_solve_cubic(0.0, -p, q, &l0, &l1, &l2);
    //printf("l0, l1, l2 = %e, %e, %e\n", l0, l1, l2);
    if (!(l0 == l1 && l0 == l2)) {
        /* eigenvalues differences */
        double a = l1 - l0;
        double b = l2 - l0;
        /* calculating A.psi and A^2.psi */
        realmatrix_complexvec(A, S->psi, S->Apsi);
        realmatrix_complexvec(A, S->Apsi, S->AApsi);
        /* parameters r0 and r1 */
        gsl_complex r0, r1, r_aux;
        /* a and b can not be zero at the same time */
        r0    = (a == 0) ? gsl_complex_rect(0.0, t) : exp1i(a, t);
        r_aux = (b == 0) ? gsl_complex_rect(0.0, t) : exp1i(b, t);
        r1 = (a == b)
           ? gsl_complex_mul(
              gsl_complex_polar(1/(b*b), b*t),
              gsl_complex_rect(-1.0, b*t)
             )
           : gsl_complex_div_real(
               gsl_complex_sub(
                 r0, r_aux
               ),
               a-b
             );
        /* scaling psi, A.psi and A^2.psi */
        gsl_vector_complex_scale(S->psi,
            gsl_complex_sub(GSL_COMPLEX_ONE,
                gsl_complex_mul_real(
                    gsl_complex_sub(r0, gsl_complex_mul_real(r1, l1)),
                    l0
                )
            )
        );
        gsl_vector_complex_scale(S->Apsi,
            gsl_complex_add(r0, gsl_complex_mul_real(r1, l2))
        );
        gsl_vector_complex_scale(S->AApsi, r1);
        /* adding all contributions */
        gsl_vector_complex_add(S->psi, S->Apsi);
        gsl_vector_complex_add(S->psi, S->AApsi);
    }
    /* final scaling */
    gsl_vector_complex_scale(S->psi,
        gsl_complex_polar(1.0, (tr3 + l0) * t)
    );
}


/* calculating q, p and tr3 parameters to exponentiate a matrix */
void get_qp(gsl_matrix *A, double *q, double *p) {
    /* q = det(A) */
    *q = MGET(A,0,0) * ( MGET(A,1,1)*MGET(A,2,2) - MGET(A,1,2)*MGET(A,2,1) )
       - MGET(A,0,1) * ( MGET(A,1,0)*MGET(A,2,2) - MGET(A,1,2)*MGET(A,2,0) )
       + MGET(A,0,2) * ( MGET(A,1,0)*MGET(A,2,1) - MGET(A,1,1)*MGET(A,2,0) );
    /* p = tr(A^2)/2 */
    *p = (MGET(A,0,0) * MGET(A,0,0) + MGET(A,1,1) * MGET(A,1,1) + MGET(A,2,2) * MGET(A,2,2))/2.0
       +  MGET(A,0,1) * MGET(A,1,0) + MGET(A,1,2) * MGET(A,2,1) + MGET(A,2,0) * MGET(A,0,2);
}


/*               exp(i x t) - 1  */
/*  exp1i(x, t) = --------------  */
/*                      x        */
gsl_complex exp1i(double x, double t) {
    return gsl_complex_div_real(
               gsl_complex_sub(
                   gsl_complex_polar(1.0, x*t), GSL_COMPLEX_ONE
               ),
               x
           );
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


/* v = sqrt(2) . G_F . N_e */
double v(double t, void *interpo) {  /* void pointer because this defines a gsl_function */
    interpol *interp = (interpol *) interpo;
    return gsl_spline_eval(interp->spline, t, interp->acc);
}


/* multiply our N_e data by sqrt(2) . G_F */
void sqrt2_GF_NA(double *Ne, int N) {
    for (int i = 0; i < N; i++)
        Ne[i] *= M_SQRT2 * G_FxN_A;
}


/* survival probability of the electron neutrino */
double surv(gsl_vector_complex *psi) {
    return gsl_complex_abs2(VGET(psi, 1));
}


/* check if norm = 1, only for debugging */
double norm(gsl_vector_complex *psi) {
    return gsl_blas_dznrm2(psi);
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
void print_matrix(gsl_matrix *M, size_t dim) {
    for (size_t i = 0; i < dim; i++)  /* OUT OF RANGE ERROR */
        for (size_t j = 0; j < dim; j++)
            printf("%e%c", gsl_matrix_get(M, i, j),
                  (j == dim-1) ? '\n' : ' ');
}
