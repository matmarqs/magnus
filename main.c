#include "main.h"
#include "param.h"


/* v = sqrt(2) . G_F . N_e */
double v(double t, void *interpo) {  /* void pointer because this defines a gsl_function */
    interpol *interp = (interpol *) interpo;
    return gsl_spline_eval(interp->spline, t, interp->acc);
}


void m2(double t, double h, Space *space) {
    double v_bar = v(t + h/2, space);
    gsl_matrix_memcpy(space->Omega2, space->W); /* W */
    gsl_matrix_scale(space->Omega2, v_bar);     /* v_bar W */
    gsl_matrix_add(space->Omega2, space->H0);   /* H0 + vbar W */
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

gsl_complex exp1(double x, double t) {
    return gsl_complex_div_real(
             gsl_complex_sub(
               gsl_complex_polar(1.0, x*t), GSL_COMPLEX_ONE
             ),
             x
           );
}

/* expA <- exp(tA) */
void matrix_exp_vec(gsl_matrix *A, double t, Space *S) {
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

    gsl_complex r0, r1;
    r0 = exp1(a, t);
    r1 = gsl_complex_div_real(
           gsl_complex_sub(
             exp1(a, t), exp1(b, t)
           ),
           a-b
         );

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

    gsl_vector_complex_add(S->psi, S->Apsi);
    gsl_vector_complex_add(S->psi, S->AApsi);

    gsl_vector_complex_scale(S->psi,
        gsl_complex_polar(1.0, (tr3 + l0) * t)
    );
}

/* check if norm = 1, only for debugging */
double norm(gsl_vector_complex *psi) {
    return gsl_blas_dznrm2(psi);
}

double surv(gsl_vector_complex *psi) {
    return gsl_complex_abs2(VGET(psi, 1));
}

int main(int argc, char *argv[]) {
    FILE *elecdens = fopen("elecdens.txt", "r");
    FILE *energy = fopen("8b-energy.txt", "r");
    FILE *distr = fopen("8b-distr.txt", "r");

                            /*************************/
                            /***   READING DATA    ***/
                            /*************************/
    double *x, *Ne;
    int N = readalloc(elecdens, &x, &Ne, 2500); /* we have 2458 lines */
    sqrt2_GF_NA(Ne, N);
    double *E, *p_E;
    double *r0, *p_r0;
    /*int num_E = */readalloc(energy, &E, &p_E, 900);   /* 844 lines */
    /*int num_r = */readalloc(distr, &r0, &p_r0, 1300); /* 1226 lines */
    fclose(elecdens); fclose(energy); fclose(distr);

                            /*************************/
                            /***   INITIALIZING    ***/
                            /*************************/
    /* generating all parameters and allocating memory */
    Space *S = init_space(x, Ne, N);

    /* declaring */
    double ti, t, t0;
    char *format = "%15.5e%15.5e%15.5e%15.5e\n";


                            /*************************/
                            /***        EDO        ***/
                            /*************************/
    t = t0;;
    double Psi[DIM] = { RE1, RE2, RE3, IM1, IM2, IM3 }; /* initial condition */
    long num_it = lround((T_FINAL - t0) / PASSO);
    //int it_width = (int) log10(num_it) + 1;
    for (long i = 0; i <= num_it; i++) {
        ti = i * PASSO + t0;
        if ((status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi)) != GSL_SUCCESS) {
            printf ("Error, return value = %d\n", status);
            break;
        }
        //printf("%*d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n", it_width,
        //        i,t,Psi[0],Psi[1],Psi[2],Psi[3],Psi[4],Psi[5],norm(Psi));
    }
    printf(format, t0, params->energ, sobrev(Psi), norm(Psi));



                            /*************************/
                            /* FREEING THE RESOURCES */
                            /*************************/
    /* space */
    free_space(space);
    /* data */
    free(x); free(Ne);
    free(E); free(p_E);
    free(r0); free(p_r0);
}

void sqrt2_GF_NA(double *Ne, int N) {
    for (int i = 0; i < N; i++)
        Ne[i] *= M_SQRT2 * G_F * N_A;
}


/* gera os parametros que precisamos e retorna o pointer dele */
Space *init_space(double *x, double *Ne, int N) {
    Space *space = malloc(sizeof(Space));
    space->interp = malloc(sizeof(interpol));
    space->interp->acc = gsl_interp_accel_alloc();
    space->interp->spline = gsl_spline_alloc(gsl_interp_steffen, N);  /* STEFFEN */
    gsl_spline_init(space->interp->spline, x, Ne, N);

    /* defining the parameters */
    double th12 = DEG_TO_RAD(THETA12),
           th23, th13;// double d_CP;

    /* checking number of neutrinos and setting parameters accordingly */
    if (NUM_NU == 2) {
        th23 = th13 = 0.0;// d_CP = 0.0;
    }
    else {
        th23 = DEG_TO_RAD(THETA23); th13 = DEG_TO_RAD(THETA13);// d_CP = DEG_TO_RAD(DELTACP);
    }

    /* calculating mixing matrix */
    double s12 = sin(th12), c12 = cos(th12),
           s23 = sin(th23), c23 = cos(th23),
           s13 = sin(th13), c13 = cos(th13);
    double
    W11=c13*c13*c12*c12,  W12=c12*s12*c13*c13,  W13=c12*c13*s13,
    W21=c12*s12*c13*c13,  W22=s12*s12*c13*c13,  W23=s12*c13*s13,
    W31=c12*s13*c13    ,  W32=s12*c13*s13    ,  W33=s13*s13;

    gsl_matrix *W = gsl_matrix_alloc(SIZE, SIZE);
    MSET(W, 0, 0, W11);  MSET(W, 0, 1, W12);  MSET(W, 0, 2, W13);
    MSET(W, 1, 0, W21);  MSET(W, 1, 1, W22);  MSET(W, 1, 2, W23);
    MSET(W, 2, 0, W31);  MSET(W, 2, 1, W32);  MSET(W, 2, 2, W33);

    /* INCOMPLETO */
    gsl_matrix *H0 = gsl_matrix_alloc(SIZE, SIZE); gsl_matrix_set_zero(H0);
    MSET(H0, 0, 0, -DM2_21);
                              MSET(H0, 1, 1, 0.0);
                                                    MSET(H0, 2, 2, DM2_32);

    space->W = W;
    space->H0 = H0;
    space->Omega2 = gsl_matrix_alloc(SIZE, SIZE);
    space->Omega4 = gsl_matrix_alloc(SIZE, SIZE);
    space->psi = gsl_vector_complex_alloc(SIZE);
    space->Apsi = gsl_vector_complex_alloc(SIZE);
    space->AApsi = gsl_vector_complex_alloc(SIZE);

    return space;
}

void free_space(Space *space) {
    /* vectors and matrices */
    gsl_matrix_free(space->H0); gsl_matrix_free(space->W);
    gsl_matrix_free(space->Omega2); gsl_matrix_free(space->Omega4);
    gsl_vector_complex_free(space->psi);
    gsl_vector_complex_free(space->Apsi);
    gsl_vector_complex_free(space->AApsi);
    /* interpol */
    gsl_interp_accel_free(space->interp->acc);
    gsl_spline_free(space->interp->spline);
    /* space components */
    free(space->interp);
    free(space);
}

/* reads data from standard input, storages them in x_ptr and y_ptr and returns N */
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
