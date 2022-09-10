#define NUM_NU  (3)     /* numero de neutrinos */
#define T_FINAL (1.0)

/* mixing */                    /* http://www.nu-fit.org/?q=node/238 */
#define DM2_21      (7.42e+5)
#define DM2_32      (2.515e+3)
#define THETA12     (33.44)
#define THETA23     (49.2)      /* todos os angulos em graus */
#define THETA13     (8.57)
#define DELTACP     (195.0)

/* precision */
#define PASSO   (1e-4)
#define EPS_ABS (1e-5)  /* erro absoluto para a rotina do GSL */
#define EPS_REL (1e-5)  /* erro relativo para a rotina do GSL */

/* initial condition */
#define RE1     (1.0)
#define RE2     (0.0)
#define RE3     (0.0)
#define IM1     (0.0)
#define IM2     (0.0)
#define IM3     (0.0)

/* constants */
#define G_F     (1.1663787e-11)     /* constante de Fermi */
#define N_A     (6.02214076e+23)    /* constante de Avogadro */
