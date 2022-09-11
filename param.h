#define NUM_NU      (3)     /* numero de neutrinos */
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
#define PASSO   (1e-4)

/* initial condition */
#define RE1     (1.0)
#define RE2     (0.0)
#define RE3     (0.0)
#define IM1     (0.0)
#define IM2     (0.0)
#define IM3     (0.0)

/* constants */
#define G_F     (1.1663787e-11)         /* constante de Fermi */
#define N_A     (6.02214076e+23)        /* constante de Avogadro */
#define G_FxN_A (7.0240967108658126)    /* Fermi x Avogadro */
#define EV_CM   (50677.30716548338) /* eV.cm = 1 / [ 100 * (c in m/s) * (hbar in eV.s) ] */
#define R_SUN   (6.957e+10) /* R_SUN = 6.957e+10 centimeters */

/*
# constantes (todas listadas aqui estao potencias de eV)
sqrt2 = 1.4142135623730951
solar = 1e+3    # plotaremos para r in [0, solar * R_sun]
gamma = 1.23984198e-04  # 1/cm = gamma * eV
# Temos que as constantes G_F e N_A so aparecem juntas da forma G_F * N_A
# G_F = 1.1663787e-23   # eV^{-2}
# N_A = 6.02214076e+23  # Avogadro
G_FxN_A = 7.0240967108658126    # 1.1663787 * 6.02214076
R_sun = (1/solar) * (1/gamma) * 6.957e+10   # R_sun/1000 em eV^{-1}; 6.957e+10 = raio solar em cm
*/
