#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void deleteRNGstates(void);
extern void draw_lambda_R(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void draw_omega_R(void *, void *, void *, void *, void *);
extern void draw_z_R(void *, void *, void *, void *, void *, void *, void *);
extern void newRNGstates(void);
extern void rinvgauss_R(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"deleteRNGstates", (DL_FUNC) &deleteRNGstates, 0},
    {"draw_lambda_R",   (DL_FUNC) &draw_lambda_R,   9},
    {"draw_omega_R",    (DL_FUNC) &draw_omega_R,    5},
    {"draw_z_R",        (DL_FUNC) &draw_z_R,        7},
    {"newRNGstates",    (DL_FUNC) &newRNGstates,    0},
    {"rinvgauss_R",     (DL_FUNC) &rinvgauss_R,     4},
    {NULL, NULL, 0}
};

void R_init_reglogit(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
