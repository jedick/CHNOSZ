#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(h2o92)(void *, void *, void *, void *);

extern void F77_NAME(ideal2)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"h2o92", (DL_FUNC) &F77_NAME(h2o92), 4},
    {"ideal2", (DL_FUNC) &F77_NAME(ideal2), 8},
    {NULL, NULL, 0}
};

/* .C calls */
extern void count_letters(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"count_letters", (DL_FUNC) &count_letters, 2},
    {NULL, NULL, 0}
};

void R_init_CHNOSZ(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
/*
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
*/
    R_useDynamicSymbols(dll, FALSE);
}

