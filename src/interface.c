#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include "interface.h"


R_CallMethodDef callMethods[] = {
  {"readSparseSamples", (DL_FUNC) &readSparseSamples, 4},
  {"samplesPerFeature", (DL_FUNC) &samplesPerFeature, 4},
  {"vcf2sparse", (DL_FUNC) &vcf2sparse, 10},
  {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

int main() {
  return 1;
}
