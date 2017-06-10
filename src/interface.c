#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include "interface.h"


R_CallMethodDef callMethods[] = {
  {"readSamplesSpRfn", (DL_FUNC) &readSamplesSpRfn, 4},
  {"samplesPerFeature", (DL_FUNC) &samplesPerFeature, 4},
  {"vcf2sparse", (DL_FUNC) &vcf2sparse, 6},
  {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

int main() {
  return 1;
}
