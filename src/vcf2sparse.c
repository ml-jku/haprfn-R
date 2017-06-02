#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include "vcf.h"

void vcf2sparse(SEXP file_nameS, SEXP prefix_pathS, SEXP no_snvsS, SEXP output_fileS) {
    const char* file_name = CHAR(STRING_ELT(file_nameS, 0));
    //const char *arg2=CHAR(STRING_ELT(arg2S,0));
    //const char *arg3=CHAR(STRING_ELT(arg3S,0));
    //const char *arg4=CHAR(STRING_ELT(arg4S,0));
	htsFile* file = hts_open(file_name, "r");
    bcf1_t* bcf;// = bcf_init();
    bcf_hdr_t* hdr;// = bcf_hdr_read(file);

    int ret = 0;
    //while ((ret = bcf_read(file, hdr, bcf)) >= 0) {
        // iterate
    //}

    //bcf_hdr_destroy(hdr);
    //bcf_destroy(bcf);
    //hts_close(file);
}
