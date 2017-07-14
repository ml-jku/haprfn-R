void vcf2sparse(SEXP file_nameS, SEXP prefix_pathS, SEXP interval_sizeS, SEXP shift_sizeS, SEXP annotateS, SEXP genotypesS, SEXP haplotypesS, SEXP missing_valuesS, SEXP output_fileS, SEXP output_prefixS);

SEXP samplesPerFeature(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS);

SEXP readSparseSamples(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS);
