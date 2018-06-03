void vcf2sparse(SEXP file_nameS, SEXP prefix_pathS, SEXP interval_sizeS, SEXP shift_sizeS,
    SEXP annotateS, SEXP genotypesS, SEXP haplotypesS, SEXP missing_valuesS,
    SEXP annotation_postfixS, SEXP genotypes_postfixS, SEXP haplotypes_postfixS, SEXP info_postfixS,
    SEXP individuals_postfixS, SEXP output_fileS, SEXP output_prefixS);

SEXP samplesPerFeature(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS);

SEXP readSparseSamples(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS);

SEXP similarityMeasure(SEXP x, SEXP y, SEXP simv, SEXP minInter);
