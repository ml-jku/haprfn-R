#' @param sparseMatrixPostfix The posftix for the input matrices.
#'   Default = "_mat.txt".
#' @param lowerBP Lower bound for filtering the inputs columns, minimal MAF 
#'   (however more than one occurence to remove private SNVs).
#'   Default = 0.0.
#' @param upperBP Upper bound for filtering the inputs columns, minimal MAF.
#'   Default = 0.05.
#' @param p Number of hidden units (biclusters). Default = 50.
#' @param l1 l1 weight decay. Default = 0.0. 
#' @param cyc Number of iterations. Default = 100.
#' @param etaW Learning rate of the W parameter
#' @param etaP Learning rate of the Psi parameter (It's probably save to set this to the same value as etaW)
#' @param minP Minimal value for Psi. Should be in 1e-5 - 1e-1
#' @param write_file Results are written to files (L in sparse format).
#'   Default = 0 (not written).
#' @param IBDsegmentLength IBD segment length in kbp. Default = 50.
#' @param Lt Percentage of largest Ls to consider for IBD segment extraction.
#'   Default = 0.1.
#' @param Zt Percentage of largest Zs to consider for IBD segment extraction.
#'   Default = 0.2.
#' @param thresCount p-value of random histogram hit. Default = 1e-5.
#' @param mintagSNVsFactor Percentage of segments overlap in IBD segments.
#'   Use 1/2 for large to 3/4 for small intervals. Default = 3/4.
#' @param pMAF Averaged and corrected minor allele frequency. Default = 0.03.
#' @param haplotypes TRUE if matrix contains haplotype data. Haplotypes are
#'   phased genotypes, so two chromosomes are expected per individual.
#'   Default = FALSE.
#' @param cut Default = 0.8.
#' @param procMinIndivids Default = 0.1.
#' @param thresPrune Default = 0.001.
#' @param simv Default = "minD".
#' @param thresA Default = NULL.
#' @param minTagSNVs Default = NULL.
#' @param minIndivid Default = 2.
#' @param avSNVsDist Default = 100.
#' @param SNVclusterLength Default = 100.
#' @param useGpu TRUE to use GPU. Use this parameter together with gpuId.
#'   Default = FALSE.
#' @param gpuId the ID of the GPU. This is used, when useGpu is TRUE.
#'   Default = -1 (invalid GPU ID).
#' @param seed random number generator seed for RFN

