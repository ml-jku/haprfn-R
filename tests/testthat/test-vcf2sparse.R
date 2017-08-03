context("vcf2sparse")

suppressPackageStartupMessages(library(fabia))
suppressPackageStartupMessages(library(hapFabia))

source("helper.R")

create_dir <- function(name) {
  directory <- paste0("test-vcf2sparse_", name, "_output/")
  path <- file.path(".", directory)
  dir.create(path, showWarning = FALSE)

  return(path)
}

test_that("Splitting matrices works", {
  base_filename = "chr22-1000"

  rfn_path <- create_dir("rfn")
  vcf2sparse(base_filename, intervalSize = 150, shiftSize = 100, outputPrefixPath = rfn_path)
  
  fabia_path <- create_dir("fabia")
  fabia_filename <- paste0(fabia_path, base_filename)
  vcftoFABIA(base_filename, outputFile = paste0(fabia_path, base_filename))
  split_sparse_matrix(fabia_filename, 150, shiftSize = 100, sparseMatrixPostfix = "_matG.txt")

  postfix <- "_0_150_matG"
  rfn_res <- readSparseSamples(paste0(rfn_path, base_filename, postfix, ".txt"))
  fabia_res <- readSamplesSpfabia(paste0(fabia_path, base_filename, postfix))

  expect_equal(rfn_res, fabia_res) 

  unlink(rfn_path)
  unlink(fabia_path)
})
