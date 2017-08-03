context("vcf2sparse")

suppressPackageStartupMessages(library(fabia))

source("helper.R")

test_that("Splitting matrices works", {
  vcf2sparse("chr22-1000", intervalSize = 150, shiftSize = 100, outputFile = "test_output")
})
