context("samplesPerFeature")

suppressPackageStartupMessages(library(fabia))

source("helper.R")

testWithParameters <- function(n, m, sp, sd) {
  test_that(paste0("It works on some examples: nrow ", n, ", ncol ", m, 
                   ", sparsity ", sp, ", standard deviation ", sd), {
    fnames <- writeMatrixFiles(n, m, sp, sd, "./")
    
    expect_equal(hapRFN::samplesPerFeature(getCsrFile(fnames)), 
                 fabia::samplesPerFeature(getFabiaFile(fnames)))

    removeMatrixFiles(fnames)
  })
}

testWithParameters(10, 10, 0.9, 1)
testWithParameters(1000, 1000, 0.9, 1)
testWithParameters(100, 200, 0.8, 40)
testWithParameters(100, 200, 0.1, 100)

fnames <- writeMatrixFiles(10, 10, 0.9, 1)

testWithSamples <- function(n, lowerB= 0 , upperB = 1000) {
  test_that(paste0("It works for selecting ", n, " random sample with lower bound",
                   lowerB, ", with upper bound ", upperB), {
    samples <- sample(1:10, n)
    expect_equal(hapRFN::samplesPerFeature(getCsrFile(fnames), samples, lowerB, upperB), 
                 fabia::samplesPerFeature(getFabiaFile(fnames), samples, lowerB, upperB)
    )
  })
}

testWithSamples(0)
testWithSamples(2, 0, 2)
testWithSamples(8, 0, 4)
testWithSamples(10, 0, 3)

removeMatrixFiles(fnames)

test_that("Works in extreme values", {
  # warn about nnz = 0 nrow = 0
  expect_warning(samplesPerFeature("zero.txt"))
  expect_warning(samplesPerFeature("null_matrix.txt"))
  expect_warning(samplesPerFeature("zero.txt", lowerB = 1000, 
                                    upperB = 1001))
})

test_that("Throws error", {
  expect_error(samplesPerFeature())
  expect_error(samplesPerFeature("nonexistent_file.txt"))
  expect_error(samplesPerFeature("malformed_sparse.txt"))
  expect_error(samplesPerFeature("malformed_sparse2.txt"))
})
