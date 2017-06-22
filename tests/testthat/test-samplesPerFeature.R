context("samplesPerFeature")

source("helper.R")

testWithParameters <- function(n, m, sp, sd) {
  test_that(paste0("It works on some examples: nrow ", n, ", ncol ", m, 
                   ", sparsity ", sp, ", standard deviation ", sd), {
    fnames <- writeMatrixFiles(n, m, sp, sd, "./")
    load(getExpectedFile(fnames))
    
    expect_equal(samplesPerFeature(getCsrFile(fnames)), expected)

    removeMatrixFiles(fnames)
  })
}

testWithParameters(10, 10, 0.9, 1)
testWithParameters(1000, 1000, 0.9, 1)
testWithParameters(100, 200, 0.8, 40)
testWithParameters(100, 200, 0.1, 100)


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
