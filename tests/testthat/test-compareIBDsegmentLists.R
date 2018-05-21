context("compareIBDsegmentLists")

suppressPackageStartupMessages(library(hapFabia))

test_that("Compare IBD segments gives the same result", {
  # Example taken from hapFabia::compareIBDsegmentLists documentation
  data(hapRes)
  IBDsegmentList1 <- hapRes$IBDsegmentList1
  IBDsegmentList2 <- hapRes$IBDsegmentList2
  comp1 <- 
    hapFabia::compareIBDsegmentLists(IBDsegmentList1, IBDsegmentList2, 
                                     simv = "minD", pTagSNVs = NULL, 
                                     pIndivid = NULL, minTagSNVs = 6,
                                     minIndivid = 2)
  comp2 <- 
    hapRFN::compareIBDsegmentLists(IBDsegmentList1, IBDsegmentList2,
                           simv = "minD", pTagSNVs = NULL, pIndivid = NULL,
                           minTagSNVs = 6, minIndivid = 2)
  expect_equal(comp1, comp2)
})
