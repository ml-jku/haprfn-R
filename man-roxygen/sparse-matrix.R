#' @details
#'   The data matrix is directly scanned by the C code and must be in sparse
#'   matrix format.
#'
#' Sparse matrix format where rows are as follows:
#' \itemize{
#'   \item number of non-zero elements (\eqn{nnz})
#'   \item number of rows (\eqn{nrow})
#'   \item row pointers, a list of integers of length \eqn{nrow + 1},
#'     where the ith element is the number of non-zero elements 
#'     until the \eqn{i - 1}th row
#'   \item column indices, a list of integers of length \eqn{nnz}, 
#'     where the \eqn{i}th element is the (zero based) column index 
#'     of the \eqn{i}th non-zero element
#'   \item values, the list of non-zero elements
#' }
#'
#' @param X The file name of the sparse matrix, in sparse format.
#' @template sparse-matrix-filter
