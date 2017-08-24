#' @param prefixPath The path to the input file. Default is the current directory.
#' @param intervalSize Size of the interval of each split.
#'   Default value = 10000.
#' @param shiftSize Size of the distance between the beginning of the
#'   intervals. Should be smaller than the interval size, otherwise
#'   features are skipped between intervals. Default value = 5000.
#' @param annotationPostfix The postfix for the annotation files.
#'   Default value = "_annot.txt".
#' @param individualsPostfix The postfix for the files containing
#'   the individuals. Default value = "_individuals.txt".
#' @param infoPostfix The postfix for the info file.
#'   Default value = "_info.txt".
