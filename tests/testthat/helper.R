library(Matrix)
set.seed(123456)

writeMatrixFiles <- function(n, m, sp = 0.9, sd = 1, basePath = "./", baseName = "samplesPerFeature") {
  mat <- generateMatrix(n, m, sp, sd)
  
  fileNames <- paste0(basePath,baseName,"_", c("csr", "fabia", "expected"), "_", n, "_", m, "_", format(sp, digits = 2), "_", format(sd, digits = 2), ".txt")
  
  matrixToCSR(mat, fileNames[1])
  matrixToFabiaSparse(mat, fileNames[2])
  
  setNames(as.list(fileNames), c("csr", "fabia", "expected"))
}

removeMatrixFiles <- function(fileNames) {
  file.remove(fileNames$csr)
  file.remove(fileNames$fabia)
}

getFabiaFile <- function(fileNames) {
  substr(fileNames$fabia, 0, nchar(fileNames$fabia) - 4)
}

getCsrFile <- function(fileNames) {
  fileNames$csr
}

getExpectedFile <- function(fileNames) {
  paste0(substr(fileNames$expected, 0, nchar(fileNames$expected) - 4), ".RData")
}

generateMatrix <- function(n, m, sp, sd = 1) {
  set.seed(123456)
  
  mat <- replicate(m, rnorm(n, sd = sd))
  mat[sample(1:(n*m), round(n*m*sp))] = 0
  mat
}


matrixToCSR <- function(x, fname) {
  sp <- Matrix(t(x), sparse = TRUE)
  
  write(ncol(sp), fname, sep = "\n")
  write(nnzero(sp), fname, append = TRUE, sep = "\n")
  write(sp@p, fname, append = TRUE, sep = " ", ncolumns = length(sp@p))
  if (length(sp@i) > 0) {
    write(sp@i, fname, append = TRUE, sep = " ", ncolumns = length(sp@i))
    write(sp@x, fname, append = TRUE, sep = " ", ncolumns = length(sp@x))
  } else {
    write(x = "", fname, append = TRUE)
    write(x = "", fname, append = TRUE)
  }
}


matrixToFabiaSparse <- function(x, fname) {
  write(nrow(x), fname, sep = "\n")
  write(ncol(x), fname, append = TRUE, sep = "\n")
  
  for (i in 1:nrow(x)) {
    sp <- Matrix(x[i, ], sparse = TRUE)
    
    write(nnzero(sp), fname, append = TRUE, sep = "\n")
    if (length(sp@i) > 0) {
      write(sp@i, fname, append = TRUE, sep = " ", ncolumns = length(sp@i))
      write(sp@x, fname, append = TRUE, sep = " ", ncolumns = length(sp@x))  
    } else {
      write(x = "", fname, append = TRUE)
      write(x = "", fname, append = TRUE)
    }
  }
}


generateResultsFile <- function(n, m, sp = 0.9, sd = 1, basePath = "./") {
  fileNames <- writeMatrixFiles(n, m, sp, sd, basePath)
  
  require(fabia)
  
  expected <-fabia::samplesPerFeature(getFabiaFile(fileNames))
  
  save(expected, file = getExpectedFile(fileNames))
  
  removeMatrixFiles(fileNames)
}

generateResultsFileReadSamplesSpRfn <- function(n, m, sp = 0.9, sd = 1, basePath = "./") {
  fileNames <- writeMatrixFiles(n, m, sp, sd, basePath, "readSamplesSpRfn")
  
  require(fabia)
  
  expected <-fabia::readSamplesSpfabia(getFabiaFile(fileNames))
  
  save(expected, file = getExpectedFile(fileNames))
  
  removeMatrixFiles(fileNames)
}

#generateResultsFile(10, 10, 0.9, 1)
#generateResultsFile(1000, 1000, 0.9, 1)
#generateResultsFile(100, 200, 0.8, 40)
#generateResultsFile(100, 200, 0.1, 100)

generateResultsFileReadSamplesSpRfn(10, 10, 0.9, 1)
#generateResultsFile(1000, 1000, 0.9, 1)
#generateResultsFile(100, 200, 0.8, 40)
#generateResultsFile(100, 200, 0.1, 100)