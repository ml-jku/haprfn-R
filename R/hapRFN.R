#' @import Rhtslib 
#' @import hapFabia
#' @import RFN
#' @useDynLib hapRFN
#'
#' @title Samples with non-zero feature.
#' 
#' @description
#'   \code{samplesPerFeature} returns the samples for each feature, where
#'   this is non-zero
#'
#' @details
#'   Supplies the samples for which a feature is not zero.
#'
#' @template descr-sparse-matrix
#'
#' @param X The name of the file, containing a matrix in sparse format.
#' @template param-samples
#' @template param-matrix-filter
#'
#' @return A list with elements
#'   \item{sL}{List with one element per feature: each element is a vector
#'     of samples where the feature is not zero.}
#'   \item{nsL}{Vector of feature length containing number of samples having
#'     a non-zero feature value.}
#'
#' @export
#'
samplesPerFeature <- 
function(X, samples = 0, lowerB = 0, upperB = 1000) {
  if (missing(X)) {
    stop("Data file name missing.")
  }
  
  X <- as.character(X)
  samples <- as.integer(sort.int(as.integer(unique(samples))))
  lowerB <- as.double(lowerB)
  upperB <- as.double(upperB)
  
  if (lowerB > upperB) {
    warning(paste0("Lower bound (", lowerB, 
                   ") is greater than the upper bound (", upperB,")."))
  }
  
  .Call("samplesPerFeature", X, samples, lowerB, upperB, PACKAGE = "hapRFN")
}

#' @title Read sparse matrix samples.
#' 
#' @description \code{readSparseSamples} Reads sparse matrix samples.
#'
#' @details
#'   Reads the sparse matrix samples filtered by \code{samples}, \code{lowerB}
#'.  and \code{upperB}.
#'
#' @template descr-sparse-matrix
#'
#' @param X The name of the file, containing a matrix in sparse format.
#' @template param-samples
#' @template param-matrix-filter
#'
#' @return Data matrix of given samples
#'
#' @export
#'
readSparseSamples <-
function(X, samples = 0, lowerB = 0, upperB = 1000) {
  if (missing(X)) {
    stop("Data file name missing. Stopped.")
  }
  
  X <- as.character(X)
  samples <- as.integer(sort.int(as.integer(unique(samples))))
  lowerB <- as.double(lowerB)
  upperB <- as.double(upperB)
  
  .Call("readSparseSamples", X, samples, lowerB, upperB, PACKAGE = "hapRFN")
}

#' @title Transforming VCF files.
#'
#' @description
#'   \code{vcf2sparse} converts and splits VCF files into sparse matrices.
#'
#' @details
#'   Reads a VCF file, converts it into a sparse matrix format, splits it into
#'   intervals and writes the result into files. It generates multiple files.
#'
#'   The code is implemented in C.
#'
#' @param fileName The name of the input VCF file without the extension. 
#'   The extension must be one of .vcf or .vcf.gz, 
#'   and the format must match the extension. Required.
#' @template param-in-file
#' @template param-interval
#' @param annotation TRUE to generate the annotation file for each interval.
#'   Annotation file names end with \code{_annot.txt}. Default value = TRUE.
#' @param haplotypes TRUE to generate the haplotype matrix for each interval.
#'   Rows are features. For each sample there are two columns, for each
#'   haplotype, having a value of \eqn{1} if the variant is present, zero
#'   otherwise. Default value = FALSE.
#' @param genotypes TRUE to generate genotype matrix for each interval.
#'   The genotype matrix is a combined version of the haplotype matrix.
#'   Each row is a feature and each column a sample, having a value of
#'   \eqn{2} if variant is present in both haplotypes, \eqn{1} if the 
#'   variant is present in one haplotype, \eqn{0} otherwise. Default
#'   value = TRUE.
#' @param missingValues Flag specifying how to deal with missing values.
#'   \describe{
#'     \item{0}{look for frequency defined in the INFO tag as AF or estimate}
#'     \item{1}{assume missing value as major allele}
#'     \item{2}{assume missing value as minor allele}
#'     \item{3}{look for frequency defined in the INFO tag as AF or abort}
#'     \item{4}{estimate based on available information}
#'     \item{5}{abort on missing value}
#'   }
#'   Default value = 0.
#' @param genotypesPostfix The postfix for the genotypes files. Default
#'   value = "_matG.txt".
#' @param haplotypesPostfix The postfix for the haplotypes files. Default
#'   value = "_math.txt".
#' @param outputFile The base name of the output file. The file postfix and
#'   extension will be appended to this file. Default is the input file name.
#' @param outputPrefixPath The path to the output file.
#'   Default is the \code{prefixPath}.
#'
#' @return None
#'
#' @export
vcf2sparse <-
function(fileName, prefixPath = NULL, intervalSize = 10000, shiftSize = 5000,
         annotation = TRUE, genotypes = TRUE, haplotypes = FALSE,
         missingValues = 0, annotationPostfix = "_annot.txt",
         genotypesPostfix = "_matG.txt", haplotypesPostfix = "_matH.txt",
         infoPostfix = "_info.txt", individualsPostfix = "_individuals.txt",
         outputFile = NULL, outputPrefixPath = NULL) {
  .Call("vcf2sparse", fileName, prefixPath, as.integer(intervalSize), 
        as.integer(shiftSize), annotation, genotypes, haplotypes, 
        as.integer(missingValues), annotationPostfix, genotypesPostfix,
        haplotypesPostfix, infoPostfix, individualsPostfix, outputFile,
        outputPrefixPath, PACKAGE = "hapRFN")
  invisible()
}

#' @title Iterate intervals.
#'
#' @description Iterate over intervals and run hapRFN.
#'
#' @param startRun The index of the first interval.
#' @param endRun The index of the last interval. Required.
#' @param annotationFile The name of the annotation file for the individuals.
#' @param fileName The base name of the input file. Required.
#' @template param-in-file
#' @template param-interval
#' @template param-samples
#' @template param-rfn
#' @param saveAsCsv Save merged IBD segment list as CSV file
#'   for every interval. Default = FALSE.
#'
#' @return None
#'
#' @export
iterateIntervals <-
function(startRun = 1, endRun, shiftSize = 5000, intervalSize = 10000, 
         annotationFile = NULL, fileName, prefixPath = "",
         sparseMatrixPostfix = "_mat.txt", annotationPostfix = "_annot.txt",
         individualsPostfix = "_individuals.txt", infoPostfix = "_info.txt",
         samples = 0, l1 = 0.0, lowerBP = 0, upperBP = 0.05, p = 50, cyc = 100,
         etaW = 0.1, etaP = 0.1, minP = 0.01, dropout = 0.0,
         noise_type = "saltpepper", write_file = 0, IBDsegmentLength = 50,
         thresCount = 1e-05, mintagSNVsFactor = 3/4, pMAF = 0.03,
         haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1,
         thresPrune = 0.001, simv = "minD", thresA = NULL, minTagSNVs = NULL,
         minIndivid = 2, avSNVsDist = 100, SNVclusterLength = 100,
         saveAsCsv = FALSE, useGpu = TRUE, gpuId = 0, seed = seed) {
  info <- .readInfo(prefixPath, fileName, infoPostfix)

  nsamples <- info$nsamples
  snvs <- info$nsnps
  
  save(nsamples, snvs, file = paste0(fileName, "_All", ".Rda"))

  for (posAll in startRun:endRun) {
    start <- (posAll - 1) * shiftSize
    end <- start + intervalSize

    if (end > snvs) {
      end <- snvs
    }
        
    labels <- .readLabels(prefixPath, fileName, individualsPostfix, 
                          annotationFile, haplotypes, nsamples)
    pRange <- .createRangeString(start, end)
    
    resHapRFN <- hapRFN(fileName = fileName, prefixPath = prefixPath, 
                        sparseMatrixPostfix = sparseMatrixPostfix,
                        annotationPostfix = annotationPostfix,
                        individualsPostfix = individualsPostfix,
                        labelsA = labels, pRange = pRange, samples = samples,
                        lowerBP = lowerBP, upperBP = upperBP, p = p, l1 = l1,
                        cyc = cyc, etaW = etaW, etaP = etaP, minP = minP,
                        dropout = dropout, noise_type = noise_type,
                        write_file = write_file,
                        IBDsegmentLength = IBDsegmentLength,
                        thresCount = thresCount,
                        mintagSNVsFactor = mintagSNVsFactor, pMAF = pMAF,
                        haplotypes = haplotypes, cut = cut,
                        procMinIndivids = procMinIndivids,
                        thresPrune = thresPrune, simv = simv, thresA = thresA,
                        minTagSNVs = minTagSNVs, minIndivid = minIndivid,
                        avSNVsDist = avSNVsDist,
                        SNVclusterLength = SNVclusterLength, useGpu = useGpu,
                        gpuId = gpuId, seed = seed)
    
    if (saveAsCsv) {
      IBDsegmentList2excel(resHapRFN$mergedIBDsegmentList, 
                           paste0(fileName, pRange, ".csv"))  
    }
    
    annot <- c()
    save(resHapRFN, annot, file = paste0(fileName, pRange, "_resAnno", ".Rda"))
  }
}

#' @title Identifying IBD segments.
#'
#' @description Identifying short IBD segments using RFNs.
#'
#' @details 
#'   Identifying short identity by descent (IBD) segments using Rectified
#'   Factor Networks (RFNs).
#'
#' @param fileName The base name of the input file.
#' @template param-in-file
#' @param labelsA Individual names as matrix individuals x 4.
#' @param pRange For intervals indicates its range.
#' @template param-samples
#' @template param-rfn
#'
#' @template ref-rfn
#'
#' @return None
#'
#' @export
#' @importFrom stats pbinom cutree
#' @importFrom utils read.table
#' @importFrom RFN train_rfn
#' @importFrom methods new
#'
hapRFN <-
function(fileName, prefixPath = "", sparseMatrixPostfix = "_mat.txt",
         annotationPostfix = "_annot.txt",
         individualsPostfix = "_individuals.txt", infoPostfix = "_info.txt",
         labelsA = NULL, pRange = "", samples = 0, lowerBP = 0, upperBP = 0.05,
         p = 50, etaW = 0.1, etaP = 0.1, minP = 0.01, dropout = 0.0,
         noise_type = "saltpepper", l1 = 0.0, cyc = 100, write_file = 0,
         IBDsegmentLength = 50, thresCount = 1e-05, mintagSNVsFactor = 3/4,
         pMAF = 0.03, haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1,
         thresPrune = 0.001, simv = "minD", thresA = NULL, minTagSNVs = NULL,
         minIndivid = 2, avSNVsDist = 100, SNVclusterLength = 100,
         useGpu = FALSE, gpuId = -1, seed = 0) {
  message("                      ")
  message("                      ")
  message("Running hapRFN with:")
  message("   Prefix string for file name of data files --------- : ", fileName)
  message("   Path of data files ---------------------------------: ", prefixPath)
  if (haplotypes) {
    message("   Data consists of phased genotypes (haplotypes) -----")
  } else {
    message("   Data consists of unphased genotypes ----------------")
  }
  message("   Postfix string for file in sparse matrix format ----: ", sparseMatrixPostfix)
  message("   Postfix string for file containing individual names : ", individualsPostfix)
  if (is.null(labelsA)) {
    message("   Individuals annotation is not supplied -------------")
  } else {
    message("   Individuals annotation is supplied -----------------")
  }
  message("   String indicating the interval that is analyzed ----: ", pRange)
  if (length(samples) > 1) {
    message("   Number of individuals included into the analysis ---: ", length(samples))
  } else {
    message("   All individuals are included into the analysis -----: 0 = all individuals")
  }
  message("   Lower bound on MAF but more than one occurence -----: ", lowerBP)
  message("   Upper bound on MAF ---------------------------------: ", upperBP)
  message("   Number of biclusters -------------------------------: ", p)
  message("   Number of cycles -----------------------------------: ", cyc)
  message("   Results written to files ---------------------------: ", write_file)
  message("   Learning rate of the W parameter -------------------: ", etaW)
  message("   Learning rate of the Psi parameter -----------------: ", etaP)
  message("   Minimal value for Psi ------------------------------: ", minP)
  if (IBDsegmentLength > 0) {
    message("   IBD segment length in kbp --------------------------: ", IBDsegmentLength)
    message("   Average distance between SNVs in bases -------------: ", avSNVsDist)
  } else {
    if (SNVclusterLength > 0) {
      message("   Number of SNVs in histogram bin --------------------: ", 
        SNVclusterLength)
    } else {
      message("   Number of SNVs in histogram bin --------------------: 100")
    }
  }
  message("   p-value threshold for random histogram counts ------: ", thresCount)
  message("   Min. % of segments overlap in IBD segments ---------: ", mintagSNVsFactor)
  message("   Averaged and corrected minor allele frequency ------: ", pMAF)
  message("   Cutoff for merging IBD segments --------------------: ", cut)
  message("   % of cluster individuals a tagSNV must tag ---------: ", procMinIndivids)
  message("   Threshold for pruning border tagSNVs ---------------: ", thresPrune)
  message("   Similarity measure for merging clusters ------------: ", simv)
  message("   Minimum matching tagSNVs for cluster similarity ----: ", minTagSNVs)
  message("   Minimum matching individuals for cluster similarity : ", minIndivid)
  message("   Seed -----------------------------------------------: ", seed)
  message("                      ")
  message("                      ")
  
  # Compute internal parameters
  
  # compute histogram length
  if (IBDsegmentLength > 0) {
    if (avSNVsDist > 0) {
      inteA <- IBDsegmentLength * (1000/avSNVsDist)
    } else {
      inteA <- IBDsegmentLength * 10
    }
  } else {
    if (SNVclusterLength > 0) {
      inteA <- SNVclusterLength
    } else {
      inteA <- 100
    }
  }

  info <- .readInfo(prefixPath, fileName, infoPostfix)
  snvs <- info$nsnps
  individualsN <- info$nsamples
  
  if (length(samples) > 1) {
    individualsN <- length(samples)
  } else {
    individuals <- 0
  }
  
  upperBindivid <- upperBP * individualsN  # remove common SNVs
  lowerBindivid <- max(1.5, lowerBP * individualsN)  # remove private SNVs
  
  if (is.null(thresA)) {
    kk <- 1
    while ((snvs / inteA) * choose(individualsN, 2) * 
           (1 - pbinom(kk, inteA, pMAF * pMAF)) > thresCount) {
      kk <- kk + 1
    }
    thresA <- kk
  }

  if (is.null(minTagSNVs)) {
    minTagSNVs <- round(mintagSNVsFactor * thresA)
  }
  
  # End Compute internal parameters
  
  matrixFileName <- paste0(prefixPath, fileName, pRange, sparseMatrixPostfix)
  X <- .readSparseMatrix(matrixFileName)

  message("start RFN")

  # RFN call
  l <- ncol(X)
  n <- nrow(X)
  
  rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "SNV")
  colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
  rowna <- rownames(X)
  colna <- colnames(X)
  
  iin <- 1/l
  
  cyc <- as.integer(cyc)
  
  p <- as.integer(p)

  com <- which(table(X@i)/X@Dim[2] > upperBP)
  if (length(com) > 0) {
    X[com + 1, ] <- 0
  }
  
  rfn_res <- train_rfn(X = X, n_hidden = p, n_iter = cyc, etaW = etaW, 
                       etaP = etaP, minP = minP, dropout_rate = dropout,
                       noise_type = noise_type, l1_weightdecay = l1, 
                       seed = seed, use_gpu = useGpu, gpu_id = gpuId)
  
  myL <- rfn_res$W
  myPsi <- as.vector(rfn_res$P)
  myZ <- rfn_res$H
  mylapla <- as.matrix(1)
  
  vz <- iin * apply(myZ, 1, function(x) sum(x ^ 2))
  vz <- sqrt(vz + 1e-10)
  ivz <- 1 / vz
  if (length(ivz) == 1) {
    nZ <- ivz * myZ
    noL <- vz * myL
  } else {
    nZ <- ivz * myZ
    noL <- t(vz * t(myL))
  }
  
  Lz <- noL %*% nZ
  rownames(noL) <- rowna
  colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
  clnames <- colnames(noL)
  rownames(nZ) <- clnames
  colnames(nZ) <- colna
  
  ini <- as.matrix(1)
  avini <- as.vector(1)
  xavini <- as.vector(1)
  
  res <- new("Factorization", parameters = list("rfn", cyc, p), n = n, p1 = p, 
             p2 = p, l = l, center = as.vector(1), scaleData = as.vector(1),
             X = as.matrix(1), L = noL, Z = nZ, M = as.matrix(1),
             LZ = as.matrix(1), U = as.matrix(1), avini = avini,
             xavini = xavini, ini = ini, Psi = myPsi, lapla = mylapla)
  
  # Load individuals to Ls of interest: load minor alleles of the Ls
  
  sparseMatrixFilename <- paste0(prefixPath, fileName, pRange, sparseMatrixPostfix)
  sPF <- hapRFN::samplesPerFeature(X = sparseMatrixFilename, samples = samples,
                                   lowerB = lowerBindivid, upperB = upperBindivid)

  if (nchar(annotationPostfix) > 0) {
    # annot[[1]] <- chromosome
    # annot[[2]] <- phys. position
    # annot[[3]] <- snvNames
    # annot[[4]] <- snvMajor 
    # annot[[5]] <- snvMinor 
    # annot[[6]] <- quality 
    # annot[[7]] <- pass 
    # annot[[8]] <- info of vcf file 
    # annot[[9]] <- fields in vcf file
    # annot[[10]] <- frequency 
    # annot[[11]] <- 1 = changed if major allele is 
    # actually minor allele otherwise 0
    
    annot <- read.table(paste0(prefixPath, fileName, pRange, annotationPostfix), 
                        header = FALSE, sep = "\t", quote = "", as.is = TRUE)
    
    for (i in 1:length(annot)) {
      annot[[i]] <- gsub(",", ";", annot[[i]])
    }
    
    for (i in 1:length(annot)) {
      annot[[i]] <- gsub("TRUE", "T", annot[[i]])
    }
    
    annot[[2]] <- as.numeric(annot[[2]])  # physical position
    annot[[10]] <- as.numeric(annot[[10]])  # SNV frequency
    annot[[11]] <- as.numeric(annot[[11]])  # changed
  }
  
  if (is.null(labelsA)) {
    labelsAA <- read.table(paste0(prefixPath, fileName, individualsPostfix), 
                           header = FALSE, sep = " ", quote = "", as.is = TRUE)
    if (haplotypes) {
      lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
    } else {
      lA <- as.vector(labelsAA[, 2])
    }
    indiA <- cbind(as.character(lA), as.character(lA),
                   as.character(lA), as.character(lA))
  } else {
    indiA <- labelsA
  }
 
 # first haplotype extraction with offset 0
  
  off1 <- 0
  IBDsegmentList1 <- extractIBDsegments(res = res, sPF = sPF, annot = annot, chrom = "",
                                        labelsA = indiA, ps = 0.0, psZ = 0.0, inteA = inteA,
                                        thresA = thresA, mintagSNVs = minTagSNVs, off = off1,
                                        procMinIndivids = procMinIndivids, thresPrune = thresPrune)
  
  # merge IBD segment lists
  
  if (lengthList(IBDsegmentList1) > 1) {
    comp <- compareIBDsegmentLists(IBDsegmentList1 = IBDsegmentList1, IBDsegmentList2 = NULL, 
                                   simv = simv, pTagSNVs = NULL, pIndivid = NULL,
                                   minTagSNVs = minTagSNVs, minIndivid = minIndivid)
    
    if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp, h = cut)
      mergedIBDsegmentList1 <- mergeIBDsegmentLists(IBDsegmentList1 = IBDsegmentList1, 
                                                    IBDsegmentList2 = NULL,
                                                    clustIBDsegmentList = clustIBDsegmentList)
    } else {
      mergedIBDsegmentList1 <- IBDsegmentList1
    }
  } else {
    mergedIBDsegmentList1 <- IBDsegmentList1
  }
  
  
  # second IBD extraction with offset half of the interval length
  
  off2 <- inteA %/% 2
  
  IBDsegmentList2 <- extractIBDsegments(res = res, sPF = sPF, annot = annot, chrom = "", 
                                        labelsA = indiA, ps = ps, psZ = psZ, inteA = inteA,
                                        thresA = thresA, mintagSNVs = minTagSNVs, off = off2,
                                        procMinIndivids = procMinIndivids, thresPrune = thresPrune)
  
  # merge IBD segments
  
  if (lengthList(IBDsegmentList2) > 1) {
    comp <- compareIBDsegmentLists(IBDsegmentList1 = IBDsegmentList2, IBDsegmentList2 = NULL,
                                   simv = simv, pTagSNVs = NULL, pIndivid = NULL,
                                   minTagSNVs = minTagSNVs, minIndivid = minIndivid)
    
    if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp, h = cut)
      mergedIBDsegmentList2 <- mergeIBDsegmentLists(IBDsegmentList1 = IBDsegmentList2, 
                                                    IBDsegmentList2 = NULL,
                                                    clustIBDsegmentList = clustIBDsegmentList)
    } else {
      mergedIBDsegmentList2 <- IBDsegmentList2
    }
  } else {
    mergedIBDsegmentList2 <- IBDsegmentList2
  }
  
  # merge IBD segments of both extractions
  if (lengthList(mergedIBDsegmentList1) > 0) {
    if (lengthList(mergedIBDsegmentList2) > 0) {
      comp12 <- compareIBDsegmentLists(IBDsegmentList1 = mergedIBDsegmentList1, 
        IBDsegmentList2 = mergedIBDsegmentList2, simv = simv, pTagSNVs = NULL, 
        pIndivid = NULL, minTagSNVs = minTagSNVs, minIndivid = minIndivid)
      
      if (!is.null(comp12)) {
        clustIBDsegmentList <- cutree(comp12, h = cut)
        mergedIBDsegmentList <- mergeIBDsegmentLists(IBDsegmentList1 = mergedIBDsegmentList1, 
                                                     IBDsegmentList2 = mergedIBDsegmentList2,
                                                     clustIBDsegmentList = clustIBDsegmentList)
      } else {
        mergedIBDsegmentList <- mergedIBDsegmentList1
      }
    } else {
      mergedIBDsegmentList <- mergedIBDsegmentList1
    }
  } else {
    if (lengthList(mergedIBDsegmentList2) > 0) {
      mergedIBDsegmentList <- mergedIBDsegmentList2
    } else {
      mergedIBDsegmentList <- mergedIBDsegmentList1
    }
  }
  
  mergedIBDsegmentList1 <- setStatistics(mergedIBDsegmentList1)
  mergedIBDsegmentList2 <- setStatistics(mergedIBDsegmentList2)
  mergedIBDsegmentList <- setStatistics(mergedIBDsegmentList)
  
  
  return(list(mergedIBDsegmentList = mergedIBDsegmentList, res = res, sPF = sPF, annot = annot,
              IBDsegmentList1 = IBDsegmentList1, IBDsegmentList2 = IBDsegmentList2,
              mergedIBDsegmentList1 = mergedIBDsegmentList1,
              mergedIBDsegmentList2 = mergedIBDsegmentList2))
}

#' @title Identify duplicates.
#'
#' @description
#'   \code{identifyDuplicates} identifies IBD segments that are similar to
#'   each other. This function is used in combination with \code{vcf2sparse}
#'   which splits a VCF file into overlapping intervals. These intervals are
#'   then analyzed by \code{analyzeIntervals}, and the checked for duplicates
#'   by this function.
#'
#'   Results are written to the file "dups.Rda".
#'
#' @param fileName The name of the Rda result file.
#' @param startRun The index of the first interval.
#' @param endRun The index of the last interval. Required.
#' @template param-interval
#'
#' @return None
#'
#' @export
identifyDuplicates <-
function(fileName, startRun = 1, endRun, shiftSize = 5000, 
         intervalSize = 10000) {
  labelsA <- c()
  snvs <- c()
  resHapRFN <- c()

  # loads nsamples and snvs
  load(file = paste(fileName, "_All", ".Rda",sep=""))
  
  avIBDsegmentLength <- list()
  avIBDsegmentPos <- list()
  count <- 0
  
  for (posAll in startRun:endRun) {
    start <- (posAll-1)*shiftSize
    end <- start + intervalSize
    
    if (end > snvs) {
      end <- snvs
    }
    
    pRange <- .createRangeString(start, end)
    
    load(file = paste0(fileName, pRange, "_resAnno.Rda"))
    
    mergedIBDsegmentList <- resHapRFN$mergedIBDsegmentList
    
    noIBDsegments <- lengthList(mergedIBDsegmentList)
    
    if (noIBDsegments > 0) {
      count <- count + noIBDsegments
      
      avIBDsegmentPos[[posAll]] <- sapply(IBDsegments(mergedIBDsegmentList),
                                          function(x) {IBDsegmentPos(x)}, simplify=FALSE)
      avIBDsegmentLength[[posAll]] <- sapply(IBDsegments(mergedIBDsegmentList),
                                             function(x) {IBDsegmentLength(x)}, simplify=FALSE)
    }
  }

  IBDsegmentPos <- unlist(avIBDsegmentPos)
  IBDsegmentLength <- unlist(avIBDsegmentLength)
  
  IBDsegmentSim <- cbind(IBDsegmentPos, IBDsegmentLength)
  
  dups <- duplicated(IBDsegmentSim)

  un <- which(dups == FALSE)
  
  #### enumerate counts ####
  allCount <- 0
  allCount1 <- 0
  
  resD <- list()
  resDA <- list()
  
  for (posAll in startRun:endRun) {
    start <- (posAll-1)*shiftSize
    end <- start + intervalSize
    
    if (end > snvs) {
      end <- snvs
    }
    
    pRange <- .createRangeString(start, end)
    
    load(file = paste(fileName, pRange, "_resAnno.Rda", sep = ""))
    
    mergedIBDsegmentList <- resHapRFN$mergedIBDsegmentList
    
    noIBDsegments <- lengthList(mergedIBDsegmentList)
    
    if (noIBDsegments > 0) {
      for (IBDsegmentC in 1:noIBDsegments) {
        allCount <- allCount + 1
        
        resD1 <- c(allCount,IBDsegmentC,posAll)
        resDA[[allCount]] <- resD1
        
        if (!dups[allCount]) {
          allCount1 <- allCount1 + 1
          
          resD1 <- c(allCount1,allCount,IBDsegmentC,posAll)
          resD[[allCount1]] <- resD1
        }
      }
    }
  }
  
  if (length(resD) > 0) {
    rr <- unlist(resD)
    l <- 4
    
    countsA1 <- matrix(rr, nrow = allCount1, ncol = l, byrow = TRUE)
    
    colnames(countsA1) <- c("allCount1", "allCount", "IBDsegmentC", "posAll")
    
    rr <- unlist(resDA)
    l <- 3
    
    countsA2 <- matrix(rr, nrow = allCount, ncol = l, byrow = TRUE)
    
    colnames(countsA2) <- c("allCount","IBDsegmentC","posAll")
  } else {
    dups <- FALSE
    un <- 0
    countsA1 <- c(0, 0, 0, 0)
    dim(countsA1) <- c(1, 4)
    colnames(countsA1) <- c("allCount1", "allCount", "IBDsegmentC", "posAll")
    countsA2 <- c(0, 0, 0)
    dim(countsA2) <- c(1, 3)
    colnames(countsA2) <- c("allCount", "IBDsegmentC", "posAll")
  }
    
  save(dups, un, countsA1, countsA2, file = paste("dups.Rda", sep = ""))
}

#' @title Analyze IBD segments.
#'
#' @description
#'   \code{analyzeIBDsegments} loops over all intervals and provides
#'   statistics about the IBD segments.
#'
#' @details
#'   The function loops over all intervals and provides statistics about the
#'   IBD segments. The loop goes over the intervals that have been provided by
#'   \code{iterateIntervals}. Duplicates are ignore at this analysis and must
#'   be identified in a preceding step via \code{identifyDuplicates}.
#'   Other statistics and annotations can be computed if the code is changed
#'   accordingly.
#'
#'   Results are written to the file "analyzeResults.Rda".
#'
#' @param fileName The name of the Rda result file.
#' @param startRun The index of the first interval.
#' @param endRun The index of the last interval. Required.
#' @template param-interval
#'
#' @return A list with elements
#'   \item{startRun}{The parameter startRun.}
#'   \item{endRun}{The parameter endRun.}
#'   \item{noIBDsegments}{The number of IBD segments.}
#'   \item{avIBDsegmentPos}{Vector of the genomic locations of the IBD segments.}
#'   \item{avIBDsegmentLengthSNV}{Vector of lengths in SNVs of the IBD segments.}
#'   \item{avIBDsegmentLength}{Vector of lengths in bp of the IBD segments.}
#'   \item{avnoIndivid}{Vector of number of samples belonging to the IBD segments.}
#'   \item{avnoTagSNVs}{Vector of number of tag SNVs marking the IBD segments.}
#'   \item{avnoFreq}{Vector of frequencies of tagSNVs in the whole dataset.}
#'   \item{avnoGroupFreq}{Vector of frequencies of tagSNVs in the population
#'     that is considered.}
#'   \item{avnotagSNVChange}{Vector of flags indicating a switch between minor
#'     and major alleles of tagSNVs (1 means switched, 0 not switched).}
#'   \item{avnotagSNVsPerIndividual}{Vector of number of tagSNVs per individual.}
#'   \item{avnoindividualPerTagSNV}{Vector of number of individuals that posses
#'     the minor allele per tagSNV.}
#'   \item{avIBDsegmentPosS}{Summary of avIBDsegmentPos.}
#'   \item{avIBDsegmentLengthSNVS}{Summary of avIBDsegmentLengthSNV.}
#'   \item{avIBDsegmentLengthS}{Summary of avIBDsegmentLength.}
#'   \item{avnoIndividS}{Summary of avnoIndivid.}
#'   \item{avnoTagSNVsS}{Summary of avnoTagSNVs.}
#'   \item{avnoFreqS}{Summary of avnoFreq.}
#'   \item{avnoGroupFreqS}{Summary of avnoGroupFreq.}
#'   \item{avnotagSNVChangeS}{Summary of avnotagSNVChange.}
#'   \item{avnotagSNVsPerIndividualS}{Summary of avnotagSNVsPerIndividual.}
#'   \item{avnoindividualPerTagSNVS}{Summary of avnoindividualPerTagSNV.}
#'
#' @seealso \code{\link{iterateIntervals}} and
#'   \code{\link{identifyDuplicates}}
#' 
#' @export
analyzeIBDsegments <-
function(fileName, startRun = 1, endRun, shiftSize = 5000, 
         intervalSize = 10000) {
  countsA2 <- c()
  mergedIBDsegmentList <- list()
  dups <- c()
  snvs <- c()
  resHapRFN <- c()

  load(file = paste0(fileName, "_All", ".Rda"))
  load(file = "dups.Rda")
  
  if (startRun > 1 && length(countsA2[,3]) > 1) {
    tzz <- countsA2[which(countsA2[, 3] < startRun), ]
    offC <- max(tzz[, 1])
  } else {
    offC <-  0
  }
  
  avIBDsegmentPos <- c()
  avIBDsegmentLengthSNV <- c()
  avIBDsegmentLength <- c()
  avnoIndivid <- c()
  avnoTagSNVs <- c()
  
  avnoFreq <- c()
  avnoGroupFreq <- c()
  avnotagSNVChange <- c()
  avnotagSNVsPerIndividual <- c()
  avnoindividualPerTagSNV <- c()
  
  allCount <- 0
  allCount1 <- 0
  
  for (posAll in startRun:endRun) {
    start <- (posAll-1)*shiftSize
    end <- start + intervalSize
    
    if (end > snvs) {
      end <- snvs
    }
    
    pRange <- .createRangeString(start, end)
    
    load(file = paste0(fileName, pRange, "_resAnno.Rda"))
    
    mergedIBDsegmentList <- resHapRFN$mergedIBDsegmentList
    
    noIBDsegments <- lengthList(mergedIBDsegmentList)
    
    if (noIBDsegments > 0) {
      for (IBDsegmentC in 1:noIBDsegments) {
        
        allCount <- allCount + 1
        
        if (!dups[allCount + offC]) {
          
          allCount1 <- allCount1 + 1
          
          vt <- mergedIBDsegmentList[[IBDsegmentC]]
          
          avIBDsegmentPos <- c(avIBDsegmentPos, IBDsegmentPos(vt))
          avIBDsegmentLengthSNV <- c(avIBDsegmentLengthSNV, IBDsegmentLength(vt))
          avIBDsegmentLength <- c(avIBDsegmentLength,
                                  max(tagSNVPositions(vt)) - min(tagSNVPositions(vt)))
          avnoIndivid <- c(avnoIndivid, numberIndividuals(vt))
          avnoTagSNVs <- c(avnoTagSNVs, numbertagSNVs(vt))
          
          avnoFreq <- c(avnoFreq, tagSNVFreq(vt))
          avnoGroupFreq <- c(avnoGroupFreq, tagSNVGroupFreq(vt))
          avnotagSNVChange <- c(avnotagSNVChange, tagSNVChange(vt))
          avnotagSNVsPerIndividual <- c(avnotagSNVsPerIndividual, tagSNVsPerIndividual(vt))
          avnoindividualPerTagSNV <- c(avnoindividualPerTagSNV, individualPerTagSNV(vt))
        }
      }
    }
  }
  
  noIBDsegments <- allCount1
  
  avIBDsegmentPosS <- summary(avIBDsegmentPos)
  avIBDsegmentLengthSNVS <- summary(avIBDsegmentLengthSNV)
  avIBDsegmentLengthS <- summary(avIBDsegmentLength)
  avnoIndividS <- summary(avnoIndivid)
  avnoTagSNVsS <- summary(avnoTagSNVs)
  
  avnoFreqS <- summary(avnoFreq)
  avnoGroupFreqS <- summary(avnoGroupFreq)
  avnotagSNVChangeS <- summary(avnotagSNVChange)
  avnotagSNVsPerIndividualS <- summary(avnotagSNVsPerIndividual)
  avnoindividualPerTagSNVS <- summary(avnoindividualPerTagSNV)
  
  save(startRun, endRun, noIBDsegments, avIBDsegmentPos, avIBDsegmentLengthSNV, 
       avIBDsegmentLength, avnoIndivid, avnoTagSNVs, avnoFreq, avnoGroupFreq,
       avnotagSNVChange, avnotagSNVsPerIndividual, avnoindividualPerTagSNV, 
       avIBDsegmentPosS, avIBDsegmentLengthSNVS, avIBDsegmentLengthS, 
       avnoIndividS, avnoTagSNVsS, avnoFreqS, avnoGroupFreqS, avnotagSNVChangeS,
       avnotagSNVsPerIndividualS, avnoindividualPerTagSNVS, 
       file = paste0("analyzeResult", ".Rda"))
  
  return(list(startRun = startRun, endRun = endRun, noIBDsegments = noIBDsegments,
              avIBDsegmentPos = avIBDsegmentPos, avIBDsegmentLengthSNV = avIBDsegmentLengthSNV,
              avIBDsegmentLength = avIBDsegmentLength, avnoIndivid = avnoIndivid,
              avnoTagSNVs = avnoTagSNVs, avnoFreq = avnoFreq, avnoGroupFreq = avnoGroupFreq,
              avnotagSNVChange = avnotagSNVChange, 
              avnotagSNVsPerIndividual = avnotagSNVsPerIndividual,
              avnoindividualPerTagSNV = avnoindividualPerTagSNV,
              avIBDsegmentPosS = avIBDsegmentPosS, avIBDsegmentLengthSNVS = avIBDsegmentLengthSNVS,
              avIBDsegmentLengthS = avIBDsegmentLengthS, avnoIndividS = avnoIndividS,
              avnoTagSNVsS = avnoTagSNVsS, avnoFreqS = avnoFreqS, avnoGroupFreqS = avnoGroupFreqS,
              avnotagSNVChangeS = avnotagSNVChangeS,
              avnotagSNVsPerIndividualS = avnotagSNVsPerIndividualS,
              avnoindividualPerTagSNVS = avnoindividualPerTagSNVS))
}

#' @importFrom stats setNames
# Read info file and return a list with nsamples and nsnps
.readInfo <- 
function(prefixPath, fileName, infoPostfix) {
  setNames(as.list(as.numeric(readLines(paste0(prefixPath, fileName, infoPostfix), n = 2,
                                        warn = FALSE))), c("nsamples", "nsnps"))
}

# Read individuals file
.readIndividuals <- 
function(prefixPath, fileName, individualsPostfix) {
  read.table(paste0(prefixPath, fileName, individualsPostfix), 
      header = FALSE, sep = " ", quote = "", as.is = TRUE)
}

# Read individuals file
.readLabels <-
function(prefixPath, fileName, individualsPostfix, annotationFile, haplotypes,
         nsamples) {
  maxcol <- 4
  labels <- c()
  # If there is no annotation file
  if (is.null(annotationFile)) {
    # labelsAA id 1..n and sample names
    individuals <- .readIndividuals(prefixPath, fileName, individualsPostfix)
      
    # If there is only 1 individual name in the individuals file
    if (length(individuals[, 2]) >= 2) {
      col <- individuals[, 2]
    } else { # Else more than one individual in the file
      col <- 1:nsamples
    }
    if (haplotypes) {
      col <- rep(col, each = 2)
    }

    charcol <- as.character(col)
    matrix(rep(charcol, times = maxcol), ncol = maxcol)
  } else { # Annotation file exists
    annotations <- read.table(annotationFile, header = FALSE, sep = "\t", quote = "", as.is = TRUE)
    colnum <- ncol(annotations)

    # try to load without tabs
    if (colnum < maxcol) {
      annotations_space <- read.table(annotationFile, header = FALSE, sep = " ", quote = "",
                                      as.is = TRUE)
      if (ncol(annotations_space) > colnum) {
        annotations <- annotations_space
        colnum <- ncol(annotations_space)
      }
    }
    
    columns <- list()
    for (i in 1:min(colnum, maxcol)) {
      columns[[i]] <- annotations[, i]
    }
    if (min(colnum, maxcol) < maxcol) {
      individuals <- .readIndividuals(prefixPath, fileName, individualsPostfix)
      for (i in (min(colnum, maxcol) + 1):maxcol) {
        if (length(individuals[, 2]) >= 2) {
          columns[[i]] <- individuals[, 2]  
        } else {
          columns[[i]] <- 1:nsamples
        }
      }
    }
    if (haplotypes) {
      columns <- lapply(columns, function(x) rep(x, each=2))
    }
    for (i in 1:maxcol) {
      columns[[i]] <- as.character(columns[[i]])
      columns[[i]] <- gsub(",", ";", columns[[i]])
    }

    do.call(cbind, columns)
  }
}

#' @importFrom Matrix sparseMatrix
# Read sparse matrix from file
.readSparseMatrix <-
function(fileName) {
  con <- file(fileName, "r")
  lines <- readLines(con, warn = FALSE)

  rowPointer <- as.integer(unlist(strsplit(lines[3], " ", fixed = TRUE)))
  columnIndices <- as.integer(unlist(strsplit(lines[4], " ", fixed = TRUE)))
  values <- as.integer(unlist(strsplit(lines[5], " ", fixed = TRUE)))

  close(con)

  # Using i instead of j to transpose the matrix
  sparseMatrix(i = columnIndices, p = rowPointer, x = values, index1 = FALSE)
}

# Create range string (usually used as pRange)
.createRangeString <-
function(start, end) {
  paste0("_", format(start, scientific = FALSE), "_", format(end, scientific = FALSE))
}
