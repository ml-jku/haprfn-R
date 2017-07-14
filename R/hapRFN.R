samplesPerFeature <- function(X, samples = 0, lowerB = 0, upperB = 1000) {
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

readSparseSamples <- function(X, samples = 0, lowerB = 0, upperB = 1000) {
  if (missing(X)) {
    stop("Data file name missing. Stopped.")
  }
  
  X <- as.character(X)
  samples <- as.integer(sort.int(as.integer(unique(samples))))
  lowerB <- as.double(lowerB)
  upperB <- as.double(upperB)
  
  .Call("readSparseSamples", X, samples, lowerB, upperB, PACKAGE = "hapRFN")
}


## Missing values ->
#       0: default: look for frequency defined in INFO as AF or calculate
#       1: major allele
#       2: minor allele
#       3: look for AF field or abort
#       4: calculate
#       5: abort on missing value
vcf2sparse <- function(fileName, prefixPath = NULL, intervalSize = 10000, shiftSize = 5000, 
                       annotation = TRUE, genotypes = TRUE, haplotypes = FALSE, missingValues = 0, 
                       outputFile = NULL) {
  .Call("vcf2sparse", fileName, prefixPath, as.integer(intervalSize), as.integer(shiftSize),
        annotation, genotypes, haplotypes, as.integer(missingValues), outputFile, 
        PACKAGE = "hapRFN")
}

readInfo <- function(prefixPath, fileName, infoPostfix) {
  setNames(as.list(as.numeric(readLines(paste0(prefixPath, fileName, infoPostfix), n = 2))), c("nsamples, nsnps"))
}

readIndividuals <- function(prefixPath, fileName, individualsPostfix) {
  read.table(paste0(prefixPath, fileName, individualsPostfix), 
      header = FALSE, sep = " ", quote = "", as.is = TRUE)
}

readLabels <- function(prefixPath, fileName, individualsPostfix, annotationFile, haplotypes, nsamples) {
  maxcol <- 4
  labels <- c()
  # If there is no annotation file
  if (is.null(annotationFile)) {
    # labelsAA id 1..n and sample names
    individuals <- readIndividuals(prefixPath, fileName, individualsPostfix)
      
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
      annotations_space <- read.table(annotationFile, header = FALSE, sep = " ", quote = "", as.is = TRUE)
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
      individuals <- readIndividuals(prefixPath, fileName, individualsPostfix)
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

iterateIntervals <- function(startRun = 1, endRun, shift = 5000, intervalSize = 10000, 
  annotationFile = NULL, fileName, prefixPath = "", sparseMatrixPostfix = "_mat.txt", 
  annotPostfix = "_annot.txt", individualsPostfix = "_individuals.txt", 
  infoPostfix = "_info.txt", individuals = 0, l1 = 0.0,
  lowerBP = 0, upperBP = 0.05, p = 10, iter = 40, quant = 0.01, eps = 1e-05, alpha = 0.03, 
  cyc = 50, non_negative = 1, write_file = 0, norm = 0, lap = 100, IBDsegmentLength = 50, 
  Lt = 0.1, Zt = 0.2, thresCount = 1e-05, mintagSNVsFactor = 3/4, pMAF = 0.03, 
  haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1, thresPrune = 0.001, simv = "minD", 
  minTagSNVs = 6, minIndivid = 2, avSNVsDist = 100, SNVclusterLength = 100, saveAsCsv = FALSE,
  use_gpu = TRUE, gpu_id = 0) {

  info <- readInfo(prefixPath, fileName, infoPostfix)

  nsamples <- info$nsamples
  snvs <- info$nsnps
  
  save(nsamples, snvs, file = paste0(fileName, "_All", ".Rda"))

  for (posAll in startRun:endRun) {
    start <- (posAll - 1) * shift
    end <- start + intervalSize

    if (end > snvs) {
      end <- snvs
    }
        
    labels <- readLabels(prefixPath, fileName, individualsPostfix, annotationFile, haplotypes, nsamples)
    pRange <- paste0("_", format(start, scientific = FALSE), "_", format(end, scientific = FALSE))
    
    resHapRFN <- hapRFN(fileName = fileName, prefixPath = prefixPath, 
      sparseMatrixPostfix = sparseMatrixPostfix, annotPostfix = annotPostfix, 
      individualsPostfix = individualsPostfix, labelsA = labels, pRange = pRange, 
      individuals = individuals, lowerBP = lowerBP, upperBP = upperBP, p = p, 
      iter = iter, quant = quant, eps = eps, alpha = alpha, cyc = cyc, non_negative = non_negative, 
      write_file = write_file, norm = norm, lap = lap, IBDsegmentLength = IBDsegmentLength, 
      Lt = Lt, Zt = Zt, thresCount = thresCount, mintagSNVsFactor = mintagSNVsFactor, 
      pMAF = pMAF, haplotypes = haplotypes, cut = cut, procMinIndivids = procMinIndivids, 
      thresPrune = thresPrune, simv = simv, minTagSNVs = minTagSNVs, minIndivid = minIndivid, 
      avSNVsDist = avSNVsDist, SNVclusterLength = SNVclusterLength, use)
    
    if (saveAsCsv) {
      IBDsegmentList2excel(resHapRFN$mergedIBDsegmentList, paste0(fileName, pRange, ".csv"))  
    }
    
    annot <- c()
    save(resHapRFN, annot, file = paste0(fileName, pRange, "_resAnno", ".Rda"))
  }
}

readSparseMatrix <- function(fileName) {
  con <- file(fileName, "r")
  lines <- readLines(con, warn = FALSE)

  nnz <- as.integer(lines[1])
  nrow <- as.integer(lines[2])
  rowPointer <- as.integer(unlist(strsplit(lines[3], " ", fixed = TRUE)))
  columnIndices <- as.integer(unlist(strsplit(lines[4], " ", fixed = TRUE)))
  values <- as.integer(unlist(strsplit(lines[5], " ", fixed = TRUE)))

  sparseMatrix(j = columnIndices, p = rowPointer, x = values, index1 = FALSE)
}

hapRFN <- function(fileName, prefixPath = "", sparseMatrixPostfix = "_mat", 
  annotPostfix = "_annot.txt", individualsPostfix = "_individuals.txt", labelsA = NULL, 
  pRange = "", individuals = 0, lowerBP = 0, upperBP = 0.05, p = 10, iter = 40, 
  quant = 0.01, eps = 1e-05, l1 = 0.0, alpha = 0.03, cyc = 50, non_negative = 1, write_file = 0, 
  norm = 0, lap = 100, IBDsegmentLength = 50, Lt = 0.1, Zt = 0.2, thresCount = 1e-05, 
  mintagSNVsFactor = 3/4, pMAF = 0.03, haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1, 
  thresPrune = 0.001, simv = "minD", minTagSNVs = 6, minIndivid = 2, avSNVsDist = 100, 
  SNVclusterLength = 100, gpu = FALSE, gpuId = -1) {
  # fileName:            the file name of the sparse matrix in sparse format.
  # prefixPath:          path of the data file
  # sparseMatrixPostfix: postfix string for the sparse matrix
  # annotPostfix:        postfix string for the annotation file
  # labelsA:             individual names as matrix individuals x 4
  # prange:              for intervals indicates its range
  # individuals:             vector of individuals which should be included into the analysis; default = 0 (all individuals)
  # lowerBP:             lower bound for filtering the inputs columns, minimal MAF (however more than one occurence to remove private SNVs)
  # upperBP:             Upper bound for filtering the inputs columns, minimal MAF
  # p:                   no biclusters per iteration
  # alpha:               sparseness loadings; default = 0.03
  # iter:                number iterations
  # quant:               percentage of Ls to remove in each iteration
  # eps:                 lower bound for variational parameter lapla; default: 1e-5
  # l1:                  l1 weight decay
  # cyc:                 number of iterations; default = 50
  # non_negative:        Non-negative factors and loadings if non_negative; default = 1 (yes).
  # write_file:          results are written to files (L in sparse format), default = 0 (not written).
  # norm:                data normalization; default = 1 (no normalization).
  # lap:                 minimal value of the variational parameter; default = 100.0.
  # IBDsegmentLength:           IBD segment length in kbp
  # Lt:                  percentage of largest Ls to consider for IBD segment extraction
  # Zt:                  percentage of largest Zs to consider for IBD segment extraction
  # thresCount:          p-value of random histogram hit, default 1e-5
  # mintagSNVsFactor:       percentage of segments overlap in IBD segments; 1/2 for large to 3/4 for small intervals
  # pMAF:                averaged and corrected minor allele frequency
  # haplotypes:          haplotypes = phased genotypes -> two chromosomes per individual
  # gpu:                 TRUE to use GPU. Use this parameter together with gpuId
  # gpuId:               the ID of the GPU. This is used, when gpu is TRUE
  
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
  if (length(individuals) > 1) {
    message("   Number of individuals included into the analysis ---: ", length(individuals))
  } else {
    message("   All individuals are included into the analysis -----: 0 = all individuals")
  }
  message("   Lower bound on MAF but more than one occurence -----: ", lowerBP)
  message("   Upper bound on MAF ---------------------------------: ", upperBP)
  message("   Number of biclusters per iteration -----------------: ", p)
  message("   Sparseness coefficient of the loadings -------------: ", alpha)
  message("   Number of iterations -------------------------------: ", iter)
  message("   Percentage of Ls to remove after each iteration ----: ", quant)
  message("   Lower bound for variational parameter lapla --------: ", eps)
  message("   Number of cycles per iteration ---------------------: ", cyc)
  message("   Non-negative factors and loadings ------------------: ", non_negative)
  message("   Results written to files ---------------------------: ", write_file)
  message("   Data normalized ------------------------------------: ", norm)
  message("   Minimal value of the variational parameter ---------: ", lap)
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
  message("   % largest Ls for IBD segment extraction ------------: ", Lt)
  message("   % largest Zs for IBD segment extraction ------------: ", Zt)
  message("   p-value threshold for random histogram counts ------: ", thresCount)
  message("   Min. % of segments overlap in IBD segments ---------: ", mintagSNVsFactor)
  message("   Averaged and corrected minor allele frequency ------: ", pMAF)
  message("   Cutoff for merging IBD segments --------------------: ", cut)
  message("   % of cluster individuals a tagSNV must tag ---------: ", procMinIndivids)
  message("   Threshold for pruning border tagSNVs ---------------: ", thresPrune)
  message("   Similarity measure for merging clusters ------------: ", simv)
  message("   Minimum matching tagSNVs for cluster similarity ----: ", minTagSNVs)
  message("   Minimum matching individuals for cluster similarity : ", minIndivid)
  message("                      ")
  message("                      ")
  
  # Maybe remove this when dependency added?
  require("fabia")
  
  # Compute internal parameters
  
  # ps: quantile above which to consider Ls
  ps <- 1 - Lt
  # psZ: quantile above which to consider Zs
  psZ <- 1 - Zt
  
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
  
  ina <- as.numeric(readLines(paste0(prefixPath, fileName, pRange, sparseMatrixPostfix, ".txt"), n = 2))
  if (length(individuals) > 1) {
    individualsN <- length(individuals)
  } else {
    individualsN <- ina[1]
    individuals <- 0
  }
  snvs <- ina[2]
  
  upperBindivid = upperBP * individualsN  # remove common SNVs
  lowerBindivid = max(1.5, lowerBP * individualsN)  # remove private SNVs
  
  kk <- 1
  while ((snvs/inteA) * choose(individualsN, 2) * (1 - pbinom(kk, inteA, pMAF * 
      pMAF)) > thresCount) {
    kk <- kk + 1
  }
  
  thresA <- kk
  mintagSNVs <- round(mintagSNVsFactor * thresA)
  
  # End Compute internal parameters
    
  matrixFileName <- paste0(prefixPath, fileName, pRange, sparseMatrixPostfix, ".txt")
  X <- readSparseMatrix(matrixFileName)

  message("start RFN")
    
  # RFN call
  l = ncol(X)
  n = nrow(X)
  
  rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "SNV")
  colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
  rowna <- rownames(X)
  colna <- colnames(X)
  
  eps <- as.double(0.001)
  eps1 <- as.double(1e-10)
  iin <- 1/l
  
  cyc <- as.integer(cyc)
  
  p <- as.integer(p)
  
  com <- which(table(X@i)/X@Dim[2] > upperBP)
  X[com + 1, ] <- 0
  
  rfn_res <- train_rfn(X = X, n_hidden = p, n_iter = cyc, etaW = 0.1, etaP = 0.1, 
    minP = 0.01, l1_weightdecay = l1, seed = 0, use_gpu = gpu, gpu_id = gpuId)
  
  myL = rfn_res$W
  myPsi = as.vector(rfn_res$P)
  myZ = rfn_res$H
  mylapla = as.matrix(1)
  
  vz <- iin * apply(myZ, 1, function(x) sum(x^2))
  vz <- sqrt(vz + 1e-10)
  ivz <- 1/vz
  if (length(ivz) == 1) {
    nZ <- ivz * myZ
    noL <- vz * myL
  } else {
    nZ <- ivz * myZ
    noL <- t(vz * t(myL))
  }
  
  Lz = noL %*% nZ
  rownames(noL) <- rowna
  colnames(noL) <- colnames(noL, do.NULL = FALSE, prefix = "bicluster")
  clnames <- colnames(noL)
  rownames(nZ) <- clnames
  colnames(nZ) <- colna
  
  ini <- as.matrix(1)
  avini <- as.vector(1)
  xavini <- as.vector(1)
  
  res <- new("Factorization", parameters = list("rfn", cyc, p), n = n, p1 = p, 
    p2 = p, l = l, center = as.vector(1), scaleData = as.vector(1), X = as.matrix(1), 
    L = noL, Z = nZ, M = as.matrix(1), LZ = as.matrix(1), U = as.matrix(1), avini = avini, 
    xavini = xavini, ini = ini, Psi = myPsi, lapla = mylapla)
  
  # Load individuals to Ls of interest: load minor alleles of the Ls
  
  sPF <- samplesPerFeature(X = paste(prefixPath, fileName, pRange, sparseMatrixPostfix, sep = ""), 
                           samples = individuals, lowerB = lowerBindivid, upperB = upperBindivid)

  if (nchar(annotPostfix) > 0) {
    # annot[[1]] <- chromosome annot[[2]] <- phys. position annot[[3]] <- snvNames
    # annot[[4]] <- snvMajor annot[[5]] <- snvMinor annot[[6]] <- quality annot[[7]]
    # <- pass annot[[8]] <- info of vcf file annot[[9]] <- fields in vcf file
    # annot[[10]] <- frequency annot[[11]] <- 1 = changed if major allele is actually
    # minor allele otherwise 0
    
    annot <- read.table(paste(prefixPath, fileName, pRange, annotPostfix, sep = ""), 
      header = FALSE, sep = "\t", quote = "", as.is = TRUE, skip = 2)
    
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
    labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, sep = ""), 
      header = FALSE, sep = " ", quote = "", as.is = TRUE)
    if (haplotypes) {
      lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
    } else {
      lA <- as.vector(labelsAA[, 2])
    }
    indiA <- cbind(as.character(lA), as.character(lA), as.character(lA), as.character(lA))
  } else {
    indiA <- labelsA
  }
  
  # Here other labels might be possible: labelsA[,i]: 1=id 2=subPopulation
  # 3=population 4 = platform
  
  # indit <- read.table('phase1_integrated_calls.20101123.ALL.panel', header =
  # FALSE, sep = '\t', quote = '',as.is=TRUE) indi1 <-
  # as.vector(unlist(rbind(indit[,1],indit[,1]))) # because haplotypes individuals
  # are doubled indi2 <- as.vector(unlist(rbind(indit[,2],indit[,2]))) indi3 <-
  # as.vector(unlist(rbind(indit[,3],indit[,3]))) indi4 <-
  # as.vector(unlist(rbind(indit[,4],indit[,4]))) indi4 <- gsub(',',';',indi4)
  # indiA <- cbind(indi1,indi2,indi3,indi4)
  
  
  
  # save fabia result
  # save(res,sPF,annot,file=paste(fileName,pRange,'_res.Rda',sep=''))

  # first haplotype extraction with offset 0
  
  off1 <- 0

  IBDsegmentList1 <- extractIBDsegments(res = res, sPF = sPF, annot = annot, chrom = "", 
    labelsA = indiA, ps = ps, psZ = psZ, inteA = inteA, thresA = thresA, mintagSNVs = mintagSNVs, 
    off = off1, procMinIndivids = procMinIndivids, thresPrune = thresPrune)
  
  # merge IBD segment lists
  
  if (lengthList(IBDsegmentList1) > 1) {
    comp <- compareIBDsegmentLists(IBDsegmentList1 = IBDsegmentList1, IBDsegmentList2 = NULL, 
      simv = simv, pTagSNVs = NULL, pIndivid = NULL, minTagSNVs = minTagSNVs, 
      minIndivid = minIndivid)
    
    if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp, h = cut)
      mergedIBDsegmentList1 <- mergeIBDsegmentLists(IBDsegmentList1 = IBDsegmentList1, 
        IBDsegmentList2 = NULL, clustIBDsegmentList = clustIBDsegmentList)
    } else {
      mergedIBDsegmentList1 <- IBDsegmentList1
    }
  } else {
    mergedIBDsegmentList1 <- IBDsegmentList1
  }
  
  
  # second IBD extraction with offset half of the interval length
  
  off2 = inteA %/% 2
  
  IBDsegmentList2 <- extractIBDsegments(res = res, sPF = sPF, annot = annot, chrom = "", 
    labelsA = indiA, ps = ps, psZ = psZ, inteA = inteA, thresA = thresA, mintagSNVs = mintagSNVs, 
    off = off2, procMinIndivids = procMinIndivids, thresPrune = thresPrune)
  
  # merge IBD segments
  
  if (lengthList(IBDsegmentList2) > 1) {
    comp <- compareIBDsegmentLists(IBDsegmentList1 = IBDsegmentList2, IBDsegmentList2 = NULL, 
      simv = simv, pTagSNVs = NULL, pIndivid = NULL, minTagSNVs = minTagSNVs, 
      minIndivid = minIndivid)
    
    if (!is.null(comp)) {
      clustIBDsegmentList <- cutree(comp, h = cut)
      mergedIBDsegmentList2 <- mergeIBDsegmentLists(IBDsegmentList1 = IBDsegmentList2, 
        IBDsegmentList2 = NULL, clustIBDsegmentList = clustIBDsegmentList)
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
          IBDsegmentList2 = mergedIBDsegmentList2, clustIBDsegmentList = clustIBDsegmentList)
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
  
  # save(mergedIBDsegmentList,res,sPF,annot,IBDsegmentList1,IBDsegmentList2,mergedIBDsegmentList1,mergedIBDsegmentList2,file=paste(fileName,pRange,'_IBDsegmentList.Rda',sep=''))
  
  return(list(mergedIBDsegmentList = mergedIBDsegmentList, res = res, sPF = sPF, 
    annot = annot, IBDsegmentList1 = IBDsegmentList1, IBDsegmentList2 = IBDsegmentList2, 
    mergedIBDsegmentList1 = mergedIBDsegmentList1, mergedIBDsegmentList2 = mergedIBDsegmentList2))
}
