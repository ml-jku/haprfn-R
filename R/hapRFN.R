samplesPerFeature <- function(X, samples = 0, lowerB = 0, upperB = 1000) {
  if (missing(X)) {
    stop("Data file name missing. Stopped.")
  }
  
  samples <- as.integer(sort.int(as.integer(unique(samples))))
  lowerB <- as.double(lowerB)
  upperB <- as.double(upperB)
  
  res <- .Call("samplesPerFeature", X, samples, lowerB, upperB, PACKAGE = "hapRFN")
  
  # TODO: Return nil value on error, or a standard empty?
  
  # if (is.null(res)) { return list(sL = list(), nsL = 0) }
  
  return(res)
}

readSamplesSpRfn <- function(X, samples = 0, lowerB = 0, upperB = 1000) {
  if (missing(X)) {
    stop("Data file name missing. Stopped.")
  }
  
  samples <- as.integer(sort.int(as.integer(unique(samples))))
  lowerB <- as.double(lowerB)
  upperB <- as.double(upperB)
  
  res <- .Call("readSamplesSpRfn", X, samples, lowerB, upperB, PACKAGE = "hapRFN")
  
  # TODO: Return nil value on error, or a standard empty?
  
  # if (is.null(res)) { return list(X = as.matrix(1)) }
  
  return(res)
}

vcf2sparse <- function(fileName, prefixPath = NULL, intervalSize = 10000, shiftSize = 5000, 
                       annotation = TRUE, outputFile = NULL) {
  message("Running 'vcf2sparse' on ", fileName)
  message("   Path to file ----------------------- : ", prefixPath)
  if (!is.null(outputFile)) {
    message("   Output file prefix------------------ : ", outputFile)
  } else {
    message("   Output file prefix given by input ----")
  }

  .Call("vcf2sparse", fileName, prefixPath, as.integer(intervalSize), as.integer(shiftSize), annotation, outputFile, 
        PACKAGE = "hapRFN")

  message("")
  message("Convert End.")
}

iterateIntervals <- function(startRun = 1, endRun, shift = 5000, intervalSize = 10000, 
  annotationFile = NULL, fileName, prefixPath = "", sparseMatrixPostfix = "_mat", 
  annotPostfix = "_annot.txt", individualsPostfix = "_individuals.txt", individuals = 0, 
  lowerBP = 0, upperBP = 0.05, p = 10, iter = 40, quant = 0.01, eps = 1e-05, alpha = 0.03, 
  cyc = 50, non_negative = 1, write_file = 0, norm = 0, lap = 100, IBDsegmentLength = 50, 
  Lt = 0.1, Zt = 0.2, thresCount = 1e-05, mintagSNVsFactor = 3/4, pMAF = 0.03, 
  haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1, thresPrune = 0.001, simv = "minD", 
  minTagSNVs = 6, minIndivid = 2, avSNVsDist = 100, SNVclusterLength = 100) {

  labelsA <- c()
  annot <- c()
  
  ina <- as.numeric(readLines(paste(prefixPath, fileName, sparseMatrixPostfix, 
                    ".txt", sep = ""), n = 2))
  # These information are missing, have to change this.
  individualsN <- ina[1]
  snvs <- ina[2]
  
  save(individualsN, snvs, file = paste(fileName, "_All", ".Rda", sep = ""))
  
  for (posAll in startRun:endRun) {
    start <- (posAll - 1) * shifs
    end <- start + intervalSize

    if (end > snvs) {
      end <- snvs
    }
    
    pRange <- paste("_", format(start, scientific = FALSE), "_", format(end, 
      scientific = FALSE), sep = "")
    
	# If there is no annotation file
    if (is.null(annotationFile)) {
	
      # Read individuals file
      labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
        sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
        
      # If there is only 1 individual in the individuals file
      if (length(labelsAA[, 2]) < 2) {
        if (haplotypes) { # If phased genotypes
          lA <- as.vector(unlist(rbind(1:individualsN, 1:individualsN)))
        } else { # Else unphased genotypes
          lA <- as.vector(1:individualsN)
        }
      } else { # Else more than one individual in the file
        if (haplotypes) {
          lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
        } else {
          lA <- as.vector(labelsAA[, 2])
        }
      }
      indiA <- cbind(as.character(lA), as.character(lA), as.character(lA), 
        as.character(lA))
    } else { # Annotation file exists
      indit <- read.table(annotationFile, header = FALSE, sep = "\t", quote = "", 
        as.is = TRUE)
      lind <- length(indit)
      if (lind < 4) {
        inditA <- read.table(annotationFile, header = FALSE, sep = " ", quote = "", 
          as.is = TRUE)
        if (length(inditA) > lind) {
          indit <- inditA
          lind <- length(inditA)
        }
      }
      
      if (haplotypes) {
        if (lind > 0) {
          indi1 <- as.vector(unlist(rbind(indit[, 1], indit[, 1])))  # because haplotypes individuals are doubled
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(unlist(rbind(1:individualsN, 1:individualsN)))
          } else {
            lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
          }
          indi1 <- as.character(lA)
        }
        if (lind > 1) {
          indi2 <- as.vector(unlist(rbind(indit[, 2], indit[, 2])))
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(unlist(rbind(1:individualsN, 1:individualsN)))
          } else {
            lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
          }
          indi2 <- as.character(lA)
        }
        if (lind > 2) {
          indi3 <- as.vector(unlist(rbind(indit[, 3], indit[, 3])))
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(unlist(rbind(1:individualsN, 1:individualsN)))
          } else {
            lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
          }
          indi3 <- as.character(lA)
        }
        if (lind > 3) {
          indi4 <- as.vector(unlist(rbind(indit[, 4], indit[, 4])))
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(unlist(rbind(1:individualsN, 1:individualsN)))
          } else {
            lA <- as.vector(unlist(rbind(labelsAA[, 2], labelsAA[, 2])))
          }
          indi4 <- as.character(lA)
        }
      } else {
        if (lind > 0) {
          indi1 <- as.vector(indit[, 1])
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(1:individualsN)
          } else {
            lA <- as.vector(labelsAA[, 2])
          }
          indi1 <- as.character(lA)
        }
        if (lind > 1) {
          indi2 <- as.vector(indit[, 2])
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(1:individualsN)
          } else {
            lA <- as.vector(labelsAA[, 2])
          }
          indi2 <- as.character(lA)
        }
        if (lind > 2) {
          indi3 <- as.vector(indit[, 3])
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(1:individualsN)
          } else {
            lA <- as.vector(labelsAA[, 2])
          }
          indi3 <- as.character(lA)
        }
        if (lind > 3) {
          indi4 <- as.vector(indit[, 4])
        } else {
          labelsAA <- read.table(paste(prefixPath, fileName, individualsPostfix, 
            sep = ""), header = FALSE, sep = " ", quote = "", as.is = TRUE)
          if (length(labelsAA[, 2]) < 2) {
            lA <- as.vector(1:individualsN)
          } else {
            lA <- as.vector(labelsAA[, 2])
          }
          indi4 <- as.character(lA)
        }
      }
      indi1 <- gsub(",", ";", indi1)
      indi2 <- gsub(",", ";", indi2)
      indi3 <- gsub(",", ";", indi3)
      indi4 <- gsub(",", ";", indi4)
      indiA <- cbind(indi1, indi2, indi3, indi4)
      labelsA <- indiA
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
    
    resHapFabia <- hapFabiaNew(fileName = fileName, prefixPath = prefixPath, 
      sparseMatrixPostfix = sparseMatrixPostfix, annotPostfix = annotPostfix, 
      individualsPostfix = individualsPostfix, labelsA = labelsA, pRange = pRange, 
      individuals = individuals, lowerBP = lowerBP, upperBP = upperBP, p = p, 
      iter = iter, quant = quant, eps = eps, alpha = alpha, cyc = cyc, non_negative = non_negative, 
      write_file = write_file, norm = norm, lap = lap, IBDsegmentLength = IBDsegmentLength, 
      Lt = Lt, Zt = Zt, thresCount = thresCount, mintagSNVsFactor = mintagSNVsFactor, 
      pMAF = pMAF, haplotypes = haplotypes, cut = cut, procMinIndivids = procMinIndivids, 
      thresPrune = thresPrune, simv = simv, minTagSNVs = minTagSNVs, minIndivid = minIndivid, 
      avSNVsDist = avSNVsDist, SNVclusterLength = SNVclusterLength)
    
    IBDsegmentList2excel(resHapFabia$mergedIBDsegmentList, paste(fileName, pRange, 
      ".csv", sep = ""))
    
    save(resHapFabia, annot, file = paste(fileName, pRange, "_resAnno", ".Rda", 
      sep = ""))
  }
}

hapFabiaNew <- function(fileName, prefixPath = "", sparseMatrixPostfix = "_mat", 
  annotPostfix = "_annot.txt", individualsPostfix = "_individuals.txt", labelsA = NULL, 
  pRange = "", individuals = 0, lowerBP = 0, upperBP = 0.05, p = 10, iter = 40, 
  quant = 0.01, eps = 1e-05, alpha = 0.03, cyc = 50, non_negative = 1, write_file = 0, 
  norm = 0, lap = 100, IBDsegmentLength = 50, Lt = 0.1, Zt = 0.2, thresCount = 1e-05, 
  mintagSNVsFactor = 3/4, pMAF = 0.03, haplotypes = FALSE, cut = 0.8, procMinIndivids = 0.1, 
  thresPrune = 0.001, simv = "minD", minTagSNVs = 6, minIndivid = 2, avSNVsDist = 100, 
  SNVclusterLength = 100) {
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
  
  
  message("                      ")
  message("                      ")
  message("Running hapFabia with:")
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
  
  ina <- as.numeric(readLines(paste(prefixPath, fileName, pRange, sparseMatrixPostfix, 
    ".txt", sep = ""), n = 2))
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
  
  # build DgC Matrix
  
  filename1 <- paste0(prefixPath, fileName, pRange, sparseMatrixPostfix, ".txt")
  
  filenamelist <- list(filename1)
  
  # What is this?
  res1 <- readsparsematrixlist(filenamelist)
  
  # What is this?
  X <- fillmat2(res1)
  
  
  message("start RFN")
  
  # Fabia call res <-
  # spfabia(X=paste(prefixPath,fileName,pRange,sparseMatrixPostfix,sep=''),p=p,alpha=alpha,cyc=cyc,non_negative=non_negative,write_file=write_file,norm=norm,lap=lap,samples=individuals,iter=iter,quant=quant,lowerB=lowerBindivid,upperB=upperBindivid,eps=eps)
  
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
    minP = 0.01, seed = 0)
  
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
