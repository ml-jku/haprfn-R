aa <- function() {

library(hapFabia)
library(Matrix)
library(RFN)

#for (chr in 15:16) {

    chr <- 22

    wd <- paste0("/system/user/tigerdata/data/1000genomes_vcf/phase3v5/Unrelated/Filtered/HapRFN/chr", chr, "_gpu")
    dir.create(wd)
    setwd(wd)

    source("~/scripts/1000g_relations_functions.R")
    source("~/HapRFN.R")

    
    fileName <- paste0("1000G.chr", chr, ".integrated_phase3_v5")
	intervalSize <- 10000
	shiftSize <- 5000
    prefixPath <- paste0("/system/user/tigerdata/data/1000genomes_vcf/phase3v5/Unrelated/Filtered/HapFABIA/chr", chr, "/")
    annotFile <- "/system/user/tigerdata/data/1000genomes_vcf/phase3v5/Unrelated/HapFABIA/individuals.txt"


#####compute how many segments we have#######
    ina <- as.numeric(readLines(paste0(prefixPath, fileName, "_mat.txt"),n=2))
    noSNVs <- ina[2]
    over <- intervalSize%/%shiftSize
    N1 <- noSNVs%/%shiftSize
    endRunA <- (N1-over+2)

#####parameters###########
    startRun=1
    endRun=endRunA
#    endRun=1
    annotationFile=annotFile
    individuals=0
    lowerBP=0
    upperBP=0.05
    iter=1
    quant=0.01
    eps=1e-5
    alpha=0.03
    IBDsegmentLength=50 #in kbp
    Lt = 0.1
    Zt = 0.2
    thresCount=1e-5
    mintagSNVsFactor=3/4
    pMAF=0.03
    haplotypes=FALSE
    cut=0.8
    procMinIndivids=0.1
    thresPrune=1e-3
    minTagSNVs=6
    minIndivid=2
    avSNVsDist=100
    SNVclusterLength=100
    set.seed(188)
    
    p=50
    cyc=100
#    l1=0.0001
    l1=0

#endRun<-1

    #####analyze each segment#######
    system.time(
    iterateIntervalsNew(startRun=startRun,endRun=endRun,shift=shiftSize,intervalSize=intervalSize,
            annotationFile=annotationFile,prefixPath=prefixPath,fileName=fileName,
            individuals=individuals,lowerBP=lowerBP,upperBP=upperBP,p=p,iter=iter,quant=quant,
            eps=eps,alpha=alpha,cyc=cyc,Lt = Lt,Zt = Zt,thresCount=thresCount,
            mintagSNVsFactor=mintagSNVsFactor,IBDsegmentLength=IBDsegmentLength,pMAF=pMAF,cut=cut,
            thresPrune=thresPrune,minTagSNVs=minTagSNVs,procMinIndivids=procMinIndivids,
            minIndivid=minIndivid,avSNVsDist=avSNVsDist,SNVclusterLength=SNVclusterLength,
            haplotypes=haplotypes, l1=l1, use_gpu=TRUE, gpu_id=3)
    )
    
    identifyDuplicates(fileName=fileName,startRun=1,endRun=endRun,
        shift=shiftSize,intervalSize=intervalSize)

    anaRes <- analyzeIBDsegments(fileName=fileName,startRun=1,endRun=endRun,
        shift=shiftSize,intervalSize=intervalSize)

    save(anaRes, file=paste0("anaRes.chr", chr, ".Rda"))
}