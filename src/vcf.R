source('src/dsg.R')

## get sequence
get.seq <- function(vcf.dir, chr, bp1, bp2, ssn = NA, gtp = 'GT')
{
    suppressPackageStartupMessages(library(
        VariantAnnotation, warn.conflicts = F, quietly = T))

    rng <- GRanges(
        seqnames = sprintf("%d", chr),
        ranges = IRanges(bp1, bp2, names = ssn))
    
    svp <- ScanVcfParam(
        fixed = 'ALT',
        info = NA,
        geno = gtp,
        which = rng)

    tbx <- TabixFileList(dir(vcf.dir, 'vcf.gz$', full.names = T))
    vgz <- tbx[[sprintf("c%02d.vcf.gz", chr)]]

    sink('/dev/null')
    vcf <- try(readVcf(vgz, 'hg38', param = svp), silent = T)
    sink()
    vcf
}

## convert vcf GT to dosage
vcf2dsg <- function(vcf)
{
    suppressPackageStartupMessages(library(
        VariantAnnotation, warn.conflicts = F, quietly = T))

    ## variant map & genotype matrix
    gmp <- as.data.frame(rowRanges(vcf))
    gmx <- geno(vcf)$GT
    
    ## variant id and subject id
    vid <- rownames(gmx)
    sid <- colnames(gmx)

    ## drop bioconductor structure for the variant map
    gmp <- with(gmp, data.frame(
        vid = vid,
        chr=seqnames, bp1=start, bp2=end,
        ref=REF, alt=sapply(gmp$ALT, toString),
        stringsAsFactors = F))

    ## convert vcf.GT  format to dosage format 
    ## save the # of variants and subjects
    N <- dim(gmx)[1];
    M <- dim(gmx)[2]

    ## "a1|a2" ==> {0, 1, 2}
    gmx <- sapply(strsplit(gmx, '|'), `[`, c(1L,3L))
    gmx[gmx=='.'] <- NA
    gmx <- as.integer(gmx)
    dim(gmx) <- c(2L, N * M)
    gmx[gmx > 0L] <- 1L
    gmx <- colSums(gmx)
    dim(gmx) <- c(N, M)
    dimnames(gmx) <- list(vid=vid, sid=sid)

    ## construct dosage data structure
    dosage(gmp=gmp, gmx=gmx)
}
