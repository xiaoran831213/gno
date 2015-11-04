library(VariantAnnotation)
source('src/ped.R')

## get gene list
get.gen <- function()
{
    gen <- read.table("dat/gen.txt", header=T, as.is = T)
    names(gen) <- tolower(names(gen))
    with(gen, {bp1 <- bp1 - 5000; bp2 <- bp2 + 5000})
    gen
}

## get sequence
get.seq <- function(chr, bp1, bp2, gen)
{
    rng <- GRanges(
        seqnames = sprintf("%d", chr),
        ranges = IRanges(bp1, bp2, names = gen))
        
    svp <- ScanVcfParam(
        fixed = 'ALT',
        info = NA,
        geno = c("GT", "BD"),
        which = rng)

    tbx <- TabixFileList(dir('wgs', 'vcf.gz$', full.names = T))
    vcf <- readVcf(tbx[[sprintf("c%02d.vcf.gz", chr)]], 'hg38', param = svp)
    vcf
}

## convert vcf GT to dosage
vgt2dsg <- function(vgt)
{
    g <- sapply(strsplit(vgt, '|'), `[`, c(1,3))
    g[g=='.'] <- NA
    g <- as.integer(g)
    dim(g) <- c(2L, length(vgt))
    g[g > 0] <- 1
    g <- colSums(g)
    dim(g) <- dim(vgt)
    g
}

## process one sequence
proc.seq <- function(vcf, ped)
{
    ## get sequence map
    gmp <- as.data.frame(rowRanges(vcf))
    
    ## variant id
    vid <- rownames(gmp)
    gmp <- with(gmp, data.frame(
        vid = vid,
        chr=seqnames, bp1=start, bp2=end,
        ref=REF, alt=sapply(gmp$ALT, toString),
        stringsAsFactors = F))

    ## subject id must be in pedigree
    sid <- intersect(samples(header(vcf)), as.character(ped$iid))
    
    ## dosage from GT, and Beagal
    dsg <- vgt2dsg(geno(vcf)$GT[, sid])
    dsb <- geno(vcf)$BD[, sid]
    dimnames(dsg) <- list(vid=vid, sid=sid)
    dimnames(dsb) <- dimnames(dsg)
    
    ## check mendilian errors
    mde <- check.ME(vcf, ped)

    ## set error genotype to NA
    err_idx <- with(
        mde,
    {
        err <- cbind(vid=vid, sid=c(iid, pid, mid))
        idx <- err[, 'sid'] > '0' ## skip subjects not genotyped
        err[idx, ]
    })
    dsg[err_idx] <- NA
    dsb[err_idx] <- NA
    
    list(
        gmp=gmp, sid=sid, vid=vid,
        dsg=dsg, dsb=dsb, mde=mde)
}

## load gene list
main <- function()
{
    ped <- ped.load()
    gen <- get.gen()

    ## process the sequence
    ret <- with(gen, mapply(function(ch, b1, b2, gn)
    {
        ## fetch VCF
        vcf <- get.seq(ch, b1, b2, gn)

        seq <- proc.seq(vcf, ped)
        seq <- within(
            seq,
        {
            chr <- ch
            bp1 <- b1
            bp2 <- b2
            name <- gn
        })
        seq
    }, chr, bp1, bp2, gen, SIMPLIFY = F))

    ## return
    names(ret) <- gen$gen
    ret
}

check.ME <- function(vcf, ped)
{
    ## pick out complete 
    ped <- ped.nuclear(ped)

    ## uniqne subject id
    si <- unique(unlist(ped), nmax = nrow(ped) * 2.5)
    si <- si[si > 0]
    
    ## pick out genotypes
    gt <- geno(vcf)$GT[, as.character(si)]

    ## pick out variant id
    gi <- rownames(gt)
    gt <- unname(gt)
        
    ## replace subject id with indices of genotype matrix
    pt <- apply(ped, 2, match, si)
    
    check.fam <- function(g)
    {
        g <- lapply(strsplit(g, '|'), `[`, c(1,3))

        ## genotype of the child, father, mather
        gc <- g[[1]]
        gf <- g[[2]]
        gm <- g[[3]]
        if(NA %in% gf)
            r <- any(gc %in% gm)
        else if(NA %in% gm)
            r <- any(gc %in% gf)
        else
            r <- any(
                gc[1] %in% gf && gc[2] %in% gm,
                gc[1] %in% gm && gc[2] %in% gf)
        r
    }

    ## go through families
    rt <- apply(pt, 1, function(fam)
    {
        apply(gt[, fam], 1, check.fam)
    })
    rt <- which(!rt, arr.ind = T, useNames = F)

    ## error genotype
    ge <- t(apply(rt, 1, function(r) gt[r[1], pt[r[2],]]))
    colnames(ge) <- c('igv', 'pgv', 'mgv')

    ## variant and subject id
    rt <- data.frame(
        vid=gi[rt[,1]],
        ped[rt[,2],],
        ge,
        stringsAsFactors = F)
    rownames(rt) <- NULL
    rt
}
