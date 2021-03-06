source('src/ped.R')
source('src/mde.R')
source('src/vcf.R')
source('src/gls.R')

## check mendilian error for vcf format
MDE <- function(vcf, ped)
{
    suppressPackageStartupMessages(library(VariantAnnotation, warn.conflicts = F, quietly = T))

    ## pick out families
    ped <- ped.nuclear(ped)

    ## uniqne subject id
    si <- unique(unlist(ped))
    si <- si[si > 0]

    ## pick out genotypes
    gt <- geno(vcf)$GT[, as.character(si)]

    ## pick out variant id
    gi <- rownames(gt)
    gt <- unname(gt)
        
    ## replace subject id with indices of genotype matrix
    pt <- apply(ped, 2L, match, si)
    
    check_fam <- function(g)
    {
        g <- lapply(strsplit(g, '|'), `[`, c(1L, 3L))

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
    err <- apply(pt, 1, function(fam)
    {
        apply(gt[, fam], 1, check_fam)
    })
    err <- which(!err, arr.ind = T, useNames = F)

    ## error genotype
    ge <- t(apply(err, 1, function(r) gt[r[1], pt[r[2],]]))
    colnames(ge) <- c('igv', 'pgv', 'mgv')

    ## variant id, and subject id
    err <- data.frame(
        vid=gi[err[,1]],
        ped[err[,2],],
        ge,
        stringsAsFactors = F)
    rownames(err) <- NULL
    err

    ## collection error genotype index
    idx <- with(
        err,
    {
        idx <- cbind(vid=vid, sid=c(iid, pid, mid))
        idx[idx[, 'sid'] > '0', ]       # skip null ID
    })

    ## return mendilian error structure
    structure(list(err=err, idx=idx), class = c('MDE', 'list'))
}

## process one sequence of vcf, return dosage data
proc.seq <- function(vcf, ped)
{
    ## subject id must be in pedigree
    ## check mendilian errors if pedigree is available,
    sid <- sort(intersect(samples(header(vcf)), rownames(ped)))
    vcf <- vcf[, sid]
    ped <- ped.clr(ped[sid, ])
    mde <- MDE(vcf, ped)
    
    ## make dosage data
    dsg <- vcf2dsg(vcf)

    ## if ME check is done, set the erronous genotype to NA
    if(exists('mde'))
    {
        dsg$mde <- mde
        dsg$gmx[mde$idx] <- NA
    }
    dsg
}
