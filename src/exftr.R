#!/usr/bin/env Rscript
source('src/mde.R')
source('src/vcf.R')
source('src/gls.R')

## process one sequence of vcf, return dosage data
proc.seq <- function(vcf, ped = NULL)
{
    ## subject id must be in pedigree
    ## check mendilian errors if pedigree is available,
    if(!is.null(ped))
    {
        sid <- sort(intersect(samples(header(vcf)), rownames(ped)))
        vcf <- vcf[, sid]
        ped <- ped.clr(ped[sid, ])
        mde <- MDE(vcf, ped)
    }
    
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

## load feature table
main <- function(
    wgs, dst = '.', ped = NULL, tab = NULL, wnd = 5000L, ovr = F, ...)
{
    ## handle destination directory
    if(is.null(dst))
        dst <- '.'
    dir.create(dst, F, T)

    if(!is.null(ped))
        ped <- ped.load(ped)
    
    if(is.null(tab))
        tab <- gls.load()
    
    ## process the sequence
    ret <- apply(tab, 1L, function(.a)
    {
        .e <- environment()
        mapply(function(nm, ar, cl)
        {
            assign(nm, as(ar, cl), envir=.e)
        }, names(.a), .a, lapply(tab, class))
        
        out <- file.path(dst, sprintf("%s.rds", ssn))
        seq <- sprintf('%2d:%9d - %9d', chr, bp1, bp2)
        cat(sprintf('xt: %s -> %s, ', seq, out))

        ## cat will return NULL, which means OK.
        if(file.exists(out) && !ovr)
            return(cat("exists\n")) 

        ## fetch and process VCF
        vcf <- try(get.seq(wgs, chr, bp1 - wnd, bp2 + wnd))
        if(inherits(vcf, 'try-error'))
            return(message(vcf))
        if(nrow(vcf) == 0L)
            return(cat('empty\n'))

        ## turn vcf to dsg, check ME if ped is available
        dsg <- try(proc.seq(vcf, ped))
        if(inherits(dsg, 'try-error'))
            return(message(dsg))

        ## append additional information
        dsg <- dosage(c(dsg, .a))
        
        ## return NULL means OK
        saveRDS(dsg, out)
        return(cat('done\n'))
    })
}

.argparser <- function()
{
    library(argparser)
    p <- arg_parser('genome extractor.')
    p <- add_argument(
        p, '--sub',
        help = paste(
            "R boolean syntex to subset features from the table:\n",
            "\tchr: chromosome\n",
            "\tbp1: starting basepair\n",
            "\tbp2: ending basepair\n",
            "\tgid: gene id\n",
            "\tsmb: symbol (short name)\n",
            "\tfnm: full name\n",
            "\tifo: other informations"))
    p <- add_argument(
        p, '--rng',
        help = paste(
            "range of further selection from the subset,",
            "using format 'start,end'."),
            default=',')
    p <- add_argument(
        p, '--dst', help = 'directory to  place extracted features.',
        default = '.')
    p <- add_argument(
        p, '--wgs', help = 'directory to look for WGS in VCF format.',
        default = ".")
    p <- add_argument(
        p, '--ped', help = 'path to the pedigree file.')
    p <- add_argument(
        p, '--wnd', help = 'genome extraction window',
        default = 5000L)
    p <- add_argument(
        p, '--dry', help = 'only show the features to be extrated.',
        default = FALSE)
    p <- add_argument(
        p, '--ovr', help = 'overwrite existing output.',
        default = FALSE)
    p
}

cml <- function(...)
{
    p <- .argparser()
    argv <- c(commandArgs(trailingOnly = TRUE), ...)
    print(argv)
    if(length(argv) < 1)
    {
        print(p) ## the help
        return
    }
    opt <- parse_args(p, argv)
    
    tab <- with(
        opt,
    {
        ## get feature table first
        tab <- gls.load()
        tab <- data.frame(
            ssn = rownames(tab), tab,
            stringsAsFactors = F)
        
        ## write subsetting script
        if(!is.na(sub)) tab <- eval(parse(
            text = sprintf('subset(tab, %s)', sub)))

        ## index range
        rng <- as.integer(unlist(strsplit(rng, ',')))
        if(is.na(rng[1]) || rng[1] == "")
            rng[1] = 1L
        if(is.na(rng[2]) || rng[2] == "")
            rng[2] = nrow(tab)
        tab <- tab[rng[1]:rng[2], ]
        tab
    })

    if(opt$dry)
        return(tab)

    opt <- lapply(opt, function(o) if(is.na(o)) NULL else o)
    opt$tab <- tab
    do.call(main, args = opt)
}

cml()
## cml('--rng', '10,20', '--wgs', 'wgs', '--dst', 'd1', '--ovr', 'T', '--ped', 'dat/all.ped')
