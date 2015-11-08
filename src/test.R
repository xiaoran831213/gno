#!/usr/bin/env Rscript
source('src/ped.R')
source('src/mde.R')
source('src/vcf.R')
source('src/gls.R')

## process one sequence of vcf, return dosage data
proc.seq <- function(vcf, ped)
{
    ## subject id must be in pedigree
    sid <- sort(intersect(samples(header(vcf)), rownames(ped)))
    vcf <- vcf[, sid]
    ped <- ped.clr(ped[sid, ])
    
    ## make dosage data
    dsg <- within(
        vcf2dsg(vcf),
    {
        ## check mendilian errors
        mde <- MDE(vcf, ped)

        ## set erronous genotype to NA
        gmx[mde$idx] <- NA
    })

    ## possible memeory leak, remove it.
    rm(vcf)
    dsg
}

## load gene list
main <- function()
{
    ped <- ped.load("dat/all.ped")
    gen <- gls.load()
    gen <- subset(gen, chr < 24)
    dst <- ("dat/dsg")
    
    ## process the sequence
    ret <- with(
        gen,
    {
        mapply(function(sn, ch, b1, b2, gi, n1, n2, ds)
        {
            dst <- file.path('dat/dsg', sprintf("%s.rds", sn))
            whr <- sprintf('%2d:%9d-%9d', ch, b1, b2)
            if(file.exists(dst))
            {
                cat(sprintf('%s -> %s: exists.\n', whr, dst))
                return(NULL)
            }

            ## fetch and process VCF
            cat(sprintf('%s -> %s: extracting', whr, dst))
            
            vcf <- try(get.seq('wgs', ch, b1 - 5000L, b2 + 5000L))
            if(inherits(vcf, 'try-error'))
                return(NULL)

            cat(', processing')
            dsg <- try(proc.seq(vcf, ped))
            if(inherits(dsg, 'try-error'))
                return(NULL)
            
            dsg <- within(dsg, {gid=gi; smb=n1; fnm=n2; ifo=ds})
            saveRDS(dsg, dst)

            cat(', done.\n')
            rm(dsg, vcf)
            NULL
        }, rownames(gen), chr, bp1, bp2, gid, smb, fnm, ifo, SIMPLIFY = F)
    })

    ## return
    names(ret) <- gen$gen
    ret
}

gno.ext <- function()
{
    library(argparser)
    p <- arg_parser('genome extractor.')
    p <- add_argument(
        p, '--sub',
        help = paste(
            "R boolean syntex to take a subset of features from the table:",
            "\tchr: chromosome",
            "\tbp1: starting basepair",
            "\tbp2: ending basepair",
            "\tgid: gene id",
            "\tsmb: symbol (short name)",
            "\tfnm: full name",
            "\tifo: other informations"),
        default = "TRUE")
    p <- add_argument(
        p, '--rng',
        help = paste(
            "range of subscripts to further select from the above subset,",
            "using format 'start,end'."),
            default='')
    p <- add_argument(
        p, '--dst', help = 'where to place extract features.',
        default = '.')
    p <- add_argument(
        p, '--wgs', help = 'where to look for WGS files, in VCF format.',
        default="")
    p <- add_argument(
        p, '--ped', help = 'path to the pedigree file.',
        default="")

    argv <- commandArgs(trailingOnly = TRUE)
    if(length(argv) > 0)
    {
        opt <- parse_args(p, argv)
        print(opt)
        
        with(
            opt,
        {
            syntex <- sprintf('subset(gls.load()[:], %s)', sub, )
            rng <- as.integer(unlist(strsplit(n.s, ',')))
            ftr <- gls.load()[rng, ]

        })
        cat('xt: success\n')
    }
}

cml.mix()
