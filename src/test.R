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
    dst <- ("dat/dsg")
    
    ## process the sequence
    ret <- with(
        gen,
    {
        mapply(function(ch, b1, b2, gi, n1, n2)
        {
            dst <- file.path('dat/dsg', sprintf("%08X.rds", gi))
            whr <- sprintf('%2d:%9d-%9d', ch, b1, b2)
            if(file.exists(dst))
            {
                cat(sprintf('%s -> %s: exists\n', whr, dst))
                return(NULL)
            }

            ## fetch and process VCF
            cat(sprintf('%s -> %s: created\n', whr, dst))

            vcf <- try(get.seq('wgs', ch, b1 - 5000L, b2 + 5000L))
            if(inherits(vcf, 'try-error'))
                return(NULL)

            dsg <- try(proc.seq(vcf, ped))
            if(inherits(dsg, 'try-error'))
                return(NULL)
            
            dsg <- within(dsg, {gid=gi; sym=n1; dsc=n2})
            saveRDS(dsg, dst)
            rm(dsg)
            
            NULL
        }, chr, bp1, bp2, gid, sym, fnm, SIMPLIFY = F)
    })

    ## return
    names(ret) <- gen$gen
    ret
}

