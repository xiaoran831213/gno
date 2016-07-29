#!/usr/bin/env Rscript

## write range in bcftools' syntax
.vcf.rng <- function(chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 0L)
{
    if(is.null(chr) || is.null(bp1))
        return("")

    cmp <- as.character(1L:22L)
    names(cmp) <- cmp
    cmp <- c(cmp, `23`='X', `24`='Y', `25`='MT', `26`='XY')
    chr <- cmp[chr]

    if(is.null(bp2))
        bp2 <- bp1
    bp1 <- bp1 - wnd
    bp2 <- bp2 + wnd

    rng <- sprintf('%s:%s-%s', chr, bp1, bp2)
    rng <- do.call(paste, c(as.list(rng), sep = ','))
    rng <- paste('-r', rng)
    rng
}

## get genome map via bcftools
.vcf.map <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 0L, stderr = FALSE)
{
    rng <- .vcf.rng(chr, bp1, bp2, wnd)
    cmd <- "bcftools query -f '%CHROM %POS %ID %REF %ALT\\n'"
    cmd <- paste(cmd, rng, vcf)
    if(!stderr)
        cmd <- paste(cmd, '2>/dev/null')
    
    sed<-"sed \'s/\\(^X\\)/23/; s/\\(^Y\\)/24/\'"
    cmd<-paste(cmd, sed, sep= "|")
    pip<-pipe(cmd, "r")
    txt <- readLines(pip)
    close(pip)

    ## create map
    i0 <- integer()
    c0 <- character()
    map <- data.frame(CHR = i0, POS = i0, UID = c0, REF = c0, ALT = c0)
    hdr <- names(map)

    ## empty map
    if(length(txt) == 0L)
        return(map)

    map <- read.table(text = txt, as.is = T)
    names(map) <- hdr

    map
}

## get subjects via bcftools
.vcf.sbj <- function(vcf)
{
    cmd <- paste("bcftools query -l", vcf)
    pip <- pipe(cmd, "r");
    sbj <- scan(file = pip, what = " ", quiet = T)
    close(pip);
    sbj
}

.vcf.gvt <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 0L, stderr = FALSE)
{
    rng <- .vcf.rng(chr, bp1, bp2, wnd)
    cmd<-"bcftools query -f '[%GT ]\\n'";
    cmd<-paste(cmd, rng);
    cmd<-paste(cmd, vcf);
    if(!stderr)
        cmd <- paste(cmd, '2>/dev/null')
    
    ## construct sed filter command
    .sd <- paste(
        "s/0[|/]0/0/g",
        "s/[1-9][|/]0/1/g",
        "s/0[|/][1-9]/1/g",
        "s/[1-9][|/][1-9]/2/g",
        "s/\\.[|/]./3/g",
        "s/.[|/]\\./3/g",
        sep = ';')
    .sd <- paste('sed \'', .sd, '\'', sep = '')

    ## the final command
    cmd<-paste(cmd, .sd, sep="|");
    pip<-pipe(cmd, "r");
    gvt<-scan(pip, what=0L, na.strings = '3', quiet=T)
    close(pip)
    gvt
}

## read vcf data via bcftools, convert it to dosage format in R
readVCF <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 0L, ...)
{
    if(!file.exists(vcf))
        stop(vcf, ' not found.')
    
    ## read genome map
    map <- .vcf.map(vcf, chr, bp1, bp2, wnd)
    if(nrow(map) == 0L)
        stop('NG')
    gvr <- rownames(map)
    
    ## -------- get subject id from VCF header -------- #
    sbj <- .vcf.sbj(vcf)
    
    ## -------- get genotype matrix -------- #
    dat <- .vcf.gvt(vcf, chr, bp1, bp2, wnd)
    gmx<-matrix(
        .vcf.gvt(vcf, chr, bp1, bp2, wnd),
        nrow=length(gvr), ncol=length(sbj), byrow=T, dimnames=list(gvr=gvr, sbj=sbj))

    ## the output
    within(
        list(...),
    {
        chr <- as.integer(chr)
        bp1 <- as.integer(bp1)
        bp2 <- as.integer(bp2)
        wnd <- as.integer(wnd)
        map <- map
        gmx <- gmx
        vcf <- vcf
    })
}

## the command line parser
.argparser <- function()
{
    library(argparser, quietly = T)
    p <- arg_parser('VCF to DSG converter, RDS version.')
    aa <- add_argument

    p <- aa(p, 'vcf', help = 'the VCF input.')
    p <- aa(p, 'rds', help = 'the RDS output.')
    p <- aa(p, '--rgn', help = "common delimeted regions (chr:bp1-bp2[;var=val]).")
    p <- aa(p, '--RGN', help = "target region file.")
    p <- aa(p, '--wnd', help = 'flanking window size', default = 0L)
    p <- aa(p, '--dsc', help = 'region description, [var=val][;var=val].')
    p <- aa(p, '--ovr', help = 'overwrite existing output.', default = FALSE)
    p <- aa(p, '--vbs', help = 'verbose', default = FALSE)
    p
}

## VCF to DSG convertor, the output is store in a RDS file.
main <- function(...)
{
    p <- .argparser()
    argv <- c(commandArgs(trailingOnly = TRUE), ...)
    
    if(length(argv) < 1)
    {
        print(p)
        return()
    }

    ## parse the arguments
    opt <- parse_args(p, argv)
    ext <- FALSE

    ss <- function(x,s) unlist(strsplit(x, s), use.names = F)
    ## initialization at command line
    opt <- within(opt,
    {
        ## parese genomic region
        if(!is.na(rgn))
        {
            ## break down regions and descriptions
            rgn <- strsplit(rgn, ',')
            ## the 1st item is the loacation
            rgn <- lapply(rgn, function(r)
            {
                m <- regexec('^([0-9XY]*):([0-9]*)(|-([0-9]*))$', r)
                r <- unlist(regmatches(r, m))
                r <- c(r[2L], r[3L], r[5L])
                r
            })
            rgn <- as.data.frame(
                do.call(rbind, rgn), stringsAsFactors = F)
        }
        else
            rgn <- integer()

        ## regions from file
        if(!is.na(RGN))
        {
            RGN <- read.table(RGN, as.is = T)
            rgn <- rbind(rgn, RGN[,1:3])
            for(i in 4:length(RGN))
                assign(names(RGN)[i], RGN[[i]])
        }
        
        ## map region locations to integers
        cmp <- 1L:26L
        cmp[as.character(cmp)] <- cmp
        cmp <- c(cmp, X=23L, Y=24, MT=25, XY=26)
        chr <- as.integer(cmp[rgn[,1]])
        bp1 <- as.integer(rgn[,2])
        bp2 <- as.integer(rgn[,3])
        rm(cmp, rgn, RGN)

        ## verbosity option
        if(vbs)
            print(opt)

        ## check overwrite
        if(!ovr && file.exists(rds))
            ext <<- TRUE
        
        ## region descriptions
        if(!is.na(dsc))
        {
            dsc <- unlist(strsplit(dsc, ','), use.names = F)
            dsc <- strsplit(dsc, '=')
            val <- sapply(dsc, `[`, 2L)
            dsc <- sapply(dsc, `[`, 1L)
            for(i in 1:length(dsc))
                assign(dsc[i], val[i])
            rm(val, i, dsc)
        }
        
        rm(vbs, ovr, help, opts)
    })

    ## verbosity option
    ## print out message
    msg <- with(opt, paste('xt:', chr, bp1, bp2, vcf, rds, sep = '\t'))
    
    ## skip extraction if the output exists and overwite is turned off
    if(ext)
    {
        msg <- paste(msg, 'EX', sep = "\t")
        ret <- NULL
    }
    else
    {
        ## drop NULL, so they won't be stored in the output
        opt <- opt[!sapply(opt, is.null)]
        ret <- try(do.call(readVCF, opt), silent = T)
        
        if(class(ret) == 'try-error')
        {
            ret <- attr(ret, 'condition')['message']
            msg <- paste(msg, ret, sep = '\t')
        }
        else
        {
            saveRDS(ret, opt$rds)
            msg <- paste(msg, 'OK', sep = '\t')
        }
    }
    msg <- paste(msg, collapse = '\n')
    cat(msg)
    cat('\n')
    invisible(ret)
}

test <- function()
{
    b <- paste(
        '--rgn 1:1334909-1337426,2:2434909-2437426',
        ##'--RGN g1.txt',
        '--dsc SMB=LOC148413,FNM=abcdefg',
        '1/all.vcf.gz G0032.rds --ovr T')
    b <- unlist(strsplit(b, ' '))
}

main()
##cml()
## cml('--rng', '10,20', '--wgs', 'wgs', '--dst', 'd1', '--ovr', 'T', '--ped', 'dat/all.ped')
