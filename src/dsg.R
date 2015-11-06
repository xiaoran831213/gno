source('src/utl.R')
source('src/generics.R')

## generics

## pick or show subjects of a data
sbj <- function(gno, ...) UseMethod("sbj")

gmp <- function(gno, ...) UseMethod("gmp")

## calculate minor allele frequency
maf <- function(gno, ...) UseMethod("maf")

## genotype statistics
stt <- function(gno, ...) UseMethod("stt")
mkMis <- function(gno, ...) UseMethod("mkMis")
fxMaf <- function(gno, ...) UseMethod("fxMaf")
rmDgr <- function(gno, ...) UseMethod("rmDgr")
imp <- function(gno, ...) UseMethod("imp")

## constructor
dosage <- function(x = NULL, gmp = NULL, gmx = NULL)
{
    if(is.null(x))
    {
        x <- list(gmp = gmp, gmx = gmx)
    }
    structure(x, class = c('dosage', 'list'))
}

str.dosage <- function(dsg, ...)
{
    with(
        dsg,
    {
        gmp = sprintf(
        "gmp: %s ~ %s",
        with(head(gmp, 1), sprintf('%s:%d-%d', chr, bp1, bp2)),
        with(tail(gmp, 1), sprintf('%s:%d-%d', chr, bp1, bp2)))

        sbj = sprintf(
        "sbj: %s... %s",
        paste(head(colnames(gmx)), collapse=" "),
        paste(tail(colnames(gmx)), collapse=" "))
        
        gvr = sprintf(
        "gvr: %s... %s",
        paste(head(rownames(gmx)), collapse=" "),
        paste(tail(rownames(gmx)), collapse=" "))
        
        dim = sprintf(
        "num.gvr = %d, num.sbj = %d",
        nrow(gmx), ncol(gmx))
        paste(gmp, sbj, gvr, dim, "\n", sep = "\n")        
    })
    
}

print.dosage <- function(dsg, ...)
{
    with(
        dsg,
    {
        if(nrow(gmp) < 20L)
            print.data.frame(gmp)
        else
        {
            dfr <- rbind(
                as.matrix(head(gmp)),
                rep('.', ncol(gmp)),
                as.matrix(tail(gmp)))
            print(dfr, right = T, quote = F)
        }
    })
}

maf.dosage <- function(dsg)
{
    with(dsg, maf.matrix(gmx))
}

maf.matrix <- function(gmx)
{
    maf <- rowMeans(gmx, na.rm = T) / 2L
    pmin(maf, 1-maf)
}

maf.vector <- function(gvt)
{
    maf <- mean(gvt, na.rm = T) / 2L
    min(maf, 1-maf)
}

## randomly make missing values in genotype matrix
mkMis.matrix<-function(gmx, frq, val=NA)
{
    gmx[sample.int(length(gmx), length(gmx) * frq)] <- val
    gmx
}

## randomly make missing values in genotype matrix
mkMis.dosage<-function(dsg, frq, val=NA)
{
    within(dsg, gmx <- mkMis.matrix(gmx))
}

.stt.itm <- list('N0', 'N1', 'N2', 'NN')
stt.vector <- function(gvt, nm = TRUE)
{
    c(n0 = sum(gvt == 0L, na.rm = T),     # 1.Homo-major
      n1 = sum(gvt == 1L, na.rm = T),     # 2.Hete
      n2 = sum(gvt == 2L, na.rm = T),     # 3.Homo-minor
      nn = sum(is.na(gvt)))               # 4.missing
}

## get genotype counts, only works for biallelic dosage data for now.
stt.matrix <- function(gmx)
{
    ## variant id
    dnm <- dimnames(gmx)
    vid <- dimnames(gmx)[[1]]

    ## statistics
    ret <- cbind(
        t(apply(gmx, 1L, stt.vector)),
        maf = maf(gmx))
    ##dimnames(ret) <- list(vid=vid, itm=.stt.itm)
    ret
}
stt.dosage <- function(dsg)
{
    stt.matrix(dsg$gmx)
}

## remove degenerated variants
rmDgr.matrix <- function(gmx)
{
    ## calculate genotype statistics
    cnt <- stt(gmx);
    
    ## find degenerated variants.
    idx <- which(cnt[, 4L] - pmax.int(cnt[, 1L:3L]) > 0)
    
    gmx[idx, , drop = F]
}
rmDgr.dosage <- function(dsg)
{
    within(dsg, gmx <- rmDgr.matrix(gmx))
}

## fix variants whoes MAF is greater than 0.5 by flipping their coding
fxMaf.vector <- function(gvt)
{
    if(maf.vector(gvt) > 0.5) 2L - gvt else gvt
}

fxMaf.matrix <- function(gmx)
{
    # mask wrong dosage values.
    idx <- which(maf.matrix(gmx) > 0.5)
    gmx[idx, ] = 2L - gmx[idx, ]
    gmx
}

fxMaf.dosage <- function(dsg)
{
    within(dsg, gmx <- fxMaf.matrix(gmx))
}

## guess missing values based on the frequency of know values
imp.vector <- function(gvt)
{
    ## calculate genotype count
    cnt<-stt.vector(gvt)

    if(cnt[4L] > 0L)
        gvt[is.na(gvt)] <-
            sample.int(3L, cnt[4L], T, cnt[1L:3L]) - 1L
    gvt
}

## guess missing values based on the frequency of know values
imp.matrix <- function(gmx)
{
    ## calculate genotype statistics
    cnt<-stt.matrix(gmx);

    # guess missings for a variant, maintain type frequency
    for(i in which(cnt[, 'NN'] > 0L))
    {
        gmx[i, is.na(gmx[i,])] <-
            sample.int(3L, cnt[i, 4L], T, cnt[i, 1L:3L]) - 1L
    }
    gmx
}

imp.dosage <- function(dsg)
{
    within(dsg, gmx <- imp.matrix(gmx))
}

## pick out subject from genotype data
sbj.dosage <- function(dsg, who = NULL)
{
    if(is.null(who))
    {
        with(dsg, colnames(gmx))
    }
    else
    {
        who <- unlist(who)
        if(is.character(who))
            within(dsg, gmx <- gmx[, match(who, colnames(gmx)), drop = F])
        else
            within(dsg, gmx <- gmx[, who, drop = F])
    }
}

gvr.dosage <- function(dsg, who = NULL)
{
    if(is.null(who))
    {
        with(dsg, rownames(gmx))
    }
    else
    {
        who <- unlist(who)
        if(is.character(who))
            within(dsg, gmx <- gmx[match(who, rownames(gmx)), ,drop = F])
        else
            within(dsg, gmx <- gmx[who, , drop = F])
    }
}

gmp.dosage <- function(dsg)
{
    with(dsg, gmp)
}


plot.dosage <- function(dsg, ...)
{
    x <- dsg$gmp$bp1

    # dummy plot to lay the background
    plot(x=0, y=0, xlim = range(x), ylim=c(0, 2))
    apply(dsg$gmx, 1L, function(y)
    {
        lines(x, y, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", ...)
    })
}
