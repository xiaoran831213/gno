#xb# * -------- generics -------- * ##
## pick or show subjects
sbj <- function(dsg, ...) UseMethod("sbj")

## pick or list genomic variants
gvr <- function(dsg, ...) UseMethod("gvr")

## number of subjects
nsb <- function(dsg, ...) UseMethod("nsb")

## number of genomic variants
ngv <- function(dsg, ...) UseMethod("ngv")

## get genomic map
map <- function(dsg, ...) UseMethod("map")

## calculate minor allele frequency
maf <- function(dsg, ...) UseMethod("maf")

## genotype statistics
stt <- function(dsg, ...) UseMethod("stt")

## make missing values
mkMis <- function(dsg, ...) UseMethod("mkMis")

## fix minor allele frequency
fxMAF <- function(dsg, ...) UseMethod("fxMAf")

## remove degenerated variant
rmDgr <- function(dsg, ...) UseMethod("rmDgr")

## impute NAs
impute <- function(dsg, ...) UseMethod("impute")
## * -------- (end) -------- * ##

## constructor
dosage <- function(x = NULL, chr = NULL, bp1= NULL, bp2 = NULL, map = NULL, sbj = NULL, gmx = NULL)
{
    if(is.null(x))
    {
        x <- list()
    }
    else if(is.character(x))            # read from file
    {
        x <- readRDS(x)
    }

    ## overwrite members
    if(!is.null(chr))
        x$chr <- chr
    if(!is.null(bp1))
        x$bp1 <- bp1
    if(!is.null(bp2))
        x$bp2 <- bp2
    if(!is.null(map))
        x$map <- map
    if(!is.null(sbj))
        x$sbj <- sbj
    if(!is.null(gmx))
        x$gmx <- gmx

    x$.vw <- I
    structure(x, class = c('dosage', 'list'))
}

## number of subjects
nsb.dosage <- function(dsg, ...) length(dimnames(dsg$gmx)$sbj)

## number of genomic variants
ngv.dosage <- function(dsg, ...) length(dimnames(dsg$gmx)$gvr)

## list or pick out subject
sbj.dosage <- function(dsg, who = NULL)
{
    if(is.null(who))
        with(dsg, colnames(gmx))
    else
    {
        
        who <- unlist(who)
        if(is.character(who))
            within(dsg, gmx <- gmx[, match(who, colnames(gmx)), drop = F])
        else
            within(dsg, gmx <- gmx[, who, drop = F])
    }
}

## list or pick out genomic variants
gvr.dosage <- function(dsg, who = NULL)
{
    if(is.null(who))
        with(dsg, rownames(gmx))
    else
    {
        who <- unlist(who)
        if(is.character(who))
            within(dsg, gmx <- gmx[match(who, rownames(gmx)), ,drop = F])
        else
            within(dsg, gmx <- gmx[who, , drop = F])
    }
}

## * ---------- extend default R generics ---------- * ##
length.dosage <- function(dsg, ...) length(dsg$gmx)

dim.dosage <- function(dsg, ...) dim(with(dsg, .vw(gmx)))

t.dosage <- function(dsg, ...)
{
    ## change the view point of the data
    dsg[['.vw']] <- ifelse(identical(dsg[['.vw']], I), t, I)
    dsg
}

as.matrix.dosage <- function(dsg, ...)
{
    with(dsg, .vw(gmx))
}

as.vector.dosage <- function(dsg, ...)
{
    as.integer(with(dsg, .vw(gmx)))
}

`[.dosage` <- function(dsg, i = NULL, j = NULL, ..., drop = FALSE)
{
    d <- drop
    v <- with(dsg, .vw(gmx))

    if(is.null(i))
    {
        if(is.null(j))
            v[, , drop = d]
        else
            v[, j, drop = d]
    }
    else
    {
        if(is.null(j))
            v[i, , drop = d]
        else
            v[i, j, drop = d]
    }
}

str.dosage <- function(dsg, ...)
{
    with(
        dsg,
    {
        map = sprintf(
        "map: %s ~ %s",
        with(head(map, 1), sprintf('%s:%d-%d', chr, bp1, bp2)),
        with(tail(map, 1), sprintf('%s:%d-%d', chr, bp1, bp2)))

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
        paste(map, sbj, gvr, dim, "\n", sep = "\n")        
    })
}

print.dosage <- function(dsg, ...)
{
    with(
        dsg,
    {
        if(nrow(map) < 20L)
            print.data.frame(map)
        else
        {
            dfr <- rbind(
                as.matrix(head(map)),
                rep('.', ncol(map)),
                as.matrix(tail(map)))
            print(dfr, right = T, quote = F)
        }
    })
}

is.na.dosage <- function(dsg, ...) is.na(dsg$gmx)

## * ----------        (end extend)        ---------- * ##

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
    s <- stt(gmx);
    
    ## find variants not degenerated
    i <- which(pmax.int(s[,'n0'], s[,'n1'], s[,'n2']) < ncol(gmx) - s[,'nn'])
    
    gmx[i, , drop = F]
}
rmDgr.dosage <- function(dsg)
{
    within(dsg, gmx <- rmDgr.matrix(gmx))
}

## fix variants whoes MAF is greater than 0.5 by flipping their coding
fxMAf.vector <- function(gvt)
{
    if(maf.vector(gvt) > 0.5) 2L - gvt else gvt
}

fxMAf.matrix <- function(gmx)
{
    # mask wrong dosage values.
    idx <- which(maf.matrix(gmx) > 0.5)
    gmx[idx, ] = 2L - gmx[idx, ]
    gmx
}

fxMAf.dosage <- function(dsg)
{
    within(dsg, gmx <- fxMAF.matrix(gmx))
}

## guess missing values based on the frequency of know values
impute.vector <- function(gvt)
{
    ## calculate genotype count
    cnt<-stt.vector(gvt)

    if(cnt['nn'] > 0L)
        gvt[is.na(gvt)] <-
            sample.int(3L, cnt[4L], T, cnt[1L:3L]) - 1L
    gvt
}

## guess missing values based on the frequency of know values
impute.matrix <- function(gmx)
{
    ## calculate genotype statistics
    cnt<-stt.matrix(gmx);
    
    # guess missings for a variant, maintain type frequency
    for(i in which(cnt[, 'nn'] > 0L))
    {
        gmx[i, is.na(gmx[i,])] <-
            sample.int(3L, cnt[i, 4L], T, cnt[i, 1L:3L]) - 1L
    }
    gmx
}

impute.dosage <- function(dsg)
{
    within(dsg, gmx <- impute.matrix(gmx))
}

map.dosage <- function(dsg)
{
    with(dsg, map)
}

plot.dosage <- function(dsg, ...)
{
    x <- dsg$map$bp1

    # dummy plot to lay the background
    plot(x=0, y=0, xlim = range(x), ylim=c(0, 2))
    apply(dsg$gmx, 1L, function(y)
    {
        lines(x, y, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", ...)
    })
}
