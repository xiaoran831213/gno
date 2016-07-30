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
        x <- within(readRDS(x), rds <- x)
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

    ## the view of the data
    x$view <- I                          
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

dim.dosage <- function(dsg, ...) with(dsg, dim(view(gmx))) #

t.dosage <- function(dsg, ...)
{
    ## change the view point of the data
    within(dsg, view <- if(identical(view, I)) t else I)
}

as.matrix.dosage <- function(dsg, ...)
{
    with(dsg, view(gmx))
}

as.vector.dosage <- function(dsg, ...)
{
    as.integer(with(dsg, view(gmx)))
}

`[.dosage` <- function(dsg, i = NULL, j = NULL, ..., drop = FALSE)
{
    d <- drop

    ## consider transposed genomic matrix
    if(identical(dsg$view, t))
    {
        k <- j
        j <- i
        i <- k
    }

    ## pick out variants and subjects
    dsg <- within(dsg,
    {
        if(is.null(i))
        {
            if(is.null(j))
                gmx <- gmx[ ,  , drop = d]
            else
                gmx <- gmx[ , j, drop = d]
        }
        else
        {
            if(is.null(j))
                gmx <- gmx[i,  , drop = d]
            else
                gmx <- gmx[i, j, drop = d]
        }
    })
    dsg
}

str.dosage <- function(x, ...)
{
    with(
        x,
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

print.dosage <- function(x, ...)
{
    with(
        x,
    {
        if(nrow(map) < 20)
            print.data.frame(map)
        else
        {
            dfr <- rbind(
                as.matrix(head(map)),
                rep('.', ncol(map)),
                as.matrix(tail(map)))
            print(dfr, right = T, quote = F)
        }

        cat(paste(c('ngv', 'nsb'), c(ngv(x), nsb(x)), sep='='), '\n')
    })
}

is.na.dosage <- function(x, ...) is.na(x$gmx)
## * ----------        (end extend)        ---------- * ##

maf.dosage <- function(x)
{
    with(x, maf.matrix(gmx))
}

maf.matrix <- function(gmx)
{
    maf <- rowMeans(gmx, na.rm = T) / 2L
    pmin(maf, 1-maf)
}
maf.Matrix <- maf.matrix

maf.vector <- function(gvt)
{
    maf <- mean(gvt, na.rm = T) / 2L
    min(maf, 1-maf)
}

## randomly make missing values in genotype matrix
mkMis.matrix <- function(gmx, frq, val=NA)
{
    gmx[sample.int(length(gmx), length(gmx) * frq)] <- val
    gmx
}
mkMis.Matrix <- mkMis.matrix

## randomly make missing values in genotype matrix
mkMis.dosage<-function(x, frq, val=NA)
{
    within(x, gmx <- mkMis.matrix(gmx))
}

.stt.itm <- list('N0', 'N1', 'N2', 'NN')
stt.vector <- function(gvt, nm = TRUE)
{
    c(sum(gvt == 0L, na.rm = T),     # 1.Homo-major
      sum(gvt == 1L, na.rm = T),     # 2.Hete
      sum(gvt == 2L, na.rm = T),     # 3.Homo-minor
      sum(is.na(gvt)))               # 4.missing
}

## get genotype counts, only works for biallelic dosage data for now.
stt.matrix <- function(gmx)
{
    ## statistics
    out <- matrix(0L, nrow(gmx), 4L)
    for(i in 1L:nrow(gmx))
    {
        out[i,2] = sum(gmx[i,] == 1L, na.rm = T)     # 2.Hete
        out[i,3] = sum(gmx[i,] == 2L, na.rm = T)     # 3.Homo-minor
        out[i,4] = sum(is.na(gmx[i,]))               # 4.missing
    }
    out[,1] <- nrow(gmx) - rowSums(out[,2:4])
    as.data.frame(
        out,
        col.names = c('N0', 'N1', 'N2', 'NN'))
}
stt.Matrix <- function(gmx)
{
    ## statistics
    ## avoid copying sparse matrix to the dense one
    out <- data.frame(
        N1 = rowSums(gmx == 1L, T),
        N2 = rowSums(gmx == 2L, T),
        NN = rowSums(is.na(gmx)))
    out <- cbind(N0 = ncol(gmx) - rowSums(out), out)
    out
}

stt.dosage <- function(x)
{
    stt.matrix(x$gmx)
}

rmDgr.matrix <- function(gmx, ret = c(gmx, msk, idx), inv = FALSE)
{
    ret <- match.arg(ret, choices = c('gmx', 'msk', 'idx'))
    ## calculate genotype statistics
    s <- stt(gmx);

    ## find variants not degenerated
    m <- pmax.int(s[,1], s[,2], s[,3]) < ncol(gmx) - s[,4]
    if(inv) i <- !m

    switch(ret, gmx=gmx[m, , drop = F], msk=m, idx=which(m))
}
rmDgr.Matrix <- rmDgr.matrix

rmDgr.dosage <- function(dsg, ret = c('gmx', 'msk', 'idx'), inv = FALSE)
{
    within(dsg,
    {
        i <- rmDgr.matrix(gmx, ret = 'idx')
        gmx <- gmx[i, , drop = FALSE]
        map <- map[i, , drop = FALSE]
    })
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

fxMAf.dosage <- function(x)
{
    within(x, gmx <- fxMAF.matrix(gmx))
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

impute.matrix <- function(gmx)
{
    ## calculate genotype statistics
    cnt <- stt(gmx)
    
    # guess missings for a variant, maintain type frequency
    for(i in which(cnt[, 4] > 0L))
    {
        gmx[i, which(is.na(gmx[i,]))] <- sample.int(3L, cnt[i, 4], T, cnt[i, 1:3]) - 1L
    }
    gmx
}
impute.Matrix <- impute.matrix

impute.dosage <- function(x)
{
    within(x, gmx <- impute.matrix(gmx))
}

map.dosage <- function(x)
{
    with(x, map)
}

plot.dosage <- function(x, ...)
{
    x <- x$map$bp1

    # dummy plot to lay the background
    plot(x=0, y=0, xlim = range(x), ylim=c(0, 2))
    apply(x$gmx, 1L, function(y)
    {
        lines(x, y, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", ...)
    })
}
