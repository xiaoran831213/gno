.vcf.rng <- function(chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 0L)
{
    if(is.null(chr))
        return("")

    if(chr == 23L)
        chr = 'X'
    if (chr == 24L)
        chr = 'Y'

    if(is.null(bp1))
        bp1 <- 0L
    else
        bp1 <- bp1 - wnd
    if(is.null(bp2))
        bp2 <- .Machine$integer.max
    else
        bp2 <- bp2 - wnd

    return(sprintf('-r %s:%s-%s', chr, bp1, bp2))
}

.vcf.map <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 5000L, stderr = FALSE)
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

    if(length(txt) == 0L)               # empty map
        map <- data.frame(
            CHR = integer(0), POS = integer(0), UID = character(0), REF = character(0),
            ALT = character(0), stringsAsFactors = FALSE)
    else
    {
        map <- read.table(text = txt, header = F, as.is = T)
        colnames(map) <- c("CHR", "POS", "UID", "REF", "ALT")
        
        rownames(map) <- sprintf('V%04X', 1L:nrow(map))
    }
    map
}

.vcf.sbj <- function(vcf)
{
    cmd <- paste("bcftools query -l", vcf)
    pip <- pipe(cmd, "r");
    sbj <- scan(file = pip, what = " ", quiet = T)
    close(pip);
    sbj
}

.vcf.gvt <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 5000L, stderr = FALSE)
{
    rng <- .vcf.rng(chr, bp1, bp2, wnd)
    cmd<-"bcftools query -f '[%GT ]\\n'";
    cmd<-paste(cmd, rng);
    cmd<-paste(cmd, vcf);
    if(!stderr)
        cmd <- paste(cmd, '2>/dev/null')
    ## construct sed filter command
    sed<-"sed \' s/0[|/]0/0/g; s/[1-9][|/]0/1/g; s/0[|/][1-9]/1/g; s/[1-9][|/][1-9]/2/g; s/\\.[|/]./3/g; s/.[|/]\\./3/g\'"
    
    ## the final command
    cmd<-paste(cmd, sed, sep="|");
    pip<-pipe(cmd, "r");
    gvt<-scan(pip, what=0L, na.strings = '3', quiet=T)
    close(pip)
    gvt
}

## read vcf data, convert it to dosage
vcf2dsg <- function(vcf, chr = NULL, bp1 = NULL, bp2 = NULL, wnd = 5000L, smb = NULL, fnm = NULL)
{
    ## read genome map
    map <- .vcf.map(vcf, chr, bp1, bp2)
    gvr <- rownames(map)
    
    ## -------- get subject id from VCF header -------- #
    sbj <- .vcf.sbj(vcf)
    
    ## -------- get genotype matrix -------- #
    gmx<-matrix(
        .vcf.gvt(vcf, chr, bp1, bp2, wnd),
        nrow=length(gvr), ncol=length(sbj), byrow=T, dimnames=list(gvr=gvr, sbj=sbj))
    
    ## the output
    within(
        list(),
    {
        chr <- as.integer(chr)
        bp1 <- as.integer(bp1)
        bp2 <- as.integer(bp2)
        wnd <- as.integer(wnd)
        map <- map
        gmx <- gmx
        vcf <- vcf
        smb <- smb
        fnm <- fnm
    })
}
