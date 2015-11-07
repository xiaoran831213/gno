## load gene list
gls.load <- function(recache = FALSE)
{
    cache <- '.genome.features.rds'
    
    ## try cached gene list first
    if(file.exists(cache) && !recache)
    {
        ftr <- readRDS(cache)
        return(ftr)
    }
    
    ## download genomic features, save it to R binary
    src <- paste(
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq",
        "vertebrate_mammalian/Homo_sapiens",
        "latest_assembly_versions/GCF_000001405.31_GRCh38.p5",
        "GCF_000001405.31_GRCh38.p5_feature_table.txt.gz",
        sep = '/')
    tmp <- tempfile()
    download.file(src, tmp)
    ftr <- readLines(gzf <- gzfile(tmp))
    close(gzf)
    unlink(tmp)
    
    ## remove header comment character
    hdr <- unlist(strsplit(ftr[1], '\t'))
    hdr[1] <- 'feature'
    
    ## pick out genes and primary assemblies
    ftr <- ftr[grep('^gene\t.*\tPrimary Assembly\t', ftr)]
    
    ## create dtat table
    ftr <- read.table(
        text=ftr, col.names = hdr, nrows = length(ftr), fill = T,
        sep = '\t', as.is = T, quote = "")
    
    ## reformat, use short column names, remove some columns
    ftr <- with(
        ftr,
    {
        data.frame(
            chr = sapply(chromosome, chr2num, USE.NAMES = F),
            bp1 = as.integer(start),
            bp2 = as.integer(end),
            gid = GeneID,
            smb = symbol,
            fnm = name,
            ifo = sprintf(
                 'clz=%s,sq_typ=%s,attr=%s',
                 class, seq_type, attributes),
            stringsAsFactors = F)
    })
    ftr <- na.omit(ftr)
    rownames(ftr) <- sprintf('G%04X', 1L:nrow(ftr))
    saveRDS(ftr, cache)
    ftr
}

##string chromosome name to numbers
chr2num <- function(chr)
{
    if(chr == 'X')
        ret <- 23L
    else if(chr == 'Y')
        ret <- 24L
    else if(chr == 'XY')
        ret <- 25L
    else if(chr == 'MT')
        ret <- 26L
    else
        ret <- as.integer(chr)
    ret
}
