## load gene list
ftr.load <- function(recache = FALSE)
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
    tmp <- tempfile(fileext = 'gz')
    download.file(src, tmp, method = 'wget', mode = 'wb')
    ftr <- read.table(gzfile(tmp), header = T, sep = '\t', comment.char = "", as.is = T)
    names(ftr)[1] <- 'feature'
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@9@"]]))##:ess-bp-end:##
    
    ## update the cache then return
    saveRDS(ftr, cache)
    file.remove(tmp)

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


gls.load <- function(recache = FALSE)
{
    cache <- '.genes.rds'

    ## try cached gene list first
    if(file.exists(cache) && !recache)
    {
        gene <- readRDS(cache)
    }
    else
    {
        gene <- ftr.load(F)
        
        ## extract gene info from feature table
        gene <- subset(gene, feature=='gene')
        gene <- subset(gene, assembly_unit == 'Primary Assembly')
#            select = c(chromosome, start, end, GeneID, symbol, name))
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@6@"]]))##:ess-bp-end:##
        
        ## reformat, use short column names
        gene <- with(
            gene,
        {
            data.frame(
                chr = chromosome, #sapply(chromosome, chr2num, USE.NAMES = F),
                bp1 = as.integer(start),
                bp2 = as.integer(end),
                gid = GeneID,
                smb = symbol,
                fnm = name,
                ## ifo = sprintf(
                ##     'clz=%s,sq_typ=%s,attr=%s',
                ##     class, seq_type, attributes),
                stringsAsFactors = F)
        })
        rownames(gene) <- sprintf('G%04X', 1L:nrow(gene))
        
        ## update cache then return
        saveRDS(gene, cache)
    }
    gene
}
