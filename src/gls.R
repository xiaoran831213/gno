library(VariantAnnotation)

## load gene list
gls.load <- function()
{
    cache <- '.gene.rds'
    if(file.exists(cache))
    {
        gene <- readRDS(cache)
    }
    else
    {
        ## load gene transcript and id
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb <- as.list(TxDb.Hsapiens.UCSC.hg38.knownGene)
        tran <- txdb$transcripts
        gene <- txdb$gene
        gene <- within(
            tran,
        {
            tx_chrom <- as.character(tx_chrom)
            gene_id <- as.integer(gene$gene_id[match(tx_id, gene$tx_id)])
        })
        gene <- subset(gene, !is.na(gene_id), select = -c(tx_id, tx_name, tx_strand))
        
        ## replace character chomosome name with integers
        seqs <- paste('chr', c(1:22, c('X', 'Y', 'M')), sep = '')
        gene <- within(gene, tx_chrom <- match(tx_chrom, seqs))
        gene <- subset(gene, !is.na(tx_chrom))
        
        ## collapse duplicates, rename columns
        gene <- with(
            aggregate(cbind(bp1=tx_start, bp2=-tx_end) ~ tx_chrom + gene_id, gene, min),
        {
            data.frame(gid=gene_id, chr=tx_chrom, bp1=bp1, bp2=-bp2)
        })
        
        ## append symbol and full names
        library(org.Hs.eg.db)
        sgid <- as.character(gene$gid)
        gene <- within(
            gene,
        {
            fnm <- as.list(org.Hs.egGENENAME)[sgid]
            sym <- as.list(org.Hs.egSYMBOL)[sgid]
        })
        saveRDS(gene, cache)
    }
    gene
}
