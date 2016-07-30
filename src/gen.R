library(vcfR)
library(Matrix)

## load one gene from VCF file
gen.vcf <- function(vcf)
{
    ## load VCF
    vc <-read.vcfR(vcf)
    
    ## extract non-synonymous variants
    gt <- extract.gt(vc)

    ## sample names, counts
    sbj <- colnames(gt)

    ## the GT to dosage dictionary.
    dc <- c(
        0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
        1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
        1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
        1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
        1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)
    ## bi-allelic format to dosage format
    gt <- as.integer(gsub('[/|]', '', gt))
    gt <- dc[gt + 1L]
    gt <- Matrix(gt, length(gt)/length(sbj), length(sbj))
    colnames(gt) <- sbj

    ## genomic map
    mp <- vcfR2tidy(vc, T)$fix
    rt <- with(mp,
    {
        list(
            gen = GEN[1],               # smybol
            chr = as.integer(CHROM[1]), # chromosome
            bp0 = head(POS, 1),         # start
            bp1 = tail(POS, 1),         # end
            vcf = vcf)                  # VCF file
    })

    ## cleanup the map
    rt$map <- within(mp,
    {
        CHROM <- as.integer(CHROM)
        FILTER <- NULL
        QUAL <- NULL
        GEN <- NULL
        FNC <- as.factor(FNC)
        NSN <- NSN > 0L
    })

    ## the genotype matrix
    rt$gmx <- gt

    rt
}
