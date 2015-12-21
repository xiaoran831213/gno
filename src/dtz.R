## read dosage data in tar.gz
readDTZ <- function(dtz)
{
    ## creat temporary directory
    tmp <- tempfile('tmp', '.')
    dir.create(tmp)

    ## decompose, unpack the tarball
    untar(dtz, compressed = 'gzip', exdir = tmp)

    ## the VCF header, with meta data
    gz <- gzfile(file.path(tmp, 'hdr.txt.gz'))
    mt <- grep('##ssn', readLines(gz), value=T)
    close(gz)
    mt <- unlist(strsplit(sub('^##', '', mt), ','))
    
    ## the genomic map
    gz <- gzfile(file.path(tmp, 'map.txt.gz'))
    map <- read.table(text=readLines(gz))
    close(gz)
    
    ## the subject IDs
    gz <- gzfile(file.path(tmp, 'sid.txt.gz'))
    sid <-  scan(gz, "", quiet = T)
    close(gz)
    
    ## the variant IDs
    gz <- gzfile(file.path(tmp, 'vid.txt.gz'))
    vid <-  scan(gz, "", quiet = T)
    close(gz)
    
    ## the genomic matrix
    gz <- gzfile(file.path(tmp, 'mtx.txt.gz'))
    mtx <-  matrix(
        scan(gz, integer(), na.strings = '3', quiet = T),
        nrow = length(vid), ncol = length(sid),
        byrow = TRUE,
        dimnames = list(vid=vid, sid=sid))
    close(gz)

    ## delete the temporary directory
    unlink(tmp, recursive = T, force = T)

    ## the return list
    ret <- list(
        src=dtz,                        # source name
        ifo=mt,                         # extra information
        map=map,                        # genomic map
        gmx=mtx,                        # genotype matrix
        sid=sid,                        # subject id
        vid=vid                         # variant id
    )
    invisible(ret)                      # prevent screen flood
}
