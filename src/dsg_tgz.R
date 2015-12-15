read_tgz <- function(tgz)
{
    tmp <- tempfile('tmp', '.')
    dir.create(tmp)

    ## the return list
    ret <- list(src=sub('.tgz$', '', tgz))

    ## decompose, unpack the tarball
    untar(tgz, compressed = 'gzip', exdir = tmp)

    ## the genomic map
    gz <- gzfile(file.path(tmp, 'map.txt.gz'))
    ret$map <- read.table(gz)
    close(gz)
    
    ## the subject IDs
    gz <- gzfile(file.path(tmp, 'sid.txt.gz'))
    ret$sid <-  scan(gz, "")
    close(gz)
    
    ## the variant IDs
    gz <- gzfile(file.path(tmp, 'vid.txt.gz'))
    ret$vid <-  scan(gz, "")
    close(gz)
    
    ## the genomic matrix
    gz <- gzfile(file.path(tmp, 'mtx.txt.gz'))
    ret <- within(
        ret,
    {
        mtx <-  matrix(
            scan(gz, integer(), na.strings = '3'),
            nrow = length(vid), ncol = length(sid),
            byrow = TRUE,
            dimnames = list(vid=vid, sid=sid))
    })
    close(gz)

    file.remove(tmp)
    ret
}
