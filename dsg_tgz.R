read <- function(tgz)
{
    tmp <- tempfile('tmp', '.')
    dir.create(tmp)

    ## the return list
    ret <- list(name=sub('.tgz$', '', tgz))

    ## decompose, unpack the tarball
    untar(tgz, compressed = 'gzip', exdir = tmp)

    ## the genomic map
    con <- gzfile(file.path(tmp, 'map.txt.gz'))
    ret$map <- read.table(con)
#    close(con)

    ## the genomic matrix
    con <- gzfile(file.path(tmp, 'mtx.txt.gz'))
    ret$mtx <-  scan(con, integer(), na.strings = '3')
#    close(con)

    ret
}
