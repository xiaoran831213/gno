## make data frame out of list of lists
## lol ---- the list of lists. each sub list represent a row
## in the data frame
lol2tab <- function(lol, id.from=seq_along)
{
    if(!is.null(id.from))
    {
        lol <- mapply(function(i,l)
        {
            c(idx=i,l)
        }, id.from(lol), lol)
    }
    
    ## create the list to hold columns
    val <- list()

    ## collect names
    hdr <- Reduce(union, lapply(lol, names))

    ## remove empty names from unnamed row elements
    hdr <- setdiff(hdr, "")

    ## create column lists
    for(col in hdr)
        val[[col]] <- as.list(rep(NA, length(lol)))
    
    ## iterate through row lists, fill up column lists
    for(i in 1L:length(lol))
    {
        row <- lol[[i]]
        for(col in names(row))
        {
            val[[col]][[i]] <- row[[col]]
        }
    }

    ## turn column lists into column vectors, then data.frame
    val <- lapply(val, unlist)
    
    val <- data.frame(val, stringsAsFactors = F)
    val
}

## turn table to a list of lists, each sublists represent a row in the
## original table
tab2lol <- function(tab, name = c('rownames', 'combined', 'none'))
{
    lol <- sapply(1L:nrow(tab), function(i) as.list(tab[i,]), simplify = F)

    name <- switch(
        match.arg(name),
        rownames = rownames(tab),
        combined = do.call(paste, c(tab, sep='.')),
        none = NULL)
    names(lol) <- name
    lol
}

tab2lov <- function(tab, name = c('rownames', 'combined', 'none'))
{
    lov <- sapply(1L:nrow(tab), function(i) drop(as.matrix(tab[i,])), simplify = F)
    
    name <- switch(
        match.arg(name),
        row=rownames(tab),
        com=do.call(paste, c(tab, sep='.')),
        non=NULL)

    names(lov) <- name
    lov
}
