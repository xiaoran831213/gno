ped.load <- function(file)
{
    ped <- read.table(
        file,
        col.names = c('fid', 'iid', 'pid', 'mid', 'sex', 'phe'),
        colClasses = c(rep('integer', 5L), 'numeric'))
    rownames(ped) <- ped$iid
    ped
}

ped.save <- function(ped, file)
{
    write.table(ped, file, quote=F, sep='\t', row.names = F, col.names = F)
}

ped.parent <- function(ped)
{
    child.father <- subset(ped, pid > 0, select = c(iid, pid))
    child.mather <- subset(ped, mid > 0, select = c(iid, mid))
    names(child.father) <- c('cid', 'pid')
    names(child.mather) <- c('cid', 'pid')
    rbind(child.father, child.mather)
}

## clear non-exist parent (pid or mid not in iid)
ped.clr <- function(ped)
{
    within(
        ped,
    {
        pid[!pid %in% iid] <- 0
        mid[!mid %in% iid] <- 0
    })
}
ped.nuclear <- function(ped, both = FALSE, missing = NA)
{
    if(both)
        filter <- expression(pid > 0 & mid > 0)
    else
        filter <- expression(pid > 0 | mid > 0)
    subset(ped, eval(filter), c(iid, pid, mid))
}
