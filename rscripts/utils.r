tic = function() {
    tic = proc.time()['elapsed']
    assign('.tic', tic, envir= baseenv())
    invisible(tic)
}

toc = function() {
    toc = proc.time()['elapsed']
    print(toc - get('.tic', envir= baseenv()))
    invisible(toc)
}

mem.list = function() {
    sort( sapply(ls(envir=baseenv()),function(x){object.size(get(x))}))
}

mem.total = function() {
    total = sum(sapply(ls(envir=baseenv()),function(x){object.size(get(x))}))
    cat('total mem = ', total, '\n')
    return(total)
}
