
Gencode.Match <- function(list1, list2) {
    # return the indices in list2 that match list1
    list1 <- sub('\\.\\d+', '', list1)
    list2 <- sub('\\.\\d+', '', list2)
    return( match(list1, list2))
}

Gencode.Intersect <- function(list1, list2) {
    # return the indices in list2 that match list1
    list1 <- sub('\\.\\d+', '', list1)
    list2 <- sub('\\.\\d+', '', list2)
    return( intersect(list1, list2))
}

Gencode.strip <- function(list1) {
    return( sub('\\.\\d+', '', list1))
}
