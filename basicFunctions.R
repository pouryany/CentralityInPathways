### Function for calculating Katz-Source-Sink
newpath.centrality  <- function(adj.matrix, alpha, beta){


    eye <- diag(nrow(adj.matrix))
    cent.out  <- solve(eye - alpha * adj.matrix)
    cent.in   <- solve(eye - alpha * t(adj.matrix))
    cent.tot  <- (cent.out) + beta * (cent.in)
    return(list("SSC" = cent.tot, "Sink" =cent.in, "Source" =cent.out))

}

### Given a list of eigenvalues, find the largest real component
largest.eigen       <- function(b) {
    c <-  b %>% unlist()  %>% Im() == 0
    c %<>% subset.default(x = b,.) %>% Re() %>% max()
    return(c)
}

### Normalizes a list of values to zero mean and 1 standard deviation
zero.one.normalize    <- function(cent.list){
    if(max(cent.list) == min(cent.list)){
        return(rep(1,length(cent.list)))
    } else{
        #normalizedVec <- (cent.list - min(cent.list))/(max(cent.list) - min(cent.list))
        normalizedVec <- (cent.list - mean(cent.list))/(sd(cent.list))

        return(normalizedVec)
    }
}

### Semi laplace calculator, works with normalized connectivty matrix
semi.laplace <- function(some.Matrix){

    inv.diag <- 1/(rowSums(some.Matrix) + 0.01)
    inv.diag[!is.finite(inv.diag)] <- 0
    norm.laplace.esque <-  diag(1,length(inv.diag)) - (diag(inv.diag) %*% some.Matrix)
    #norm.laplace.esque <-  Ginv(norm.laplace.esque)
    norm.laplace.esque <-  solve(norm.laplace.esque)
    return(norm.laplace.esque)
}
