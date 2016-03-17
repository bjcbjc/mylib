
library(ggplot2)
library(reshape2)
library(grid)

numSubplot = function(nplots) {
    b = sqrt(nplots)
    l = floor(b)
    u = ceiling(b)
    r = round(b)

    if (r == u) {
        nr = u
        nc = u
    } else {
        nr = l
        nc = u
    }
    return(c(nr, nc))
}

matrixHeatmap = function(mtx,
    color= c(scales::muted('red'), 'white', scales::muted('blue')),
    clim= NULL, gridColor= NULL, xTickLabel= NULL, yTickLabel= NULL, textData= NULL, textSize=8) {

    if(is.null(clim)) {
        clim = c(min(mtx), max(mtx))
        midPoint = mean(clim)
    } else {
        midPoint = mean(clim)
        mtx[mtx < clim[1]] = clim[1]
        mtx[mtx > clim[2]] = clim[2]
    }

    if(!is.null(textData)) {
        if (!all(dim(textData) == dim(mtx))) {
            stop('textData does not have the same dim as mtx')
        }
    }
    
    # mtx is a matrix
    data = melt(t(mtx))
    if(!is.null(textData)) {
        data$mtxText = as.vector(t(textData))
        data$mtxText[is.na(data$mtxText)] = ''
    }
    colnames(data)[1:2] = c('Var1', 'Var2')
    
    p = ggplot(data, aes(x= as.factor(Var1), y= as.factor(Var2)))
    if (is.null(gridColor)) {
        p = p + geom_tile(aes(fill= value)) 
    } else {
        p = p + geom_tile(aes(fill= value), color= gridColor)
    }
    p = p + 
        scale_fill_gradient2(high= color[1], mid= color[2], low= color[3],
                             midpoint= midPoint) +
        scale_x_discrete(expand= c(0,0)) +
        scale_y_discrete(expand= c(0,0)) +
        theme(axis.ticks= element_blank(),
              axis.text.x= element_text(angle = 90, hjust= 0),
              axis.title.x= element_blank(),
              axis.title.y= element_blank(),
              legend.title= element_blank(),
              panel.background= element_rect(color= 'white', fill= 'white'))
    if (!is.null(xTickLabel)) {
        p = p %+% scale_x_discrete(labels= xTickLabel)
    }
    if (!is.null(yTickLabel)) {
        p = p %+% scale_y_discrete(labels= yTickLabel)
    }
    if (!is.null(textData)) {
        p = p %+% geom_text(aes(fill= value, label= mtxText), size= textSize)
    }
        
    return(p)
}

hclustOrder = function(data, distance='cor') {
    #cluster columns
    if (distance == 'cor') {
        idx = hclust(as.dist(1-cor(data)))$order
    } else {
        idx = hclust(dist(t(data), distance))$order
    }
    return(idx)
}
