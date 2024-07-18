
#' @title Score or loading plot for the O2PLS results
#' @importFrom ggplot2 ggplot aes theme_classic geom_point stat_ellipse
#' @importFrom ggplot2 geom_text scale_colour_brewer scale_color_manual
#' @importFrom ggplot2 xlab ylab coord_flip theme element_text
#' @importFrom ggplot2 geom_bar stat_ellipse
#' @importFrom ggrepel geom_text_repel 
#' @importFrom stats reorder
#' @param x an O2pls object
#' @param type score or loading 
#' @param var specify Xjoint
#' @param group color used for score plot
#' @param ind which components to be used for score plot or loading plot
#' @param color color used for score or loading plot
#' @param top the number of largest loading value to plot
#' @param ellipse TRUE/FALSE
#' @param order order by the value or not
#' @param pt.size point size
#' @param label plot label or not (TRUE/FALSE) 
#' @param label.size label size
#' @param repel use ggrepel to show the label or not
#' @param rotation flip the figure or not (TRUE/FALSE)
#' @param ... For consistency 
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- matrix(rnorm(50),10,5)
#' fit <- o2pls(X,Y,2,1,1)
#' plot(fit, type="score")
#' @return a ggplot2 object
#' @export
#' @author Kai Guo
plot.O2pls <- function(x, type = "score", var = "Xjoint",group = NULL, 
                       ind = c(1,2), color = NULL,
                       top = 20, ellipse = TRUE, order = FALSE,
                       pt.size = 3, label = TRUE, label.size = 4,
                       repel=TRUE,rotation =FALSE,...){
    if(type == "score"){
        if(var=="Xjoint"){
            dd <- scores(x,score="Xjoint")
            va <- x@results$varXj
        }else if(var == "Yjoint"){
            dd <- scores(x,score="Yjoint")
            va <- x@results$varYj
        }else if(var == "Xorth"){
            dd <- scores(x,score="Xorth")
            va <- x@results$varXorth
        }else if(var == "Yorth"){
            dd <- scores(x,score="Yorth")
            va <- x@results$varYorth
        }else{
            stop('Please specify the score: ["Xjoint", "Yjoint", "Xorth", "Yorth"] \n')
        }
        dd <- as.data.frame(dd)
        if(ncol(dd)==1){
            dd <- dd[, 1, drop = FALSE]
            dd[,2]<-1:nrow(dd)
            dd <- dd[,c(2,1)]
        }else{
            dd <- dd[,ind]
        }
        colnames(dd)[1:2]<-c("L1","L2")
        dd$lab <- rownames(dd)
        if(!is.null(group)){
            dd$Group <- group
            p <- ggplot(dd, aes(L1, L2, color=Group)) + 
                geom_point(size = pt.size)
        }else{
            p <- ggplot(dd, aes(L1, L2)) + geom_point(size = pt.size)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=lab), size = label.size)
            }else{
                p <- p + geom_text(aes(label=lab), size = label.size)
            }
        }
        if(!is.null(color)){
            if(length(color)!=length(unique(group))){
                p <- p + scale_colour_brewer(palette = "Set1")
            }else{
                p <- p + scale_color_manual(values = color)
            }
        }else{
            p <- p + scale_colour_brewer(palette = "Set1")
        }
        if(isTRUE(ellipse)&!is.null(group)){
            p <- p + stat_ellipse()
        }
        if(x@params$nc==1){
            p <- p + xlab(paste0("LV",ind[1],"(",round(va[1]*100,2),"%)"))+
                ylab("")
        }else{
        p <- p + xlab(paste0("LV",ind[1],"(",round(va[1]*100,2),"%)"))+
            ylab(paste0("LV",ind[2],"(",round(va[2]*100,2),"%)"))
        }
        p <- p + theme_classic(base_size =14)
    }
    if(type == "loading"){
        if(var=="Xjoint"){
            dd <- loadings(x,loading = "Xjoint")
        }else if(var == "Yjoint"){
            dd <- loadings(x,loading = "Yjoint")
        }else if(var == "Xorth"){
            dd <- loadings(x,loading = "Xorth")
        }else if(var == "Yorth"){
            dd <- loadings(x,loading = "Yorth")
        }else{
            stop('Please specify the loading: ["Xjoint", "Yjoint", "Xorth", "Yorth"] \n')
        }
        dd <- as.data.frame(dd)
        if(length(ind)>1){
            ind = ind[1]
        }
        dd <- dd[,ind,drop=F]
        colnames(dd) <- "LV"
        dd$lab <- rownames(dd)
        dd$value <- round(dd[,1],2)
        dd <- dd[order(abs(dd$LV),decreasing = T),]
        dd <- dd[1:top, ]
        if(!is.null(color)){
            color <- color[1]
        }else{
            color <- "cyan4"
        }
        if(isTRUE(order)){
            p <- ggplot(dd, aes(reorder(lab,abs(LV)),LV))+geom_bar(stat ="identity",fill=color)
        }else{
            p <- ggplot(dd, aes(lab,LV))+geom_bar(stat ="identity",fill=color)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=value))
            }else{
                p <- p + geom_text(aes(label = value),size = label.size)
            }
        }
        p <- p + theme_classic(base_size =14)
        if(isTRUE(rotation)){
            p <- p + coord_flip()+xlab("Loading") + ylab("")
        }else{
            p <- p + xlab("") + ylab("Loading") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
        }
    }
    
    p
    
}
#' @title Score, VIP or loading plot for the O2PLS results
#' @importFrom ggplot2 ggplot aes theme_classic geom_point stat_ellipse
#' @importFrom ggplot2 geom_text scale_colour_brewer scale_color_manual
#' @importFrom ggplot2 xlab ylab coord_flip theme element_text
#' @importFrom ggplot2 geom_bar stat_ellipse
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats reorder
#' @param x an o2plsda object
#' @param type score, vip or loading 
#' @param group color used for score plot
#' @param ind which components to be used for score plot or loading plot
#' @param color color used for score or loading plot
#' @param top the number of largest loading value to plot
#' @param ellipse TRUE/FALSE
#' @param order order by the value or not
#' @param pt.size point size
#' @param label plot label or not (TRUE/FALSE) 
#' @param label.size label size
#' @param repel use ggrepel to show the label or not
#' @param rotation flip the figure or not (TRUE/FALSE)
#' @param ... For consistency 
#' @return a ggplot2 object
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- matrix(rnorm(50),10,5)
#' fit <- o2pls(X,Y,2,1,1)
#' yy <- rep(c(0,1),5)
#' fit0 <- oplsda(fit,yy,2)
#' plot(fit0, type="score", group = factor(yy))
#' @export
#' @author Kai Guo
plot.o2plsda <- function(x,type = "score",group = NULL, 
                        ind = c(1,2), color = NULL,
                        top = 20, ellipse = TRUE,order=FALSE,
                        pt.size = 3, label = TRUE, label.size = 4,
                        repel=FALSE, rotation =FALSE,...){
    if(type == "score"){
        dd <- scores(x)
        if(x$ncomp==1){
            va <- x$xvar[2,1] 
        }else{
            va <- x$xvar[2,ind]
        }
        
        dd <- as.data.frame(dd)
        if(ncol(dd)==1){
            dd <- dd[, 1, drop = FALSE]
            dd[,2]<-1:nrow(dd)
            dd <- dd[,c(2,1)]
        }else{
            dd <- dd[,ind]
        }
        colnames(dd)[1:2]<-c("L1","L2")
        dd$lab <- rownames(dd)
        if(!is.null(group)){
            dd$Group <- group
            p <- ggplot(dd, aes(L1, L2, color=Group)) + 
                geom_point(size = pt.size)
        }else{
            p <- ggplot(dd, aes(L1, L2)) + geom_point(size = pt.size)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=lab), size = label.size)
            }else{
                p <- p + geom_text(aes(label=lab), size = label.size)
            }
        }
        if(!is.null(color)){
            if(length(color)!=length(unique(group))){
                p <- p + scale_colour_brewer(palette = "Set1")
            }else{
                p <- p + scale_color_manual(values = color)
            }
        }else{
            p <- p + scale_colour_brewer(palette = "Set1")
        }
        if(isTRUE(ellipse)&!is.null(group)){
            p <- p + stat_ellipse()
        }
        if(x$ncomp==1){
            p <- p + xlab(paste0("LV",ind[1],"(",round(va[1]*100,2),"%)"))+
                ylab("")
        }else{
        p <- p + xlab(paste0("LV",ind[1],"(",round(va[1]*100,2),"%)"))+
            ylab(paste0("LV",ind[2],"(",round(va[2]*100,2),"%)"))
        }
        p <- p + theme_classic(base_size =14)
    }
    if(type == "vip"){
        dd <- vip(x)
        dd <- as.data.frame(dd)
        if(length(ind)>1){
            ind = ind[1]
        }
        dd <- dd[,ind,drop=F]
        colnames(dd) <- "LV"
        dd$lab <- rownames(dd)
        dd$value <- round(dd[,1],2)
        dd <- dd[order(abs(dd$LV),decreasing = T),]
        dd <- dd[1:top, ]
        if(!is.null(color)){
            color <- color[1]
        }else{
            color <- "cyan4"
        }
        if(isTRUE(order)){
            p <- ggplot(dd, aes(reorder(lab,LV),LV))+geom_bar(stat ="identity",fill=color)
        }else{
            p <- ggplot(dd, aes(lab,LV))+geom_bar(stat ="identity",fill=color)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=value))
            }else{
                p <- p + geom_text(aes(label = value),size = label.size)
            }
        }
        p <- p + theme_classic(base_size =14)
        if(isTRUE(rotation)){
            p <- p + coord_flip()+xlab("VIP") + ylab("")
        }else{
            p <- p + xlab("") + ylab("VIP") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
        }
    }
    if(type == "loading"){
        dd <- loadings(x,loading = "Xloading")
        dd <- as.data.frame(dd)
        if(length(ind)>1){
            stop("Please specify one loading component\n")
        }
        dd <- dd[,ind,drop=F]
        colnames(dd) <- "LV"
        dd$lab <- rownames(dd)
        dd$value <- round(dd[,1],2)
        dd <- dd[order(abs(dd$LV),decreasing = T),]
        dd <- dd[1:top, ]
        if(!is.null(color)){
            color <- color[1]
        }else{
            color <- "cyan4"
        }
        if(isTRUE(order)){
            p <- ggplot(dd, aes(reorder(lab,abs(LV)),LV))+geom_bar(stat ="identity",fill=color)
        }else{
            p <- ggplot(dd, aes(lab,LV))+geom_bar(stat ="identity",fill=color)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=value))
            }else{
                p <- p + geom_text(aes(label = value),size = label.size)
            }
        }
        p <- p + theme_classic(base_size =14)
        if(isTRUE(rotation)){
            p <- p + coord_flip()+xlab("Loading") + ylab("")
        }else{
            p <- p + xlab("") + ylab("Loading") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
        }
    }
    p
}

#' @title Score, VIP or loading plot for the plsda results
#' @importFrom ggplot2 ggplot aes theme_classic geom_point stat_ellipse
#' @importFrom ggplot2 geom_text scale_colour_brewer scale_color_manual
#' @importFrom ggplot2 xlab ylab coord_flip theme element_text
#' @importFrom ggplot2 geom_bar stat_ellipse
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats reorder
#' @param x an plsda object
#' @param type score, vip or loading 
#' @param group color used for score plot
#' @param ind which components to be used for score plot or loading plot
#' @param color color used for score or loading plot
#' @param top the number of largest loading value to plot
#' @param ellipse TRUE/FALSE
#' @param order order by the value or not
#' @param pt.size point size
#' @param label plot label or not (TRUE/FALSE) 
#' @param label.size label size
#' @param repel use ggrepel to show the label or not
#' @param rotation flip the figure or not (TRUE/FALSE)
#' @param ... For consistency 
#' @examples 
#' X <- matrix(rnorm(500),10,50)
#' Y <- rep(c("a","b"),each=5)
#' fit0 <- plsda(X,Y,2)
#' plot(fit0, type = "score", group = factor(Y))
#' @return a ggplot2 object
#' @export
#' @author Kai Guo
plot.plsda <- function(x,type = "score",group = NULL, 
                         ind = c(1,2), color = NULL,
                         top = 20, ellipse = TRUE,order=FALSE,
                         pt.size = 3, label = TRUE, label.size = 4,
                         repel=FALSE, rotation =FALSE,...){
    if(type == "score"){
        dd <- scores(x)
        va <- x$xvar[ind]
        dd <- as.data.frame(dd)
        if(ncol(dd)==1){
            dd <- dd[, 1, drop = FALSE]
            dd[,2]<-1:nrow(dd)
            dd <- dd[,c(2,1)]
        }else{
            dd <- dd[,ind]
        }
        colnames(dd)[1:2]<-c("L1","L2")
        dd$lab <- rownames(dd)
        if(!is.null(group)){
            dd$Group <- group
            p <- ggplot(dd, aes(L1, L2, color=Group)) + 
                geom_point(size = pt.size)
        }else{
            p <- ggplot(dd, aes(L1, L2)) + geom_point(size = pt.size)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=lab), size = label.size)
            }else{
                p <- p + geom_text(aes(label=lab), size = label.size)
            }
        }
        if(!is.null(color)){
            if(length(color)!=length(unique(group))){
                p <- p + scale_colour_brewer(palette = "Set1")
            }else{
                p <- p + scale_color_manual(values = color)
            }
        }else{
            p <- p + scale_colour_brewer(palette = "Set1")
        }
        if(isTRUE(ellipse)&!is.null(group)){
            p <- p + stat_ellipse()
        }

        p <- p + xlab(paste0("LV",ind[1],"(",round(va[1]*100,2),"%)"))+
            ylab(paste0("LV",ind[2],"(",round(va[2]*100,2),"%)"))
        p <- p + theme_classic(base_size =14)
    }
    if(type == "vip"){
        dd <- x$vip
        dd <- as.data.frame(dd)
        if(length(ind)>1){
            ind = ind[1]
        }
        dd <- dd[,ind,drop=F]
        colnames(dd) <- "LV"
        dd$lab <- rownames(dd)
        dd$value <- round(dd[,1],2)
        dd <- dd[order(abs(dd$LV),decreasing = T),]
        dd <- dd[1:top, ]
        if(!is.null(color)){
            color <- color[1]
        }else{
            color <- "cyan4"
        }
        if(isTRUE(order)){
            p <- ggplot(dd, aes(reorder(lab,LV),LV))+geom_bar(stat ="identity",fill=color)
        }else{
            p <- ggplot(dd, aes(lab,LV))+geom_bar(stat ="identity",fill=color)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=value))
            }else{
                p <- p + geom_text(aes(label = value),size = label.size)
            }
        }
        p <- p + theme_classic(base_size =14)
        if(isTRUE(rotation)){
            p <- p + coord_flip()+xlab("VIP") + ylab("")
        }else{
            p <- p + xlab("") + ylab("VIP") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
        }
    }
    if(type == "loading"){
        dd <- loadings(x)
        dd <- as.data.frame(dd)
        if(length(ind)>1){
            stop("Please specify one loading component\n")
        }
        dd <- dd[,ind,drop=F]
        colnames(dd) <- "LV"
        dd$lab <- rownames(dd)
        dd$value <- round(dd[,1],2)
        dd <- dd[order(abs(dd$LV),decreasing = T),]
        dd <- dd[1:top, ]
        if(!is.null(color)){
            color <- color[1]
        }else{
            color <- "cyan4"
        }
        if(isTRUE(order)){
            p <- ggplot(dd, aes(reorder(lab,abs(LV)),LV))+geom_bar(stat ="identity",fill=color)
        }else{
            p <- ggplot(dd, aes(lab,LV))+geom_bar(stat ="identity",fill=color)
        }
        if(isTRUE(label)){
            if(isTRUE(repel)){
                p <- p + geom_text_repel(aes(label=value))
            }else{
                p <- p + geom_text(aes(label = value),size = label.size)
            }
        }
        p <- p + theme_classic(base_size =14)
        if(isTRUE(rotation)){
            p <- p + coord_flip()+xlab("Loading") + ylab("")
        }else{
            p <- p + xlab("") + ylab("Loading") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1))
        }
    }
    p
}
