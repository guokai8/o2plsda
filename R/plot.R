
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


#' Plot method for Sparse PLS-DA results
#' 
#' @param x A sparse_plsda object (corrected class name)
#' @param type Character. Type of plot: "score", "loading", "vip", "selection"
#' @param component Integer or vector. Which component(s) to plot
#' @param group Factor. Grouping variable (usually the Y classes)
#' @param top Integer. Number of top variables to display
#' @param color Custom colors for groups
#' @param ellipse Logical. Add confidence ellipses to score plots
#' @param ... Additional plotting parameters
#' 
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' p <- 150
#' 
#' X <- matrix(rnorm(n * p), n, p)
#' # Add class-specific signals
#' classes <- factor(rep(c("A", "B", "C"), length.out = n))
#' X[classes == "A", 1:20] <- X[classes == "A", 1:20] + 1.5
#' X[classes == "B", 21:40] <- X[classes == "B", 21:40] + 1.5
#' X[classes == "C", 41:60] <- X[classes == "C", 41:60] + 1.5
#' 
#' colnames(X) <- paste0("Gene_", 1:p)
#' 
#' # Fit sparse PLS-DA (using corrected class name)
#' sparse_plsda_fit <- sparse_plsda(X, classes, nc = 2, keepX = c(30, 25))
#' class(sparse_plsda_fit) <- "sparse_plsda"  # Ensure correct class
#' 
#' # Score plot
#' plot(sparse_plsda_fit, type = "score", group = classes, ellipse = TRUE)
#' 
#' # Loading plot
#' plot(sparse_plsda_fit, type = "loading", component = 1, top = 20)
#' 
#' # VIP plot
#' plot(sparse_plsda_fit, type = "vip", component = 1, top = 25)
#' 
#' @export
plot.sparse_plsda <- function(x, type = "score", component = c(1, 2),
                              group = NULL, top = 20, color = NULL, 
                              ellipse = TRUE, ...) {
    
    # Import required libraries
    library(ggplot2)
    library(dplyr)
    
    type <- match.arg(type, c("score", "loading", "vip", "selection", "classification"))
    
    if(type == "score") {
        return(plot_sparse_plsda_scores(x, component, group, color, ellipse, ...))
    } else if(type == "loading") {
        return(plot_sparse_plsda_loadings(x, component, top, ...))
    } else if(type == "vip") {
        return(plot_sparse_plsda_vip(x, component, top, ...))
    } else if(type == "selection") {
        return(plot_sparse_plsda_selection(x, ...))
    } else if(type == "classification") {
        return(plot_sparse_plsda_classification(x, component, group, ...))
    }
}

#' Plot Sparse PLS-DA scores
#' @keywords internal
plot_sparse_plsda_scores <- function(x, component, group, color, ellipse, ...) {
    # Extract scores (adapt to actual object structure)
    if(!is.null(x@scores)) {
        scores <- x@scores
    } else if(!is.null(x@score)) {
        scores <- x@score
    } else {
        stop("Cannot find scores in sparse_plsda object")
    }
    
    if(length(component) == 1) {
        # 1D score plot
        plot_data <- data.frame(
            Score = scores[, component],
            Index = 1:nrow(scores),
            Sample = rownames(scores) %||% paste0("Sample_", 1:nrow(scores))
        )
        
        if(!is.null(group)) {
            plot_data$Group <- group
            p <- ggplot(plot_data, aes(x = Index, y = Score, color = Group)) +
                geom_point(size = 3) +
                geom_line(alpha = 0.5, aes(group = Group))
        } else {
            p <- ggplot(plot_data, aes(x = Index, y = Score)) +
                geom_point(size = 3) +
                geom_line(alpha = 0.5)
        }
        
        p <- p + labs(x = "Sample Index", 
                      y = paste0("Component ", component),
                      title = "Sparse PLS-DA Scores")
        
    } else {
        # 2D score plot
        plot_data <- data.frame(
            PC1 = scores[, component[1]],
            PC2 = scores[, component[2]],
            Sample = rownames(scores) %||% paste0("Sample_", 1:nrow(scores))
        )
        
        if(!is.null(group)) {
            plot_data$Group <- group
            p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
                geom_point(size = 3, alpha = 0.8)
            
            if(ellipse) {
                p <- p + stat_ellipse(level = 0.95, linewidth = 1)
            }
        } else {
            p <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
                geom_point(size = 3, alpha = 0.8)
        }
        
        p <- p + labs(x = paste0("Component ", component[1]),
                      y = paste0("Component ", component[2]),
                      title = "Sparse PLS-DA Scores")
    }
    
    # Apply colors
    if(!is.null(color) && !is.null(group)) {
        p <- p + scale_color_manual(values = color)
    } else if(!is.null(group)) {
        p <- p + scale_color_brewer(type = "qual", palette = "Set1")
    }
    
    p <- p + theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
}

#' Plot Sparse PLS-DA loadings
#' @keywords internal
plot_sparse_plsda_loadings <- function(x, component, top, ...) {
    # Extract loadings
    if(!is.null(x@loadings)) {
        loadings <- x@loadings
    } else if(!is.null(x@Xloading)) {
        loadings <- x@Xloading
    } else {
        stop("Cannot find loadings in sparse_plsda object")
    }
    
    if(length(component) > 1) {
        component <- component[1]
        warning("Multiple components specified, using first component only")
    }
    
    loading_vals <- loadings[, component]
    var_names <- rownames(loadings) %||% paste0("Var_", 1:nrow(loadings))
    
    # Get top variables
    top_indices <- order(abs(loading_vals), decreasing = TRUE)[1:min(top, length(loading_vals))]
    
    plot_data <- data.frame(
        Variable = var_names[top_indices],
        Loading = loading_vals[top_indices],
        Absolute = abs(loading_vals[top_indices])
    ) %>%
        arrange(desc(Absolute))
    
    plot_data$Variable <- factor(plot_data$Variable, levels = plot_data$Variable)
    
    p <- ggplot(plot_data, aes(x = Variable, y = Loading)) +
        geom_col(aes(fill = Loading > 0), alpha = 0.8) +
        scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
                          name = "Direction", labels = c("Negative", "Positive")) +
        coord_flip() +
        labs(x = "Variables",
             y = paste0("Loading Values - Component ", component),
             title = "Sparse PLS-DA Loadings") +
        theme_classic(base_size = 12) +
        theme(axis.text.y = element_text(size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
}

#' Plot Sparse PLS-DA VIP scores
#' @keywords internal
plot_sparse_plsda_vip <- function(x, component, top, ...) {
    # Extract VIP scores
    if(!is.null(x@vip)) {
        vip_scores <- x@vip
    } else {
        stop("Cannot find VIP scores in sparse_plsda object")
    }
    
    if(length(component) > 1) {
        component <- component[1]
        warning("Multiple components specified, using first component only")
    }
    
    vip_vals <- vip_scores[, component]
    var_names <- rownames(vip_scores) %||% paste0("Var_", 1:nrow(vip_scores))
    
    # Get top variables
    top_indices <- order(vip_vals, decreasing = TRUE)[1:min(top, length(vip_vals))]
    
    plot_data <- data.frame(
        Variable = var_names[top_indices],
        VIP = vip_vals[top_indices]
    ) %>%
        arrange(desc(VIP))
    
    plot_data$Variable <- factor(plot_data$Variable, levels = plot_data$Variable)
    plot_data$Important <- plot_data$VIP > 1  # VIP > 1 considered important
    
    p <- ggplot(plot_data, aes(x = Variable, y = VIP)) +
        geom_col(aes(fill = Important), alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D1E5F0"),
                          name = "Important", labels = c("VIP ≤ 1", "VIP > 1")) +
        coord_flip() +
        labs(x = "Variables",
             y = paste0("VIP Scores - Component ", component),
             title = "Variable Importance in Projection (VIP)") +
        theme_classic(base_size = 12) +
        theme(axis.text.y = element_text(size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
}

#' Plot Sparse PLS-DA variable selection
#' @keywords internal
plot_sparse_plsda_selection <- function(x, ...) {
    # Extract selection information
    if(!is.null(x@selected_vars)) {
        selected_vars <- x@selected_vars
        total_vars <- nrow(x@loadings)
    } else {
        # Calculate from loadings
        loadings <- x@loadings %||% x@Xloading
        selected_vars <- which(rowSums(abs(loadings)) > 1e-10)
        total_vars <- nrow(loadings)
    }
    
    n_selected <- length(selected_vars)
    n_not_selected <- total_vars - n_selected
    percent_selected <- round((n_selected / total_vars) * 100, 1)
    
    # Create summary plot
    plot_data <- data.frame(
        Category = c("Selected", "Not Selected"),
        Count = c(n_selected, n_not_selected),
        Percentage = c(percent_selected, 100 - percent_selected)
    )
    
    p <- ggplot(plot_data, aes(x = "", y = Count, fill = Category)) +
        geom_col(alpha = 0.8) +
        geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
                  position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
        scale_fill_manual(values = c("Selected" = "#2166AC", "Not Selected" = "#D1E5F0")) +
        coord_polar(theta = "y") +
        labs(title = paste0("Variable Selection Summary\n", n_selected, " out of ", total_vars, " variables selected")) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              legend.position = "bottom")
    
    return(p)
}

#' Plot Sparse PLS-DA classification boundaries
#' @keywords internal
plot_sparse_plsda_classification <- function(x, component, group, ...) {
    # This would require the classification results from the model
    # Extract scores and centroids
    scores <- x@scores %||% x@score
    
    if(is.null(group)) {
        stop("Group variable required for classification plot")
    }
    
    # Calculate class centroids
    plot_data <- data.frame(
        PC1 = scores[, component[1]],
        PC2 = scores[, component[2]],
        Group = group
    )
    
    centroids <- plot_data %>%
        group_by(Group) %>%
        summarise(
            PC1_center = mean(PC1),
            PC2_center = mean(PC2),
            .groups = 'drop'
        )
    
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_point(data = centroids, aes(x = PC1_center, y = PC2_center, color = Group),
                   size = 6, shape = 17) +
        stat_ellipse(level = 0.95, linewidth = 1) +
        scale_color_brewer(type = "qual", palette = "Set1") +
        labs(x = paste0("Component ", component[1]),
             y = paste0("Component ", component[2]),
             title = "Sparse PLS-DA Classification") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
}

#' Enhanced plot method for TuneResult objects
#' 
#' @param x A TuneResult object
#' @param type Character. Type of plot: "heatmap", "line", "comparison", "optimal"
#' @param metric Character. Performance metric to visualize
#' @param ... Additional plotting parameters
#' 
#' @examples
#' # Example with tuning results
#' \dontrun{
#' # Assuming you have a TuneResult object from tune_o2pls()
#' tune_results <- tune_o2pls(X, Y, nc_range = 1:3, nx_range = 0:2, ny_range = 0:2)
#' 
#' # Heatmap of results
#' plot(tune_results, type = "heatmap")
#' 
#' # Line plot showing trends
#' plot(tune_results, type = "line")
#' 
#' # Comparison across parameters
#' plot(tune_results, type = "comparison")
#' }
#' 
#' @export
plot.TuneResult <- function(x, type = "heatmap", metric = "mean_score", ...) {
    
    library(ggplot2)
    library(dplyr)
    
    type <- match.arg(type, c("heatmap", "line", "comparison", "optimal", "surface"))
    
    if(type == "heatmap") {
        return(plot_tuning_heatmap_impl(x, metric, ...))
    } else if(type == "line") {
        return(plot_tuning_line_impl(x, metric, ...))
    } else if(type == "comparison") {
        return(plot_tuning_comparison_impl(x, ...))
    } else if(type == "optimal") {
        return(plot_optimal_params_impl(x, ...))
    } else if(type == "surface") {
        return(plot_tuning_surface_impl(x, metric, ...))
    }
}

#' Implementation of tuning heatmap
#' @keywords internal
plot_tuning_heatmap_impl <- function(x, metric, ...) {
    # Extract results from TuneResult object
    results_list <- x@all_results
    
    # Convert to data frame
    results_df <- do.call(rbind, lapply(results_list, function(res) {
        data.frame(
            nc = res$nc,
            nx = res$nx,
            ny = res$ny,
            score = res[[metric]],
            stringsAsFactors = FALSE
        )
    }))
    
    # Create separate heatmaps for each ny value
    ny_values <- unique(results_df$ny)
    
    if(length(ny_values) == 1) {
        # Single heatmap
        p <- ggplot(results_df, aes(x = factor(nx), y = factor(nc), fill = score)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = round(score, 3)), color = "white", size = 3, fontface = "bold") +
            scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen",
                                 midpoint = median(results_df$score, na.rm = TRUE),
                                 name = x@measure) +
            labs(x = "X-Orthogonal Components (nx)",
                 y = "Joint Components (nc)",
                 title = paste0("Parameter Tuning Results (ny = ", ny_values[1], ")")) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
        # Mark optimal point
        optimal <- x@optimal
        p <- p + geom_point(data = data.frame(nx = factor(optimal$nx), nc = factor(optimal$nc)),
                            aes(x = nx, y = nc), color = "blue", size = 4, shape = 1, 
                            stroke = 2, inherit.aes = FALSE)
        
    } else {
        # Multiple heatmaps
        p <- ggplot(results_df, aes(x = factor(nx), y = factor(nc), fill = score)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = round(score, 3)), color = "white", size = 2.5) +
            facet_wrap(~paste("ny =", ny), scales = "free") +
            scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen",
                                 midpoint = median(results_df$score, na.rm = TRUE),
                                 name = x@measure) +
            labs(x = "X-Orthogonal Components (nx)",
                 y = "Joint Components (nc)",
                 title = "Parameter Tuning Results") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }
    
    return(p)
}

#' Implementation of tuning line plot
#' @keywords internal
plot_tuning_line_impl <- function(x, metric, ...) {
    results_list <- x@all_results
    
    results_df <- do.call(rbind, lapply(results_list, function(res) {
        data.frame(
            nc = res$nc,
            nx = res$nx,
            ny = res$ny,
            score = res[[metric]]
        )
    }))
    
    # Line plot showing trends across nc for different nx/ny combinations
    results_df$config <- paste0("nx=", results_df$nx, ", ny=", results_df$ny)
    
    p <- ggplot(results_df, aes(x = nc, y = score, color = config)) +
        geom_line(linewidth = 1, alpha = 0.8) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_brewer(type = "qual", palette = "Set3", name = "Configuration") +
        labs(x = "Joint Components (nc)",
             y = x@measure,
             title = "Parameter Tuning Trends") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right")
    
    # Mark optimal point
    optimal <- x@optimal
    optimal_score <- optimal[[metric]]
    p <- p + geom_point(data = data.frame(nc = optimal$nc, score = optimal_score),
                        aes(x = nc, y = score), color = "red", size = 5, shape = 8,
                        inherit.aes = FALSE)
    
    return(p)
}

#' Implementation of parameter comparison
#' @keywords internal
plot_tuning_comparison_impl <- function(x, ...) {
    results_list <- x@all_results
    
    results_df <- do.call(rbind, lapply(results_list, function(res) {
        data.frame(
            nc = res$nc,
            nx = res$nx,
            ny = res$ny,
            mean_score = res$mean_score,
            sd_score = res$sd_score %||% 0
        )
    }))
    
    # Create comparison across different parameter types
    param_summary <- data.frame(
        Parameter = c("nc", "nx", "ny"),
        Optimal_Value = c(x@optimal$nc, x@optimal$nx, x@optimal$ny),
        Score_at_Optimal = rep(x@optimal$mean_score, 3)
    )
    
    p <- ggplot(param_summary, aes(x = Parameter, y = Optimal_Value)) +
        geom_col(aes(fill = Parameter), alpha = 0.8) +
        geom_text(aes(label = Optimal_Value), vjust = -0.5, fontface = "bold") +
        scale_fill_brewer(type = "qual", palette = "Set2") +
        labs(x = "Parameter Type",
             y = "Optimal Value",
             title = "Optimal Parameter Configuration",
             subtitle = paste0("Best ", x@measure, ": ", round(x@optimal$mean_score, 4))) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              legend.position = "none")
    
    return(p)
}

#' Enhanced plot method for stability_selection objects
#' 
#' @param x A stability_selection object
#' @param component Integer. Which component to visualize
#' @param type Character. Type of plot: "heatmap", "barplot", "threshold", "summary"
#' @param top_n Integer. Number of top variables to show
#' @param ... Additional plotting parameters
#' 
#' @examples
#' \dontrun{
#' # Assuming you have a stability_selection object
#' stability_res <- stability_selection(X, Y, nc = 2, n_bootstrap = 50)
#' 
#' # Stability heatmap
#' plot(stability_res, type = "heatmap", component = 1)
#' 
#' # Stability barplot
#' plot(stability_res, type = "barplot", component = 1, top_n = 30)
#' 
#' # Threshold analysis
#' plot(stability_res, type = "threshold")
#' }
#' 
#' @export
plot.stability_selection <- function(x, component = 1, type = "heatmap", 
                                     top_n = 50, ...) {
    
    library(ggplot2)
    library(dplyr)
    
    type <- match.arg(type, c("heatmap", "barplot", "threshold", "summary", "network"))
    
    if(type == "heatmap") {
        return(plot_stability_heatmap_impl(x, component, top_n, ...))
    } else if(type == "barplot") {
        return(plot_stability_barplot_impl(x, component, top_n, ...))
    } else if(type == "threshold") {
        return(plot_stability_threshold_impl(x, ...))
    } else if(type == "summary") {
        return(plot_stability_summary_impl(x, ...))
    } else if(type == "network") {
        return(plot_stability_network_impl(x, component, ...))
    }
}

#' Implementation of stability heatmap
#' @keywords internal
plot_stability_heatmap_impl <- function(x, component, top_n, ...) {
    # Extract selection probabilities
    selection_prob <- x@selection_probabilities
    
    if(length(dim(selection_prob)) == 3) {
        # 3D array: variables x keepX x components
        prob_matrix <- selection_prob[, , component]
    } else {
        # 2D matrix: variables x keepX
        prob_matrix <- selection_prob
    }
    
    # Select top variables by maximum probability
    max_probs <- apply(prob_matrix, 1, max, na.rm = TRUE)
    top_vars <- order(max_probs, decreasing = TRUE)[1:min(top_n, length(max_probs))]
    
    plot_data <- prob_matrix[top_vars, , drop = FALSE]
    
    # Convert to long format for ggplot
    melted_data <- expand.grid(
        Variable = 1:nrow(plot_data),
        KeepX = 1:ncol(plot_data)
    )
    melted_data$Probability <- as.vector(plot_data)
    melted_data$Variable_Name <- rep(rownames(plot_data), ncol(plot_data))
    melted_data$KeepX_Value <- rep(x$keepX_range, each = nrow(plot_data))
    melted_data$Stable <- melted_data$Probability >= x$threshold
    
    p <- ggplot(melted_data, aes(x = factor(KeepX_Value), y = reorder(Variable_Name, Probability))) +
        geom_tile(aes(fill = Probability), color = "white", linewidth = 0.1) +
        geom_point(data = subset(melted_data, Stable), 
                   aes(x = factor(KeepX_Value), y = Variable_Name), 
                   color = "red", size = 1, shape = 8) +
        scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                             midpoint = 0.5, name = "Selection\nProbability",
                             limits = c(0, 1)) +
        labs(title = paste("Stability Selection - Component", component),
             x = "Number of Variables to Keep (keepX)",
             y = "Variables (ordered by max probability)",
             caption = paste0("Red stars indicate stable variables (≥", x$threshold, " threshold)")) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 6),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.caption = element_text(size = 8))
    
    return(p)
}

#' Implementation of stability barplot
#' @keywords internal
plot_stability_barplot_impl <- function(x, component, top_n, ...) {
    # Extract stable variables for the component
    if(!is.null(x@stable_variables[[paste0("Component_", component)]])) {
        stable_vars <- x@stable_variables[[paste0("Component_", component)]]
        
        if(nrow(stable_vars) == 0) {
            stop("No stable variables found for component ", component)
        }
        
        # Take top variables by probability
        top_stable <- head(stable_vars[order(stable_vars$Max_Probability, decreasing = TRUE), ], top_n)
        
        p <- ggplot(top_stable, aes(x = reorder(Variable_Name, Max_Probability), y = Max_Probability)) +
            geom_col(aes(fill = Max_Probability >= x$threshold), alpha = 0.8) +
            geom_hline(yintercept = x$threshold, linetype = "dashed", color = "red", linewidth = 1) +
            scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D1E5F0"),
                              name = "Stable", labels = c("Below threshold", "Above threshold")) +
            coord_flip() +
            labs(x = "Variables",
                 y = "Selection Probability",
                 title = paste0("Stable Variables - Component ", component),
                 subtitle = paste0("Threshold: ", x$threshold, " (", nrow(top_stable), " variables shown)")) +
            theme_classic(base_size = 12) +
            theme(axis.text.y = element_text(size = 8),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5))
        
    } else {
        # No stable variables found
        p <- ggplot() +
            annotate("text", x = 0.5, y = 0.5, 
                     label = paste0("No stable variables found for component ", component, 
                                    "\nwith threshold ", x$threshold),
                     hjust = 0.5, vjust = 0.5, size = 6) +
            xlim(0, 1) + ylim(0, 1) +
            theme_void() +
            ggtitle(paste0("Stability Results - Component ", component))
    }
    
    return(p)
}

#' Implementation of stability threshold analysis
#' @keywords internal
plot_stability_threshold_impl <- function(x, ...) {
    # Analyze how many variables are stable at different thresholds
    thresholds <- seq(0.1, 0.9, 0.1)
    threshold_results <- data.frame(
        Threshold = thresholds,
        N_Stable = numeric(length(thresholds))
    )
    
    # Count stable variables at each threshold
    selection_prob <- x@selection_probabilities
    if(length(dim(selection_prob)) == 3) {
        # Average across components
        avg_prob <- apply(selection_prob, c(1, 2), mean, na.rm = TRUE)
        max_prob_per_var <- apply(avg_prob, 1, max, na.rm = TRUE)
    } else {
        max_prob_per_var <- apply(selection_prob, 1, max, na.rm = TRUE)
    }
    
    for(i in 1:length(thresholds)) {
        threshold_results$N_Stable[i] <- sum(max_prob_per_var >= thresholds[i], na.rm = TRUE)
    }
    
    # Add current threshold
    current_threshold <- x@threshold
    current_n_stable <- sum(max_prob_per_var >= current_threshold, na.rm = TRUE)
    
    p <- ggplot(threshold_results, aes(x = Threshold, y = N_Stable)) +
        geom_line(linewidth = 2, color = "#2166AC") +
        geom_point(size = 3, color = "#2166AC") +
        geom_vline(xintercept = current_threshold, linetype = "dashed", color = "red", linewidth = 1) +
        geom_point(data = data.frame(Threshold = current_threshold, N_Stable = current_n_stable),
                   aes(x = Threshold, y = N_Stable), color = "red", size = 5, shape = 18) +
        annotate("text", x = current_threshold, y = current_n_stable + max(threshold_results$N_Stable) * 0.1,
                 label = paste0("Current: ", current_threshold, "\n(", current_n_stable, " vars)"),
                 hjust = 0.5, color = "red", fontface = "bold") +
        labs(x = "Stability Threshold",
             y = "Number of Stable Variables",
             title = "Stability Threshold Analysis",
             subtitle = "How threshold choice affects number of selected variables") +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5))
    
    return(p)
}

#' Implementation of stability summary
#' @keywords internal
plot_stability_summary_impl <- function(x, ...) {
    # Create summary across all components
    n_components <- length(x@stable_variables)
    
    if(n_components == 0) {
        stop("No stable variables found in any component")
    }
    
    summary_data <- data.frame(
        Component = character(),
        N_Stable = numeric(),
        Mean_Probability = numeric(),
        Max_Probability = numeric()
    )
    
    for(i in 1:n_components) {
        comp_name <- names(x@stable_variables)[i]
        comp_data <- x@stable_variables[[i]]
        
        if(!is.null(comp_data) && nrow(comp_data) > 0) {
            summary_data <- rbind(summary_data, data.frame(
                Component = comp_name,
                N_Stable = nrow(comp_data),
                Mean_Probability = mean(comp_data$Max_Probability, na.rm = TRUE),
                Max_Probability = max(comp_data$Max_Probability, na.rm = TRUE)
            ))
        }
    }
    
    if(nrow(summary_data) == 0) {
        stop("No stable variables found")
    }
    
    # Create multi-panel summary plot
    p1 <- ggplot(summary_data, aes(x = Component, y = N_Stable)) +
        geom_col(fill = "#2166AC", alpha = 0.8) +
        geom_text(aes(label = N_Stable), vjust = -0.5, fontface = "bold") +
        labs(x = "Component", y = "Number of Stable Variables",
             title = "Stable Variables per Component") +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- ggplot(summary_data, aes(x = Component, y = Mean_Probability)) +
        geom_col(fill = "#D6604D", alpha = 0.8) +
        geom_text(aes(label = round(Mean_Probability, 2)), vjust = -0.5, fontface = "bold") +
        geom_hline(yintercept = x$threshold, linetype = "dashed", color = "red") +
        labs(x = "Component", y = "Mean Selection Probability",
             title = "Average Stability per Component") +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine plots
    library(gridExtra)
    combined_plot <- grid.arrange(p1, p2, ncol = 2,
                                  top = textGrob("Stability Selection Summary", 
                                                 gp = gpar(fontsize = 14, fontface = "bold")))
    
    return(combined_plot)
}

#' Implementation of stability network plot
#' @keywords internal
plot_stability_network_impl <- function(x, component, ...) {
    # Create network based on co-selection patterns
    # This requires the igraph package
    if(!requireNamespace("igraph", quietly = TRUE)) {
        stop("igraph package required for network plots")
    }
    
    library(igraph)
    
    # Extract selection probabilities for the component
    selection_prob <- x@selection_probabilities
    if(length(dim(selection_prob)) == 3) {
        prob_matrix <- selection_prob[, , component]
    } else {
        prob_matrix <- selection_prob
    }
    
    # Get stable variables
    max_probs <- apply(prob_matrix, 1, max, na.rm = TRUE)
    stable_vars <- which(max_probs >= x$threshold)
    
    if(length(stable_vars) < 2) {
        stop("Need at least 2 stable variables for network plot")
    }
    
    # Calculate co-selection matrix
    stable_prob_matrix <- prob_matrix[stable_vars, , drop = FALSE]
    
    # Create correlation matrix based on selection patterns
    cor_matrix <- cor(t(stable_prob_matrix))
    
    # Create adjacency matrix (only strong correlations)
    adj_matrix <- ifelse(abs(cor_matrix) > 0.7 & cor_matrix != 1, 1, 0)
    
    # Create igraph object
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    V(g)$name <- rownames(stable_prob_matrix)
    V(g)$stability <- max_probs[stable_vars]
    
    # Set up layout
    layout <- layout_with_fr(g)
    
    # Create plot
    plot(g, 
         vertex.color = colorRampPalette(c("lightblue", "darkblue"))(100)[round(V(g)$stability * 100)],
         vertex.size = 10 + V(g)$stability * 20,
         vertex.label = V(g)$name,
         vertex.label.cex = 0.8,
         edge.color = "gray",
         edge.width = 2,
         layout = layout,
         main = paste0("Co-selection Network - Component ", component),
         sub = paste0("Variables with stability ≥ ", x$threshold, " (edges show co-selection > 0.7)"))
    
    # Add legend
    legend("topright", 
           legend = c("Low Stability", "High Stability"),
           fill = c("lightblue", "darkblue"),
           title = "Stability Level")
    
    return(invisible(g))
}


# Helper function for null coalescing
`%||%` <- function(x, y) {
    if(is.null(x)) y else x
}

# Helper function to extract sparsity info (if not available from other files)
sparsity_info <- function(x) {
    if(inherits(x, "sparse_o2pls")) {
        return(x$sparsity)
    } else {
        # Create basic sparsity info
        selected_X <- which(rowSums(abs(x$results$Xloading)) > 1e-10)
        selected_Y <- which(rowSums(abs(x$results$Yloading)) > 1e-10)
        
        return(list(
            selected_X = length(selected_X),
            total_X = nrow(x$results$Xloading),
            selected_Y = length(selected_Y), 
            total_Y = nrow(x$results$Yloading),
            sparsity_X = 1 - length(selected_X) / nrow(x$results$Xloading),
            sparsity_Y = 1 - length(selected_Y) / nrow(x$results$Yloading)
        ))
    }
}

#' Complete example workflow incorporating all new plotting functions
#' @examples
#' \dontrun{
#' # =============================================================================
#' # COMPREHENSIVE PLOTTING WORKFLOW
#' # =============================================================================
#' 
#' library(o2plsda)
#' library(ggplot2)
#' library(gridExtra)
#' 
#' # 1. Generate example data
#' set.seed(42)
#' n <- 100
#' p_X <- 200
#' p_Y <- 150
#' 
#' X <- matrix(rnorm(n * p_X), n, p_X)
#' Y <- matrix(rnorm(n * p_Y), n, p_Y)
#' 
#' # Add structured signal
#' signal <- matrix(rnorm(n * 3), n, 3)
#' X[, 1:30] <- signal %*% matrix(rnorm(3 * 30), 3, 30) + matrix(rnorm(n * 30, sd = 0.4), n, 30)
#' Y[, 1:25] <- signal %*% matrix(rnorm(3 * 25), 3, 25) + matrix(rnorm(n * 25, sd = 0.4), n, 25)
#' 
#' colnames(X) <- paste0("Gene_", 1:p_X)
#' colnames(Y) <- paste0("Metabolite_", 1:p_Y)
#' groups <- factor(rep(c("Control", "Disease", "Treatment"), length.out = n))
#' 
#' # 2. Parameter tuning
#' tune_results <- tune_o2pls(X, Y, nc_range = 1:3, nx_range = 0:2, ny_range = 0:2, 
#'                           folds = 3, nrepeat = 3, progressBar = FALSE)
#' 
#' # Plot tuning results
#' plot(tune_results, type = "heatmap")
#' plot(tune_results, type = "line")
#' plot(tune_results, type = "comparison")
#' 
#' # 3. Fit models
#' regular_fit <- o2pls(X, Y, nc = 2, nx = 1, ny = 1)
#' sparse_fit <- sparse_o2pls(X, Y, nc = 2, nx = 1, ny = 1, keepX = c(40, 30), keepY = c(30, 25))
#' 
#' # 4. Sparse O2PLS plots
#' plot(sparse_fit, type = "score", block = "X", group = groups, ellipse = TRUE)
#' plot(sparse_fit, type = "loading", block = "X", component = 1, top = 25)
#' plot(sparse_fit, type = "sparsity", block = "both")
#' plot(sparse_fit, type = "selection")
#' plot(sparse_fit, type = "biplot", component = c(1, 2), group = groups, top = 20)
#' 
#' # 5. Classification analysis
#' sparse_plsda_fit <- sparse_plsda(X, groups, nc = 2, keepX = c(35, 25))
#' class(sparse_plsda_fit) <- "sparse_plsda"  # Ensure correct class
#' 
#' plot(sparse_plsda_fit, type = "score", group = groups, ellipse = TRUE)
#' plot(sparse_plsda_fit, type = "loading", component = 1, top = 20)
#' plot(sparse_plsda_fit, type = "vip", component = 1, top = 25)
#' plot(sparse_plsda_fit, type = "selection")
#' 
#' # 6. Stability selection
#' stability_results <- stability_selection(X, groups, nc = 2, keepX_range = seq(20, 60, 20),
#'                                          n_bootstrap = 30, threshold = 0.7, 
#'                                          method = "sparse_plsda", parallel = FALSE)
#' 
#' plot(stability_results, type = "heatmap", component = 1, top_n = 40)
#' plot(stability_results, type = "barplot", component = 1, top_n = 25)
#' plot(stability_results, type = "threshold")
#' plot(stability_results, type = "summary")
#' 
#' # 7. Model comparison
#' compare_models(regular_fit, sparse_fit, type = "performance", 
#'               names = c("Regular O2PLS", "Sparse O2PLS"))
#' 
#' # 8. Create comprehensive dashboard
#' p1 <- plot(sparse_fit, type = "score", block = "X", group = groups) + ggtitle("A) Sparse O2PLS Scores")
#' p2 <- plot(sparse_fit, type = "loading", block = "X", component = 1, top = 15) + ggtitle("B) Top Loadings")
#' p3 <- plot(sparse_fit, type = "sparsity", block = "both") + ggtitle("C) Sparsity Pattern")
#' p4 <- plot(sparse_plsda_fit, type = "vip", component = 1, top = 15) + ggtitle("D) VIP Scores")
#' 
#' dashboard <- grid.arrange(p1, p2, p3, p4, ncol = 2,
#'                          top = textGrob("Comprehensive O2PLS Analysis Dashboard", 
#'                                        gp = gpar(fontsize = 16, fontface = "bold")))
#' 
#' print(dashboard)
#' }

# =============================================================================
# SPARSE O2PLS S4 CLASS PLOTTING FUNCTIONS
# =============================================================================

#' Plot method for SparseO2pls S4 objects
#' 
#' Creates various plots for sparse O2PLS models stored as S4 objects, including 
#' score plots, loading plots, sparsity patterns, and variable selection summaries.
#' 
#' @param x A SparseO2pls S4 object
#' @param type Character. Type of plot: "score", "loading", "sparsity", "selection", 
#'   "biplot", "contribution", "comparison", or "diagnostic"
#' @param component Integer or vector. Which component(s) to plot (default: c(1,2))
#' @param block Character. Which data block: "X", "Y", or "both" (default: "X")
#' @param group Factor. Grouping variable for score plots
#' @param top Integer. Number of top variables to show in loading plots (default: 20)
#' @param threshold Numeric. Threshold for variable selection display (default: 1e-10)
#' @param color Character vector. Custom colors for groups
#' @param ellipse Logical. Add confidence ellipses to score plots (default: TRUE)
#' @param title Character. Custom plot title
#' @param ... Additional plotting parameters
#' 
#' @return ggplot2 object or plotly object (for 3D plots)
#' 
#' @examples
#' # Generate example data for SparseO2pls plotting
#' set.seed(456)
#' n <- 120
#' p_X <- 180
#' p_Y <- 120
#' 
#' # Create structured data
#' X <- matrix(rnorm(n * p_X), n, p_X)
#' Y <- matrix(rnorm(n * p_Y), n, p_Y)
#' 
#' # Add correlated signal
#' signal <- matrix(rnorm(n * 4), n, 4)
#' X[, 1:40] <- signal %*% matrix(rnorm(4 * 40), 4, 40) + 
#'              matrix(rnorm(n * 40, sd = 0.4), n, 40)
#' Y[, 1:30] <- signal %*% matrix(rnorm(4 * 30), 4, 30) + 
#'              matrix(rnorm(n * 30, sd = 0.4), n, 30)
#' 
#' # Add meaningful names
#' colnames(X) <- paste0("Gene_", 1:p_X)
#' colnames(Y) <- paste0("Protein_", 1:p_Y)
#' rownames(X) <- rownames(Y) <- paste0("Sample_", 1:n)
#' 
#' # Create grouping variable
#' treatment_groups <- factor(rep(c("Control", "Low_Dose", "High_Dose", "Recovery"), 
#'                               each = n/4))
#' 
#' # Fit sparse O2PLS and convert to S4 object
#' sparse_fit_list <- sparse_o2pls(X, Y, nc = 3, nx = 1, ny = 1,
#'                                 keepX = c(50, 40, 30), keepY = c(35, 25, 20))
#' 
#' # Convert to S4 SparseO2pls object
#' sparse_s4 <- new("SparseO2pls",
#'                  X = X, Y = Y,
#'                  params = sparse_fit_list$params,
#'                  results = sparse_fit_list$results,
#'                  sparsity = sparse_fit_list$sparsity)
#' 
#' # =============================================================================
#' # SCORE PLOTS
#' # =============================================================================
#' 
#' # 2D score plot for X block with treatment groups
#' p1 <- plot(sparse_s4, type = "score", block = "X", 
#'            component = c(1, 2), group = treatment_groups, ellipse = TRUE)
#' print(p1)
#' 
#' # 2D score plot for Y block
#' p2 <- plot(sparse_s4, type = "score", block = "Y", 
#'            component = c(1, 3), group = treatment_groups,
#'            color = c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00"))
#' print(p2)
#' 
#' # 3D score plot
#' p3 <- plot(sparse_s4, type = "score", block = "X",
#'            component = c(1, 2, 3), group = treatment_groups)
#' print(p3)
#' 
#' # =============================================================================
#' # LOADING PLOTS
#' # =============================================================================
#' 
#' # Top gene loadings for component 1
#' p4 <- plot(sparse_s4, type = "loading", block = "X",
#'            component = 1, top = 30)
#' print(p4)
#' 
#' # Top protein loadings for component 2
#' p5 <- plot(sparse_s4, type = "loading", block = "Y",
#'            component = 2, top = 25,
#'            title = "Top Protein Loadings - Component 2")
#' print(p5)
#' 
#' # Multi-component loading comparison
#' p6 <- plot(sparse_s4, type = "loading", block = "X",
#'            component = c(1, 2, 3), top = 35)
#' print(p6)
#' 
#' # =============================================================================
#' # SPARSITY AND SELECTION PLOTS
#' # =============================================================================
#' 
#' # Sparsity pattern for both blocks
#' p7 <- plot(sparse_s4, type = "sparsity", block = "both")
#' print(p7)
#' 
#' # Variable selection summary
#' p8 <- plot(sparse_s4, type = "selection")
#' print(p8)
#' 
#' # =============================================================================
#' # ADVANCED PLOTS
#' # =============================================================================
#' 
#' # Biplot combining scores and loadings
#' p9 <- plot(sparse_s4, type = "biplot", 
#'            component = c(1, 2), group = treatment_groups, top = 20)
#' print(p9)
#' 
#' # Variable contribution analysis
#' p10 <- plot(sparse_s4, type = "contribution", block = "X",
#'             component = 1, top = 25)
#' print(p10)
#' 
#' # Diagnostic plots
#' p11 <- plot(sparse_s4, type = "diagnostic")
#' print(p11)
#' 
#' @export
plot.SparseO2pls <- function(x, type = "score", component = c(1, 2), 
                             block = "X", group = NULL, top = 20,
                             threshold = 1e-10, color = NULL, ellipse = TRUE,
                             title = NULL, ...) {
    
    # Load required libraries
    library(ggplot2)
    library(dplyr)
    
    # Validate inputs
    type <- match.arg(type, c("score", "loading", "sparsity", "selection", 
                              "biplot", "contribution", "comparison", "diagnostic"))
    block <- match.arg(block, c("X", "Y", "both"))
    
    # Generate appropriate plot based on type
    if(type == "score") {
        return(plot_sparse_s4_scores(x, component, block, group, color, ellipse, title, ...))
    } else if(type == "loading") {
        return(plot_sparse_s4_loadings(x, component, block, top, title, ...))
    } else if(type == "sparsity") {
        return(plot_sparse_s4_sparsity(x, block, title, ...))
    } else if(type == "selection") {
        return(plot_sparse_s4_selection(x, title, ...))
    } else if(type == "biplot") {
        return(plot_sparse_s4_biplot(x, component, group, top, color, title, ...))
    } else if(type == "contribution") {
        return(plot_sparse_s4_contribution(x, component, block, top, title, ...))
    } else if(type == "comparison") {
        return(plot_sparse_s4_comparison(x, title, ...))
    } else if(type == "diagnostic") {
        return(plot_sparse_s4_diagnostic(x, title, ...))
    }
}

# =============================================================================
# IMPLEMENTATION FUNCTIONS FOR SPARSE O2PLS S4 PLOTS
# =============================================================================

#' Plot SparseO2pls scores (S4 implementation)
#' @keywords internal
plot_sparse_s4_scores <- function(x, component, block, group, color, ellipse, title, ...) {
    # Extract scores from S4 object
    if(block == "X") {
        scores <- x@results$Xscore
        block_name <- "X Block"
        var_info <- paste0(length(x@sparsity$selected_X), "/", ncol(x@X), " variables selected")
    } else if(block == "Y") {
        scores <- x@results$Yscore
        block_name <- "Y Block"
        var_info <- paste0(length(x@sparsity$selected_Y), "/", ncol(x@Y), " variables selected")
    } else {
        stop("block must be 'X' or 'Y' for score plots")
    }
    
    # Handle different component specifications
    if(length(component) == 1) {
        # 1D score plot
        plot_data <- data.frame(
            Score = scores[, component],
            Index = 1:nrow(scores),
            Sample = rownames(scores) %||% paste0("Sample_", 1:nrow(scores))
        )
        
        if(!is.null(group)) {
            plot_data$Group <- group
            p <- ggplot(plot_data, aes(x = Index, y = Score, color = Group)) +
                geom_point(size = 3, alpha = 0.8) +
                geom_smooth(aes(group = Group), method = "loess", se = FALSE, alpha = 0.7)
        } else {
            p <- ggplot(plot_data, aes(x = Index, y = Score)) +
                geom_point(size = 3, alpha = 0.8) +
                geom_smooth(method = "loess", se = FALSE, alpha = 0.7)
        }
        
        p <- p + labs(x = "Sample Index", 
                      y = paste0("Component ", component),
                      title = title %||% paste0("Sparse O2PLS Scores - ", block_name),
                      subtitle = var_info)
        
    } else if(length(component) == 2) {
        # 2D score plot
        plot_data <- data.frame(
            PC1 = scores[, component[1]],
            PC2 = scores[, component[2]],
            Sample = rownames(scores) %||% paste0("Sample_", 1:nrow(scores))
        )
        
        if(!is.null(group)) {
            plot_data$Group <- group
            p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
                geom_point(size = 3, alpha = 0.8)
            
            if(ellipse) {
                p <- p + stat_ellipse(level = 0.95, linewidth = 1, alpha = 0.6)
            }
        } else {
            p <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
                geom_point(size = 3, alpha = 0.8)
        }
        
        # Calculate variance explained if available
        if(!is.null(x@results$component_variance)) {
            var_exp <- x@results$component_variance
            x_lab <- paste0("Component ", component[1], " (", 
                            round(var_exp$X[component[1]] * 100, 1), "%)")
            y_lab <- paste0("Component ", component[2], " (", 
                            round(var_exp$X[component[2]] * 100, 1), "%)")
        } else {
            x_lab <- paste0("Component ", component[1])
            y_lab <- paste0("Component ", component[2])
        }
        
        p <- p + labs(x = x_lab, y = y_lab,
                      title = title %||% paste0("Sparse O2PLS Scores - ", block_name),
                      subtitle = var_info)
        
    } else if(length(component) == 3) {
        # 3D score plot using plotly
        if(!requireNamespace("plotly", quietly = TRUE)) {
            stop("plotly package required for 3D plots")
        }
        library(plotly)
        
        plot_data <- data.frame(
            PC1 = scores[, component[1]],
            PC2 = scores[, component[2]], 
            PC3 = scores[, component[3]],
            Sample = rownames(scores) %||% paste0("Sample_", 1:nrow(scores))
        )
        
        if(!is.null(group)) {
            plot_data$Group <- group
            p <- plot_ly(plot_data, x = ~PC1, y = ~PC2, z = ~PC3,
                         color = ~Group, text = ~Sample,
                         type = "scatter3d", mode = "markers",
                         marker = list(size = 5))
        } else {
            p <- plot_ly(plot_data, x = ~PC1, y = ~PC2, z = ~PC3,
                         text = ~Sample, type = "scatter3d", mode = "markers",
                         marker = list(size = 5, color = "#2166AC"))
        }
        
        p <- p %>% layout(
            scene = list(
                xaxis = list(title = paste("Component", component[1])),
                yaxis = list(title = paste("Component", component[2])),
                zaxis = list(title = paste("Component", component[3]))
            ),
            title = title %||% paste0("3D Sparse O2PLS Scores - ", block_name)
        )
        
        return(p)
    }
    
    # Apply custom colors if provided
    if(!is.null(color) && !is.null(group)) {
        p <- p + scale_color_manual(values = color)
    } else if(!is.null(group)) {
        p <- p + scale_color_brewer(type = "qual", palette = "Set1")
    }
    
    p <- p + theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    return(p)
}

#' Plot SparseO2pls loadings (S4 implementation)
#' @keywords internal
plot_sparse_s4_loadings <- function(x, component, block, top, title, ...) {
    # Extract loadings from S4 object
    if(block == "X") {
        loadings <- x@results$Xloading
        block_name <- "X Block"
        var_names <- colnames(x@X) %||% paste0("X_", 1:ncol(x@X))
        selected_vars <- x@sparsity$selected_X
    } else if(block == "Y") {
        loadings <- x@results$Yloading
        block_name <- "Y Block"
        var_names <- colnames(x@Y) %||% paste0("Y_", 1:ncol(x@Y))
        selected_vars <- x@sparsity$selected_Y
    } else {
        stop("block must be 'X' or 'Y' for loading plots")
    }
    
    if(length(component) == 1) {
        # Single component loading plot
        loading_vals <- loadings[, component]
        
        # Only show selected variables (non-zero loadings)
        non_zero_idx <- which(abs(loading_vals) > 1e-10)
        if(length(non_zero_idx) == 0) {
            stop("No selected variables found for component ", component)
        }
        
        # Get top variables by absolute loading value
        top_indices <- non_zero_idx[order(abs(loading_vals[non_zero_idx]), decreasing = TRUE)]
        top_indices <- head(top_indices, min(top, length(top_indices)))
        
        plot_data <- data.frame(
            Variable = var_names[top_indices],
            Loading = loading_vals[top_indices],
            Absolute = abs(loading_vals[top_indices]),
            Selected = top_indices %in% selected_vars
        ) %>%
            arrange(desc(Absolute))
        
        # Reorder factor levels for plotting
        plot_data$Variable <- factor(plot_data$Variable, levels = plot_data$Variable)
        
        p <- ggplot(plot_data, aes(x = Variable, y = Loading)) +
            geom_col(aes(fill = interaction(Loading > 0, Selected)), alpha = 0.8) +
            scale_fill_manual(values = c("FALSE.FALSE" = "#CCCCCC", "TRUE.FALSE" = "#CCCCCC",
                                         "FALSE.TRUE" = "#D6604D", "TRUE.TRUE" = "#2166AC"),
                              name = "Type", 
                              labels = c("Negative (Selected)", "Positive (Selected)", 
                                         "Negative (Not Selected)", "Positive (Not Selected)")) +
            coord_flip() +
            labs(x = "Variables",
                 y = paste0("Loading Values - Component ", component),
                 title = title %||% paste0("Sparse O2PLS Loadings - ", block_name),
                 subtitle = paste0(length(non_zero_idx), " selected out of ", length(loading_vals), " variables")) +
            theme_classic(base_size = 12) +
            theme(axis.text.y = element_text(size = 8),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
        
    } else {
        # Multiple component comparison
        loading_data <- data.frame(
            Variable = rep(var_names, length(component)),
            Component = rep(paste0("Comp_", component), each = nrow(loadings)),
            Loading = as.vector(loadings[, component]),
            Selected = rep(1:nrow(loadings) %in% selected_vars, length(component))
        )
        
        # Filter to only selected variables and top by importance
        loading_data <- loading_data[loading_data$Selected, ]
        
        # Get top variables across all components
        var_importance <- aggregate(abs(Loading) ~ Variable, loading_data, max)
        top_vars <- var_importance[order(var_importance[,2], decreasing = TRUE), ][1:min(top, nrow(var_importance)), 1]
        loading_data <- loading_data[loading_data$Variable %in% top_vars, ]
        
        p <- ggplot(loading_data, aes(x = Variable, y = Loading, fill = Component)) +
            geom_col(position = "dodge", alpha = 0.8) +
            coord_flip() +
            scale_fill_brewer(type = "qual", palette = "Set2") +
            labs(x = "Variables",
                 y = "Loading Values",
                 title = title %||% paste0("Sparse O2PLS Loadings Comparison - ", block_name),
                 subtitle = paste0("Top ", length(top_vars), " selected variables across components")) +
            theme_classic(base_size = 12) +
            theme(axis.text.y = element_text(size = 8),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
    }
    
    return(p)
}

#' Plot SparseO2pls sparsity pattern (S4 implementation)
#' @keywords internal
plot_sparse_s4_sparsity <- function(x, block, title, ...) {
    library(reshape2)
    
    if(block == "both") {
        # Create combined sparsity plot
        X_loadings <- x@results$Xloading
        Y_loadings <- x@results$Yloading
        
        # Create binary selection matrices
        X_selected <- (abs(X_loadings) > 1e-10) * 1
        Y_selected <- (abs(Y_loadings) > 1e-10) * 1
        
        # Combine data
        sparsity_data <- rbind(
            data.frame(
                Variable = rep(1:nrow(X_selected), ncol(X_selected)),
                Component = rep(1:ncol(X_selected), each = nrow(X_selected)),
                Selected = as.vector(X_selected),
                Block = "X"
            ),
            data.frame(
                Variable = rep(1:nrow(Y_selected), ncol(Y_selected)),
                Component = rep(1:ncol(Y_selected), each = nrow(Y_selected)),
                Selected = as.vector(Y_selected),
                Block = "Y"
            )
        )
        
        # Add sparsity information
        sparsity_X <- round(x@sparsity$sparsity_X * 100, 1)
        sparsity_Y <- round(x@sparsity$sparsity_Y * 100, 1)
        
        p <- ggplot(sparsity_data, aes(x = Component, y = Variable, fill = factor(Selected))) +
            geom_tile(color = "white", linewidth = 0.1) +
            facet_wrap(~Block, scales = "free_y", 
                       labeller = labeller(Block = c("X" = paste0("X Block (", sparsity_X, "% sparse)"),
                                                     "Y" = paste0("Y Block (", sparsity_Y, "% sparse)")))) +
            scale_fill_manual(values = c("0" = "white", "1" = "#2166AC"),
                              name = "Selected", labels = c("No", "Yes")) +
            labs(x = "Component", y = "Variable Index",
                 title = title %||% "Variable Selection Pattern Across Components") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  strip.text = element_text(face = "bold"))
        
    } else {
        # Single block sparsity plot
        if(block == "X") {
            loadings <- x@results$Xloading
            block_name <- "X Block"
            sparsity_pct <- round(x@sparsity$sparsity_X * 100, 1)
            selected_vars <- x@sparsity$selected_X
        } else {
            loadings <- x@results$Yloading
            block_name <- "Y Block"
            sparsity_pct <- round(x@sparsity$sparsity_Y * 100, 1)
            selected_vars <- x@sparsity$selected_Y
        }
        
        selected <- (abs(loadings) > 1e-10) * 1
        sparsity_data <- melt(selected)
        colnames(sparsity_data) <- c("Variable", "Component", "Selected")
        
        p <- ggplot(sparsity_data, aes(x = Component, y = Variable, fill = factor(Selected))) +
            geom_tile(color = "white", linewidth = 0.1) +
            scale_fill_manual(values = c("0" = "white", "1" = "#2166AC"),
                              name = "Selected", labels = c("No", "Yes")) +
            labs(x = "Component", y = "Variable Index",
                 title = title %||% paste0("Variable Selection Pattern - ", block_name),
                 subtitle = paste0(sparsity_pct, "% sparsity (", length(selected_vars), " variables selected)")) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
    }
    
    return(p)
}

#' Plot SparseO2pls variable selection summary (S4 implementation)
#' @keywords internal
plot_sparse_s4_selection <- function(x, title, ...) {
    # Extract sparsity information from S4 object
    selected_X <- length(x@sparsity$selected_X)
    total_X <- ncol(x@X)
    selected_Y <- length(x@sparsity$selected_Y)
    total_Y <- ncol(x@Y)
    
    # Create summary data
    summary_data <- data.frame(
        Block = c("X Block", "Y Block"),
        Selected = c(selected_X, selected_Y),
        Total = c(total_X, total_Y),
        Percentage = c(round((selected_X/total_X)*100, 1), 
                       round((selected_Y/total_Y)*100, 1))
    )
    
    summary_data$Not_Selected <- summary_data$Total - summary_data$Selected
    
    # Reshape for stacked bar plot
    plot_data <- data.frame(
        Block = rep(summary_data$Block, 2),
        Count = c(summary_data$Selected, summary_data$Not_Selected),
        Type = rep(c("Selected", "Not Selected"), each = nrow(summary_data))
    )
    
    p <- ggplot(plot_data, aes(x = Block, y = Count, fill = Type)) +
        geom_col(position = "stack", alpha = 0.8) +
        geom_text(data = summary_data, 
                  aes(x = Block, y = Total/2, label = paste0(Selected, " vars\n(", Percentage, "%)")),
                  inherit.aes = FALSE, fontface = "bold", size = 4, color = "white") +
        scale_fill_manual(values = c("Selected" = "#2166AC", "Not Selected" = "#D1E5F0")) +
        labs(x = "Data Block", y = "Number of Variables",
             title = title %||% "Sparse O2PLS Variable Selection Summary",
             subtitle = paste0("Total sparsity: X = ", round(x@sparsity$sparsity_X*100, 1), 
                               "%, Y = ", round(x@sparsity$sparsity_Y*100, 1), "%")) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10),
              legend.position = "top")
    
    return(p)
}

#' Plot SparseO2pls diagnostic plots (S4 implementation)
#' @keywords internal
plot_sparse_s4_diagnostic <- function(x, title, ...) {
    # Create comprehensive diagnostic dashboard
    library(gridExtra)
    
    # 1. Variance explained plot
    if(!is.null(x@results$R2X) && !is.null(x@results$R2Y)) {
        var_data <- data.frame(
            Block = c("X", "Y"),
            Joint = c(x@results$R2X_joint %||% 0, x@results$R2Y_joint %||% 0),
            Orthogonal = c(x@results$R2X_orth %||% 0, x@results$R2Y_orth %||% 0),
            Total = c(x@results$R2X, x@results$R2Y)
        )
        
        var_long <- data.frame(
            Block = rep(var_data$Block, 2),
            Component = rep(c("Joint", "Orthogonal"), each = 2),
            Variance = c(var_data$Joint, var_data$Orthogonal)
        )
        
        p1 <- ggplot(var_long, aes(x = Block, y = Variance, fill = Component)) +
            geom_col(position = "stack", alpha = 0.8) +
            scale_fill_brewer(type = "qual", palette = "Set2") +
            labs(x = "Block", y = "R² Value", title = "Variance Explained") +
            theme_classic(base_size = 10) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    } else {
        p1 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Variance data not available") +
            theme_void() + ggtitle("Variance Explained")
    }
    
    # 2. Sparsity levels
    sparsity_data <- data.frame(
        Block = c("X", "Y"),
        Sparsity = c(x@sparsity$sparsity_X, x@sparsity$sparsity_Y) * 100
    )
    
    p2 <- ggplot(sparsity_data, aes(x = Block, y = Sparsity, fill = Block)) +
        geom_col(alpha = 0.8) +
        geom_text(aes(label = paste0(round(Sparsity, 1), "%")), vjust = -0.5, fontface = "bold") +
        scale_fill_manual(values = c("X" = "#2166AC", "Y" = "#D6604D")) +
        labs(x = "Block", y = "Sparsity (%)", title = "Sparsity Levels") +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "none")
    
    # 3. Component-wise selection
    nc <- x@params$nc
    comp_selection <- data.frame(
        Component = rep(1:nc, 2),
        Block = rep(c("X", "Y"), each = nc),
        Selected = c(x@params$keepX[1:nc], x@params$keepY[1:nc])
    )
    
    p3 <- ggplot(comp_selection, aes(x = factor(Component), y = Selected, fill = Block)) +
        geom_col(position = "dodge", alpha = 0.8) +
        scale_fill_manual(values = c("X" = "#2166AC", "Y" = "#D6604D")) +
        labs(x = "Component", y = "Variables Selected", title = "Selection per Component") +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # 4. Model parameters summary
    param_data <- data.frame(
        Parameter = c("Joint Components", "X-Orthogonal", "Y-Orthogonal"),
        Value = c(x@params$nc, x@params$nx %||% 0, x@params$ny %||% 0)
    )
    
    p4 <- ggplot(param_data, aes(x = Parameter, y = Value, fill = Parameter)) +
        geom_col(alpha = 0.8) +
        geom_text(aes(label = Value), vjust = -0.5, fontface = "bold") +
        scale_fill_brewer(type = "qual", palette = "Set3") +
        labs(x = "Parameter", y = "Count", title = "Model Configuration") +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine all diagnostic plots
    combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                                  top = textGrob(title %||% "Sparse O2PLS Diagnostic Summary", 
                                                 gp = gpar(fontsize = 14, fontface = "bold")))
    
    return(combined_plot)
}

#' Plot SparseO2pls biplot (S4 implementation)
#' @keywords internal
plot_sparse_s4_biplot <- function(x, component, group, top, color, title, ...) {
    # Extract scores and loadings
    scores <- x@results$Xscore[, component[1:2]]
    loadings <- x@results$Xloading[, component[1:2]]
    
    # Get selected variables only
    selected_idx <- x@sparsity$selected_X
    if(length(selected_idx) == 0) {
        stop("No selected variables for biplot")
    }
    
    # Get top loading variables from selected ones
    loading_magnitudes <- sqrt(rowSums(loadings[selected_idx, ]^2))
    top_selected_idx <- selected_idx[order(loading_magnitudes, decreasing = TRUE)][1:min(top, length(selected_idx))]
    
    # Scale loadings for plotting
    loading_scale <- max(abs(scores)) / max(abs(loadings[top_selected_idx, ])) * 0.8
    scaled_loadings <- loadings[top_selected_idx, ] * loading_scale
    
    # Create plot data
    score_data <- data.frame(
        PC1 = scores[, 1],
        PC2 = scores[, 2],
        Type = "Samples"
    )
    
    if(!is.null(group)) {
        score_data$Group <- group
    }
    
    loading_data <- data.frame(
        PC1 = scaled_loadings[, 1],
        PC2 = scaled_loadings[, 2],
        Variable = rownames(loadings)[top_selected_idx] %||% 
            colnames(x@X)[top_selected_idx] %||% 
            paste0("Var_", top_selected_idx),
        Type = "Variables"
    )
    
    # Create biplot
    p <- ggplot() +
        geom_point(data = score_data, aes(x = PC1, y = PC2, color = Group),
                   size = 3, alpha = 0.7) +
        geom_segment(data = loading_data, 
                     aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.1, "inches")),
                     color = "red", alpha = 0.8, linewidth = 1) +
        geom_text(data = loading_data,
                  aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Variable),
                  size = 3, color = "red", fontface = "bold") +
        labs(x = paste0("Component ", component[1]),
             y = paste0("Component ", component[2]),
             title = title %||% "Sparse O2PLS Biplot",
             subtitle = paste0("Top ", nrow(loading_data), " selected variables shown")) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    if(!is.null(color) && !is.null(group)) {
        p <- p + scale_color_manual(values = color)
    } else if(!is.null(group)) {
        p <- p + scale_color_brewer(type = "qual", palette = "Set1")
    }
    
    return(p)
}

#' Plot SparseO2pls variable contribution (S4 implementation)
#' @keywords internal
plot_sparse_s4_contribution <- function(x, component, block, top, title, ...) {
    # Extract appropriate data
    if(block == "X") {
        loadings <- x@results$Xloading
        scores <- x@results$Xscore
        var_names <- colnames(x@X) %||% paste0("X_", 1:ncol(x@X))
        selected_vars <- x@sparsity$selected_X
        block_name <- "X Block"
    } else if(block == "Y") {
        loadings <- x@results$Yloading
        scores <- x@results$Yscore
        var_names <- colnames(x@Y) %||% paste0("Y_", 1:ncol(x@Y))
        selected_vars <- x@sparsity$selected_Y
        block_name <- "Y Block"
    } else {
        stop("block must be 'X' or 'Y' for contribution plots")
    }
    
    if(length(component) > 1) {
        component <- component[1]
        warning("Multiple components specified, using first component only")
    }
    
    # Calculate variable contributions (loading * score variance)
    loading_vals <- loadings[, component]
    score_var <- var(scores[, component])
    contributions <- abs(loading_vals) * sqrt(score_var)
    
    # Focus on selected variables
    selected_contributions <- contributions[selected_vars]
    selected_names <- var_names[selected_vars]
    
    if(length(selected_contributions) == 0) {
        stop("No selected variables found for component ", component)
    }
    
    # Get top contributors
    top_indices <- order(selected_contributions, decreasing = TRUE)[1:min(top, length(selected_contributions))]
    
    plot_data <- data.frame(
        Variable = selected_names[top_indices],
        Contribution = selected_contributions[top_indices],
        Loading = loading_vals[selected_vars[top_indices]]
    ) %>%
        arrange(desc(Contribution))
    
    plot_data$Variable <- factor(plot_data$Variable, levels = plot_data$Variable)
    
    p <- ggplot(plot_data, aes(x = Variable, y = Contribution)) +
        geom_col(aes(fill = Loading > 0), alpha = 0.8) +
        scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
                          name = "Loading Direction", labels = c("Negative", "Positive")) +
        coord_flip() +
        labs(x = "Variables",
             y = paste0("Contribution to Component ", component),
             title = title %||% paste0("Variable Contributions - ", block_name),
             subtitle = paste0("Top ", nrow(plot_data), " contributing selected variables")) +
        theme_classic(base_size = 12) +
        theme(axis.text.y = element_text(size = 8),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    return(p)
}

#' Plot SparseO2pls comparison with regular O2PLS (S4 implementation)
#' @keywords internal
plot_sparse_s4_comparison <- function(x, title, ...) {
    # This function compares sparse vs regular O2PLS performance
    # For demonstration, we'll show the trade-off between sparsity and variance explained
    
    library(gridExtra)
    
    # 1. Variance explained comparison
    var_data <- data.frame(
        Block = c("X", "Y"),
        Sparse_R2 = c(x@results$R2X %||% 0, x@results$R2Y %||% 0),
        Sparsity = c(x@sparsity$sparsity_X, x@sparsity$sparsity_Y) * 100
    )
    
    p1 <- ggplot(var_data, aes(x = Sparsity, y = Sparse_R2, color = Block)) +
        geom_point(size = 4) +
        geom_text(aes(label = Block), vjust = -1, fontface = "bold") +
        scale_color_manual(values = c("X" = "#2166AC", "Y" = "#D6604D")) +
        labs(x = "Sparsity (%)", y = "R² Value",
             title = "Sparsity vs. Variance Explained Trade-off") +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "none")
    
    # 2. Selection efficiency across components
    nc <- x@params$nc
    if(nc > 1) {
        efficiency_data <- data.frame(
            Component = rep(1:nc, 2),
            Block = rep(c("X", "Y"), each = nc),
            KeepRatio = c(x@params$keepX[1:nc] / ncol(x@X), 
                          x@params$keepY[1:nc] / ncol(x@Y)) * 100
        )
        
        p2 <- ggplot(efficiency_data, aes(x = Component, y = KeepRatio, color = Block)) +
            geom_line(linewidth = 1.5) +
            geom_point(size = 3) +
            scale_color_manual(values = c("X" = "#2166AC", "Y" = "#D6604D")) +
            labs(x = "Component", y = "Variables Kept (%)",
                 title = "Selection Strategy Across Components") +
            theme_classic(base_size = 11) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    } else {
        p2 <- ggplot() + 
            annotate("text", x = 0.5, y = 0.5, 
                     label = "Single component model\nNo comparison available") +
            theme_void() + ggtitle("Selection Strategy")
    }
    
    # 3. Model complexity summary
    complexity_data <- data.frame(
        Measure = c("Total Variables", "Selected Variables", "Sparsity Level", "Components"),
        X_Block = c(ncol(x@X), length(x@sparsity$selected_X), 
                    round(x@sparsity$sparsity_X * 100, 1), x@params$nc),
        Y_Block = c(ncol(x@Y), length(x@sparsity$selected_Y), 
                    round(x@sparsity$sparsity_Y * 100, 1), x@params$nc)
    )
    
    complexity_long <- data.frame(
        Measure = rep(complexity_data$Measure, 2),
        Block = rep(c("X", "Y"), each = nrow(complexity_data)),
        Value = c(complexity_data$X_Block, complexity_data$Y_Block)
    )
    
    # Filter out sparsity for this plot (different scale)
    complexity_filtered <- complexity_long[complexity_long$Measure != "Sparsity Level", ]
    
    p3 <- ggplot(complexity_filtered, aes(x = Measure, y = Value, fill = Block)) +
        geom_col(position = "dodge", alpha = 0.8) +
        scale_fill_manual(values = c("X" = "#2166AC", "Y" = "#D6604D")) +
        labs(x = "Measure", y = "Count",
             title = "Model Complexity Summary") +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine plots
    if(nc > 1) {
        combined_plot <- grid.arrange(p1, p2, p3, ncol = 2, nrow = 2,
                                      layout_matrix = rbind(c(1, 2), c(3, 3)),
                                      top = textGrob(title %||% "Sparse O2PLS Performance Analysis", 
                                                     gp = gpar(fontsize = 14, fontface = "bold")))
    } else {
        combined_plot <- grid.arrange(p1, p3, ncol = 2,
                                      top = textGrob(title %||% "Sparse O2PLS Performance Analysis", 
                                                     gp = gpar(fontsize = 14, fontface = "bold")))
    }
    
    return(combined_plot)
}

