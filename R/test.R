#' Create Venn Diagrams with Customizable Appearance
#'
#' Generates Venn diagrams to visualize overlaps between multiple sets of elements.
#' Supports 2-4 sets and provides extensive customization options for appearance.
#'
#' @param ... One or more vectors containing elements to compare (max 4)
#' @param out_name Character, output filename prefix (without extension), NULL means don't save
#' @param PPT Character, if not NULL, saves plot as PowerPoint with this name
#' @param label Character vector, custom names for each input set
#' @param label.cex Numeric, size of set labels, default is 8
#' @param text.cex Numeric, size of intersection counts, default is 8
#' @param fill.alpha Numeric, transparency of fill colors (0-1), default is 1
#' @param fill.col Character vector, fill colors for each set
#' @param stroke.alpha Numeric, transparency of stroke lines (0-1), default is 1
#' @param stroke.size Numeric, thickness of stroke lines, default is 1
#' @param mar Numeric, margin around plot, default is 1
#' @param show.elements Logical, whether to show elements instead of counts, default is FALSE
#' @param label.sep Character, separator for element labels when show.elements=TRUE, default is newline
#' @param auto.scale Logical, whether to auto-scale diagram, default is FALSE
#' @param x.lim Numeric vector, x-axis limits, NULL means automatic
#' @param y.lim Numeric vector, y-axis limits, NULL means automatic
#' @param add_border Logical, whether to add plot border, default is FALSE
#' @param hjust Numeric, horizontal justification of labels, NULL means automatic
#' @param stroke.linetype Character, line type for set boundaries, default is "longdash"
#' @param col.stroke Character, color for stroke lines, default is "black"
#' @param show.percentage Logical, whether to show percentages instead of counts, default is FALSE
#' @param w Numeric, output width in inches, default is 7
#' @param h Numeric, output height in inches, default is 5
#'
#' @return Returns a ggplot object containing the Venn diagram. For >4 sets, returns intersection matrix.
#'
#' @examples
#' set1 <- letters[1:10]
#' set2 <- letters[5:15]
#' venn(set1, set2, label = c("Group A", "Group B"))
#'
#' @import ggvenn
#' @import ggplot2
#' @import export
#' @export
venn <- function (..., out_name = NULL, PPT = NULL, label = NULL, 
    label.cex = 8, text.cex = 8, fill.alpha = 1, 
    fill.col = c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0"), 
    stroke.alpha = 1, stroke.size = 1, mar = 1, 
    show.elements = FALSE, label.sep = "\n", auto.scale = FALSE, 
    x.lim = NULL, y.lim = NULL, add_border = FALSE, hjust = NULL, stroke.linetype = "longdash",
    col.stroke = "black", show.percentage = FALSE, w = 7, h = 5) 
{
    data_list = list(...)
    if(is.null(label)){
        names(data_list) = getName(...)
    } else {
        names(data_list) = label
    }

    Len = length(data_list)
    vennM = matrix("", Len, Len, dimnames = list(getName(...), 
        getName(...)))
    for (i in 1:Len) {
        for (j in i:Len) {
            vennM[i, j] = length(intersect(data_list[[i]], data_list[[j]]))
        }
    }
    resM = apply(vennM, 1, as.numeric)
    dimnames(resM) = dimnames(vennM)
    if (Len > 4) {
        message("More than 4 objects, will not plot venn.")
        return(resM)
    }

    loadp(ggvenn)
    p <- suppressWarnings(
        ggvenn(data = data_list, 
               show_percentage = show.percentage,
               show_elements   = show.elements, 
               label_sep       = label.sep, 
               auto_scale      = auto.scale,
               set_name_size   = label.cex,
               text_size       = text.cex, 
               stroke_alpha    = stroke.alpha,
               stroke_size     = stroke.size, 
               fill_alpha      = fill.alpha, 
               fill_color      = fill.col,
               stroke_linetype = stroke.linetype
        ))
    
    if(is.null(x.lim)){
        if(length(data_list) == 4){
            nchar_left = nchar(names(data_list))[1]
            nchar_right = nchar(names(data_list))[4]
            x.lim = c(-1 - 0.25 * nchar_left, 1 + 0.25 * nchar_right)
        } else {
            x.lim = c(-1.8, 1.8)        
        }
    }

    if(is.null(y.lim)){
        if(length(data_list) == 2){
            y.lim = c(-1.2, 1.5)
        } else if(length(data_list) == 3){
            y.lim = c(-2.1, 2)        
        } else {
            y.lim = c(-1.8, 1.5)
        }
    }
   
    if(add_border) p <- p + theme_bw()
    if(!is.null(hjust)) p$layers[[3]]$data$hjust = hjust
    p <- p + scale_x_discrete(expand = c(0.05, 0.05))
    p <- p + xlim(x.lim) + ylim(y.lim)

    if(!is.null(out_name)) save.graph(p, file = out_name, w, h)

    if (!is.null(PPT)) {
        require(export)
        graph2ppt(p, file = paste0(out_name, ".pptx"), width = w, height = h)
    }
    return(p)
}

#' Perform Principal Component Analysis (PCA) and Visualization
#'
#' Conducts PCA on an expression matrix and generates a scatter plot of PC1 vs PC2
#' with optional confidence ellipses and extensive customization options.
#'
#' @param expM Numeric matrix where rows represent genes/features and columns represent samples
#' @param Group Character vector specifying group information for each sample
#' @param out_name Character, output filename (without extension), NULL means don't save
#' @param text.size Numeric, text size in plot, default is 20
#' @param pt.size Numeric, point size, default is 1
#' @param ellipse Logical, whether to add confidence ellipses, default is TRUE
#' @param w Numeric, output width in inches, default is 5
#' @param h Numeric, output height in inches, default is 5
#' @param legend.ncol Integer, number of columns in legend, NULL means auto
#' @param color Character vector, group colors, default is c("firebrick3", "steelblue")
#'
#' @return A ggplot2 object containing PCA visualization
#'
#' @examples
#' data <- matrix(rnorm(1000), nrow=100)
#' groups <- rep(c("A","B"), each=5)
#' p <- run_PCA(data, groups)
#' print(p)
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @export
run_PCA <- function(expM, Group, out_name = NULL, text.size = 20, pt.size = 1, 
    ellipse = TRUE, w = 5, h = 5, legend.ncol = NULL, color = c("firebrick3", "steelblue")
    ){
    expM = expM[apply(expM, 1, sd) != 0, ]
    pca = prcomp(t(expM), scale = TRUE)
    pca.data = data.frame(pca$x)

    loadp(ggplot2)
    p <- ggplot(pca.data, aes(x = PC1, y = PC2, color = Group))
    p <- p + geom_point(size = pt.size)
    p <- p + geom_hline(yintercept = 0, lty = 4, col = "black", lwd = 0.6)
    p <- p + geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.6)
    p <- p + guides(color = guide_legend(override.aes = list(size = 2)))
    p <- p + gg_style(text.size, theme = "bw") + theme(legend.position = "top")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_color_manual(values = c("firebrick3", "steelblue"))
    if(ellipse) p <- p + stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95)
    if(!is.null(legend.ncol))
        p <- p + guides(color = guide_legend(ncol = legend.ncol))
    if(!is.null(out_name)) save.graph(p, file = out_name, w, h)
    return(p)
}

#' Perform Principal Component Analysis (PCA) and Visualization
#'
#' Conducts PCA on an expression matrix and generates a scatter plot of PC1 vs PC2
#' with optional confidence ellipses and extensive customization options.
#'
#' @export
gg_style <- function(title.size   = 10,
                    title.face = "bold",
                    title.hjust = 0.5,
                    title.vjust = 0,
                    x.title.size = NA,
                    x.title.face = NA,
                    x.title.vjust = 0.5,
                    x.text.size  = NA,
                    x.text.angle = 0,
                    x.text.hjust = 0.5,
                    x.text.vjust = 1,
                    y.title.size = NA,
                    y.title.vjust = 2,
                    y.text.size  = NA,
                    face         = "plain",
                    y.title.face = NA,
                    y.text.angle = 0,
                    y.text.hjust = 1,
                    legend.position = "right",
                    plot.margin = c(5.5, 5.5, 5.5, 15),
                    color = NULL,
                    fill  = NULL,
                    theme = "classic"
){
    loadp(ggplot2)
    checkNA <- function(x, y, ratio = NA){
        
        if(is.na(x))
            if(is.na(ratio)){
                x = y
            } else{
                x = y * ratio
            }
        return(x)
    }
    
    x.title.size = checkNA(x.title.size, title.size, 0.8)
    y.title.size = checkNA(y.title.size, title.size, 0.8)
    x.text.size  = checkNA(x.text.size, title.size, 0.6)
    y.text.size  = checkNA(y.text.size, title.size, 0.6)
    x.title.face = checkNA(x.title.face, face)
    y.title.face = checkNA(y.title.face, face)
    if(x.text.angle){
        if(x.text.hjust != 45){
            x.text.hjust = 1
        }
    }
    if(y.text.angle)
        y.text.hjust = 1
    if(x.text.angle == 90){
        x.text.vjust = 0.5
    }

    mytheme = switch(theme,
            "bw" = theme_bw(),
            "classic"  = theme_classic(),
            "box" = theme_bw() + theme(panel.grid =element_blank())
            )
    p <- mytheme
    p <- p + theme(
               plot.title   = element_text(color = "black", size = title.size, face = title.face, hjust = title.hjust, vjust = title.vjust),
               axis.title.x = element_text(color = "black", size = x.title.size, face = x.title.face, vjust = x.title.vjust),
               axis.title.y = element_text(color = "black", size = y.title.size, face = y.title.face, vjust = y.title.vjust),
               axis.text.x  = element_text(color = "black", size = x.text.size, angle = x.text.angle, hjust = x.text.hjust, vjust = x.text.vjust),
               axis.text.y  = element_text(color = "black", size = y.text.size, angle = y.text.angle, hjust = y.text.hjust),
               plot.margin  = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt"),
               legend.title    = element_text(size = unit(x.title.size, "pt")),
               legend.key.size = unit(x.text.size, "pt"),
               legend.text     = element_text(size = unit(x.text.size, "pt")),
               legend.position = legend.position)
    if(theme == "blank"){
        p <- p + theme(axis.line = element_blank(), 
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), plot.background = element_blank())
    }
    if(!is.null(color)) p <- p + scale_color_manual(values = color)
    if(!is.null(fill))  p <- p + scale_fill_manual(values = fill)
    return(p)
} 

