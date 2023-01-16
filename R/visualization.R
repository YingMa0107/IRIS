####################################################################################################
## Package : IRIS
## Version : 1.0.1
## Date    : 2022-11-13 09:10:08
## Title   : Integrative and Reference-Informed Spatial Domain Detection for Spatial Transcriptomics  
## Authors : Ying Ma
####################################################################################################

####################################################################################################
#' Visualize the spatial domains
#'
#' @param IRIS_domain Data frame, spatial domains infered by IRIS across multiple slices. The column "IRIS_domain" contains the domain label. With each row representing one spatial location from each slice.
#' @param spatial_location Data frame, spatial location information for each spot in IRIS_domain in each tissue slice. The row names must be consistent with the row names in the data frame IRIS_domain. The data contains two columns, with the first column representing the x coordinates, with the second column representing the y coordinates
#' @param colors Vector of color names that you want to use, if NULL, we will use the default color. 
#' @param numCols Numeric, number of figure panels, it depends on the number of tissue slices you want to visualize.
#'
#' @import ggplot2 
#' @importFrom grDevices colorRampPalette
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'

IRIS.visualize.domain <- function(IRIS_domain,spatial_location,colors = NULL,numCols){
colorValBasic = c("#eec7b6","#e79d8d","#8ed7e6","#219ebc","#2A9D8F",
"#CC9A81","#81b29a","#A9CBCA","#f4f1de","#B18D97",
"#FAD471","#8e9aaf","#cbc0d3","#ffb703","#F7BEEC",
"#5965B0","#E16E6E")
if(sum(rownames(IRIS_domain) == rownames(spatial_location)) != nrow(IRIS_domain)){
    stop("The rownames of IRIS_domain data does not match with the rownames of spatial_location data")
}
if(is.null(colors)){
    if(length(unique(IRIS_domain$IRIS_domain)) > length(colorVal)){
    colors = colorRampPalette(colorValBasic)(length(unique(IRIS_domain$IRIS_domain)))
    }else{
        colors = colorValBasic
    }
}else{
    colors = colors
}

colnames(spatial_location) = c("plotx","ploty")
data_plot = cbind(IRIS_domain,spatial_location)
data_plot$IRIS_domain = factor(data_plot$IRIS_domain,levels = sort(unique(data_plot$IRIS_domain)))
data_plot$Slice = factor(data_plot$Slice,levels = sort(unique(data_plot$Slice)))
p = suppressMessages(ggplot(data_plot) + 
    geom_point(aes(x = plotx, y = ploty, colour = IRIS_domain),alpha = 1.0,size = 1.0,stroke = 0, shape = 16) +
    scale_colour_manual(values =  colors) + 
    #coord_fixed()+
    #facet_wrap(~Slice,nrow = floor(sqrt(Nslice)),scales = "free") + 
    facet_wrap(~Slice,scales = "free",ncol = numCols)+
    #geom_label(data = data, aes(x = 300, y = -80,label = ARI,fill = colorARI),hjust = "middle",colour = "#4D4E4D",size = 7,fontface = "bold")+
    #scale_fill_manual(values =  c("#e07a5f","#f4f1de")) + 
    theme(plot.margin = margin(0.35, 0.35, 0.35, 0.35, "cm"),
        panel.background = element_rect(colour = "white", fill="white"),
        plot.background = element_rect(colour = "white", fill="white"),
        panel.border = element_rect(colour = "grey39", fill=NA, size=0.5),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
          legend.text=element_text(size=12),
          legend.title=element_text(size = 13,face="bold"),
          legend.position = 'bottom',
          strip.text = element_text(size = 16,face="bold"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.key.size = unit(1.0, 'cm'))+
    guides(color = guide_legend(override.aes = list(size = 7,alpha=1),nrow=2),fill = 'none')+
    labs(color='IRIS domain'))
return(p)
}

#' Visualize the spatial distribution of cell type proportion
#'
#' @param proportion Data frame, cell type proportion estimated by IRIS in either original resolution or enhanced resolution.
#' @param spatial_location Data frame, spatial location information.
#' @param ct.visualize Vector of selected cell type names that are interested to visualize
#' @param colors Vector of color names that you want to use, if NULL, we will use the default color scale c("lightblue","lightyellow","red")
#' @param numCols Numeric, number of columns in the figure panel, it depends on the number of cell types you want to visualize.
#'
#' @import ggplot2 
#' @importFrom reshape2 melt
#' @return Returns a ggplot2 figure. 
#'
#' @export
#'

IRIS.visualize.eachProp <- function(proportion,spatial_location,ct.visualize = ct.visualize,colors = c("lightblue","lightyellow","red"),numCols){
if(is.null(colors)){
    colors = c("lightblue","lightyellow","red")
}else{
    colors = colors
}
res_IRIS = as.data.frame(proportion)
res_IRIS = res_IRIS[,order(colnames(res_IRIS))]
location = as.data.frame(spatial_location)
if(sum(rownames(res_IRIS)==rownames(location))!= nrow(res_IRIS)){
   stop("The rownames of proportion data does not match with the rownames of spatial location data")
}
ct.select = ct.visualize
res_IRIS = res_IRIS[,ct.select]
res_IRIS_scale = as.data.frame(apply(res_IRIS,2,function(x){
    (x - min(x)) / (max(x) - min(x))
} ))
res_IRIS_scale$x = as.numeric(location$x)
res_IRIS_scale$y = as.numeric(location$y)
mData = melt(res_IRIS_scale,id.vars = c("x","y"))
colnames(mData)[3] <- "Cell_Type"
b = c(0,1)
p = suppressMessages(ggplot(mData, aes(x, y)) + 
geom_point(aes(colour = value),alpha = 1.0,size = 1.0,stroke = 0, shape = 16) +
scale_color_gradientn(colours = colors) + 
#scale_color_viridis_c(option = 2)+
scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
facet_wrap(~Cell_Type,ncol = numCols)+ 
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    #legend.position=c(0.14,0.76),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank(),
    legend.title=element_text(size = 13,face="bold"),
    legend.text=element_text(size = 12),
    strip.text = element_text(size = 16,face="bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.key.size = unit(1.0, 'cm')))
return(p)
}



