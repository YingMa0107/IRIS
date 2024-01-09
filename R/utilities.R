####################################################################################################
## Package : IRIS
## Version : 1.0.1
## Date    : 2022-11-13 09:10:08
## Title   : Accurate and Scalable Spatial Domain Detection via Integrated Reference-Informed Segmentation for Spatial transcriptomics 
## Authors : Ying Ma
####################################################################################################
#' Each IRIS object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot countList A list of spatial count data matrix. Each spatial count data with rows as the genes and column as measured spatial locations
#' @slot locationList A list of spatial location matrix. Each spatial location with rows as the measured spatial locations matched with those in the spatial_countMat_list and columns as theccccczzzz
#' @slot spatialDomain A meta data frame of spatial domain labels detected by IRIS.
#' @slot project The name of the project, default is "IRIS".
#' @slot IRIS_Prop A meta data frame of concatenated V matrix estimated by IRIS.
#' @slot internal_info The neccessary information that are used by IRIS computation.
#'
setClass("IRIS", 
	slots = list(
	countList = "list",
	locationList = "list",
	spatialDomain = "data.frame",
	project = "character",
	IRIS_Prop = "data.frame",
	internal_info = "list")
	)
#' Each IRISfree object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot countList A list of spatial count data matrix. Each spatial count data with rows as the genes and column as measured spatial locations
#' @slot locationList A list of spatial location matrix. Each spatial location with rows as the measured spatial locations matched with those in the spatial_countMat_list and columns as theccccczzzz
#' @slot spatialDomain A meta data frame of spatial domain labels detected by IRIS.
#' @slot project The name of the project, default is "IRIS".
#' @slot IRIS_Prop A meta data frame of concatenated V matrix estimated by IRIS.
#' @slot internal_info The neccessary information that are used by IRIS computation.
#'
setClass("IRISfree", 
	slots = list(
	countList = "list",
	locationList = "list",
	spatialDomain = "data.frame",
	project = "character",
	IRIS_Prop = "data.frame",
	internal_info = "list")
	)



#' Quality control of scRNA-seq count data
#'
#' @param sc_count Raw scRNAseq count data, each column is a cell and each row is a gene.
#' @param sc_meta data frame, metaData with "ct.varname" specify the cell type annotation information and "sample.varname" specify the sample information
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#' @param min.cells numeric, we filtered out the non-expressed cells.
#' @param min.genes numeric we filtered out the non-expressed genes
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return Return the filtered scRNA-seq data and meta data stored in a S4 class (SingleCellExperiment)
#'
#' @export
#'
sc_QC <- function(sc_count,sc_meta,ct.varname,sample.varname = NULL, min.cells = 0,min.genes = 0){
    if(is(sc_count,"matrix")){
	sc_count  <- as(as.matrix(sc_count), "sparseMatrix")
	}else if(is(sc_count,"vector")){
		sc_count  <- as(t(as.matrix(sc_count)), "sparseMatrix")
		}else if(is(sc_count,"sparseMatrix")){
			sc_count  <- sc_count
			}else{
				stop("scRNASeq counts has to be of following forms: vector,matrix or sparseMatrix")
}
if (missing(x = sc_count)) {
	stop("Please provide scRNASeq count data")
	} else if (is.null(sample.varname) || missing(sample.varname)) {
		sample.varname = "Sample"
		sc_meta = as.data.frame(sc_meta)
		sc_meta$sampleID = "Sample"
		} else if (any(rownames(x = sc_count) == '')) {
			stop("Feature names of sc_count matrix cannot be empty", call. = FALSE)
			} else if(sum(rownames(sc_meta) == colnames(sc_count)) != ncol(sc_count)){
				stop("Cell name in sc_count count data does not match with the rownames of sc_meta")
				} else if(ncol(sc_count)!=nrow(sc_meta)){
					stop("The number of cells in sc_count and sc_meta should be consistent! (sc_count -- p x c; sc_meta -- c x 2 or c x m)")
}
if (is.null(ct.varname)){
	stop("Please provide the column name indicating the cell type information in the sc_count")
}
ct.select <- unique(sc_meta[, ct.varname])
if (is.null(ct.select)){
	stop("No cell type information in the column of ct.varname in the sc_meta")
}
# Filter based on min.features
if (min.genes >= 0) {
    nfeatures <- colSums(x = sc_count )
    keep = which(nfeatures > min.genes)
    sc_count <- sc_count[, keep]
    sc_meta <- sc_meta[keep,]
}
# filter genes on the number of cells expressing
if (min.cells >= 0) {
    num.cells <- rowSums(x = sc_count > 0)
    sc_count <- sc_count[which(x = num.cells > min.cells), ]
}
fdata = as.data.frame(rownames(sc_count))
rownames(fdata) = rownames(sc_count)
ct.select <- as.character(ct.select[!is.na(ct.select)])
####### at least two cells
ct.select = names(table(sc_meta[,ct.varname]))[table(sc_meta[,ct.varname]) > 1]
keepCell = as.character(sc_meta[,ct.varname]) %in% ct.select
sc_count = sc_count[,keepCell]
sc_meta = sc_meta[keepCell,]
keepGene = rowSums(sc_count) > 0
fdata = as.data.frame(fdata[keepGene,])
sc_count = sc_count[keepGene,]
sce <- SingleCellExperiment(list(counts=sc_count),
colData=as.data.frame(sc_meta),
rowData=as.data.frame(fdata))
return(sce)
}

#' Create the IRIS object
#'
#' @param spatial_countMat_list A list of raw spatial resolved transcriptomics data. In each spatial transcriptomics data, it is a sparse matrix with each column is a spatial location, and each row is a gene. 
#' @param spatial_location_list A list of spatial location data frames, with two columns representing the x and y coordinates of the spatial location (The column names must contain "x" and "y"). The rownames of this data frame should match eaxctly with the columns of the corresponding spatial_count in the spatial_countMat_list.
#' @param sc_count Reference scRNA-seq count data, each column is a cell and each row is a gene.
#' @param sc_meta A data frame of meta data providing information on each cell in the sc_count, with each row representing the cell type and/or sample information of a specific cell. The row names of this data frame should match exactly with the column names of the sc_count data
#' @param ct.varname character, the name of the column in metaData that provides the cell type annotation information for IRIS to perform informed domain detection.
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#' @param markerList a list of marker genes, with each element of the list being the vector of cell type specific marker genes, default is NULL
#' @param version A character indicating whether to use IRIS, IRISfree or IRISsingle. The default version is "IRIS". We highly recommend using the original IRIS version.
#' @param minCountGene Minimum counts for each gene in the spatial count data
#' @param minCountSpot Minimum counts for each spatial location in the spatial count data

#'
#' @importFrom SummarizedExperiment assays rowData
#' @import methods
#' @return Returns IRIS object with filtered spatial count list, location list and reference scRNA-seq dataset.
#'
#' @export
#'
createIRISObject <- function(spatial_countMat_list,spatial_location_list,sc_count,sc_meta,ct.varname,sample.varname,markerList = NULL,version = "IRIS",minCountGene = 100,minCountSpot =5){  

#### QC on spatial dataset
cat(paste0("## QC on spatially-resolved dataset! ...\n"))
if (length(spatial_countMat_list) != length(spatial_location_list)) {
			stop("The length of spatial_countMat_list does not equal to the length of spatial_location_list", call. = FALSE)
}
Nslices = length(spatial_countMat_list)
sliceNames = names(spatial_countMat_list)
if(is.null(sliceNames)){
	if(length(spatial_countMat_list) >= 2){
	cat(paste0("## There are no names for the tissue slice, we will use Slice 1, ..., Slice",length(sliceNames)," to represent! ...\n"))
	}
	sliceNames = paste0("Slice",1:length(spatial_countMat_list))
}
countList = list()
locationList = list()
for(islice in 1:Nslices){
	spatial_count = spatial_countMat_list[[islice]]
	spatial_location = spatial_location_list[[islice]]
	if(is(spatial_count,"matrix")){
		spatial_count  <- as(as.matrix(spatial_count), "sparseMatrix")
	}else if(is(spatial_count,"vector")){
		spatial_count  <- as(t(as.matrix(spatial_count)), "sparseMatrix")
		}else if(is(spatial_count,"sparseMatrix")){
			spatial_count <- spatial_count
			}else{
				stop(paste0("The ",islice,"-th spatial resolved transcriptomic counts has to be of following forms: vector,matrix or sparseMatrix"))
			}
	if (any(rownames(x = spatial_count) == '')) {
			stop(paste0("Gene names of the ",islice,"-th spatial count matrix cannot be empty"), call. = FALSE)
	}
	if(is.null(spatial_location)){
	stop(paste0("Please provide the matched spatial location data frame for the ",islice,"-th spatial count matrix"))
    }
    if(ncol(spatial_count)!=nrow(spatial_location)){
	stop(paste0("The number of the columns in the ",islice,"-th spatial_count and the number of rows in the ",islice,"-th spatial_location should be consistent! (spatial_count -- p x n_slice; spatial_location -- n_slice x 2)"))
	}# end fi
	## check data order should consistent
	if(!identical(colnames(spatial_count), rownames(spatial_location))){
	stop(paste0("The column names of the ",islice,"-th spatial_count and row names of the ",islice,"-th spatial_location should be consistent! (spatial_count -- p x n; spatial_location -- n x 2)"))
	}# end fi
	spatial_count = spatial_count[rowSums(spatial_count > 0) > minCountSpot,]
    spatial_count = spatial_count[,(colSums(spatial_count) >= minCountGene)]
    spatial_location = spatial_location[rownames(spatial_location) %in% colnames(spatial_count),]
    spatial_location = spatial_location[match(colnames(spatial_count),rownames(spatial_location)),] 
    countList[[islice]] <- spatial_count 
    locationList[[islice]] <- spatial_location
}
names(countList) = names(spatial_countMat_list)
rm(spatial_countMat_list)
rm(spatial_location_list)
gc()
common_names = Reduce(intersect, lapply(countList, row.names))
if (length(common_names) == 0) {
			stop("There are no common gene names in the spatial_countMat_list", call. = FALSE)
}
## QC on scRNA-seq count
if(version == "IRIS"){
	#### create SingleCellExperiment Object and QC on the scRNA-seq reference data
	cat(paste0("## QC on scRNA-seq dataset! ...\n"))
	sc_eset = sc_QC(sc_count,sc_meta,ct.varname,sample.varname)
	ct.select = as.character(unique(colData(sc_eset)[,ct.varname]))
	common_genes = markerList = NULL
	rm(sc_count)
	rm(sc_meta)
	gc()
}else if(version == "IRISfree"){
	if(is.null(markerList)){
		stop("Please provide marker genes for at least two cell types in the markerList")
	}else if(!is.null(markerList) & length(markerList) == 1){
		stop("The length of markerList is 1, please provide marker genes for at least two cell types in the markerList")
	}else if(!is.null(markerList) & length(markerList) > 1){
	numK = length(markerList)
	marker = unique(unlist(markerList))
	#marker = toupper(marker)
	cat(paste0("## Number of unique marker genes: ",length(marker)," for ",numK," cell types ...\n"))
	commonGene = intersect(toupper(common_names),toupper(marker))
	#### remove mitochondrial and ribosomal genes
	commonGene  = commonGene[!(commonGene %in% commonGene[grep("mt-",commonGene)])]
	if(length(commonGene) < numK * 10){
		cat(paste0("## STOP! The average number of unique marker genes for each cell type is less than 10 ...\n"))
	}
	marker = unique(unlist(markerList))
	numK = length(markerList)
	common_genes = intersect(common_names,marker)
	##### rm the mito- genes 
	common_genes  = common_genes[!(common_genes %in% common_genes[grep("mt-",common_genes)])]
	sc_eset = sample.varname = ct.select = ct.varname = NULL

}
}
object <- new(
		Class = version,
		countList = countList,
		locationList = locationList,
		project = version,
		internal_info = list(ct.varname = ct.varname,ct.select = ct.select,sample.varname = sample.varname,sc_eset = sc_eset,sliceNames = sliceNames, markerList = markerList, marker = common_genes)
		)
gc()
return(object)
}


