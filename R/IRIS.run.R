####################################################################################################
## Package : IRIS
## Version : 1.0.1
## Date    : 2024-01-04 09:10:08
## Title   : Accurate and Efficient Integrative Reference-Informed Spatial Domain Detection for Spatial Transcriptomics   
## Authors : Ying Ma
####################################################################################################

#' Accurate and Efficient Integrative Reference-Informed Spatial Domain Detection for Spatial Transcriptomics
#'
#' @param IRIS_object IRIS object create by the createIRISObject function
#' @param numCluster A numeric value indicating the number of clusters
#' @useDynLib IRIS
#' @import RcppArmadillo
#' @import rliger
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rdirichlet
#' @importFrom fields rdist
#' @importFrom SummarizedExperiment assays rowData colData
#' @return Returns a IRIS object with domain labels stored in IRIS_object@spatial_domain for each tissue slice.
#'
#' @export
#' 
IRIS_spatial <- function(IRIS_object,numCluster){
version = IRIS_object@project
ct.select = IRIS_object@internal_info$ct.select
ct.varname = IRIS_object@internal_info$ct.varname
sample.varname = IRIS_object@internal_info$sample.varname
sc_eset = IRIS_object@internal_info$sc_eset
sliceNames = IRIS_object@internal_info$sliceNames
countList = IRIS_object@countList
locationList = IRIS_object@locationList
K = numCluster
nSlices = length(countList)
if(version == "IRIS"){
cat(paste0("## Extract information from the scRNA-seq reference data! ...\n"))
Basis = createscRef(assays(sc_eset)$counts, colData(sc_eset),ct.select, ct.varname, sample.varname)
Basis = Basis[,match(ct.select,colnames(Basis))]
common_genes_sp = Reduce(intersect, lapply(countList, row.names))
common_genes = intersect(common_genes_sp,rownames(Basis))
##### rm the mito- genes 
common_genes  = common_genes[!(common_genes %in% common_genes[grep("mt-",common_genes)])]
info_genes = selectInfo(Basis,sc_eset,common_genes,ct.select,ct.varname) 
rm(sc_eset)
IRIS_object@internal_info$sc_eset <- NULL
gc()
}else if(version == "IRISfree"){
  markerList = IRIS_object@internal_info$markerList
  info_genes = IRIS_object@internal_info$marker
}
##### create the files that IRIS used for
countList = lapply(countList, function(x) { x[row.names(x) %in% info_genes,] })
countList = lapply(countList, function(x) { x[sort(rownames(x)),] })
countList = lapply(countList, function(x) { x[rowSums(x) > 0,] }) 
countList = lapply(countList, function(x) { x[,colSums(x) > 0] }) 
common_info_genes = Reduce(intersect, lapply(countList, row.names))
countList = lapply(countList, function(x) { x[row.names(x) %in% common_info_genes,] })
countList = lapply(countList, function(x) { x[rowSums(x) > 0,] }) 
countList = lapply(countList, function(x) { x[,colSums(x) > 0] }) 

###### normalize the count data set
countList = lapply(countList, function(x) { sweep(x,2,colSums(x),"/") })

###### check the location list
for(islice in 1:length(countList)){
  locationList[[islice]] = locationList[[islice]][rownames(locationList[[islice]]) %in% colnames(countList[[islice]]),]
  locationList[[islice]] = locationList[[islice]][match(colnames(countList[[islice]]),rownames(locationList[[islice]])),]
}

if(version == "IRIS"){
###### retain the common_info_genes in the B matrtix
B = Basis[rownames(Basis) %in% common_info_genes,]
B = B[order(rownames(B)),]
rm(Basis)
gc()
###### initialize V list
initVList = list()
for(islice in 1:nSlices){
set.seed(islice)
V = as.matrix(rdirichlet(ncol(countList[[islice]]), rep(1,ncol(B))))
V = t(V)
rownames(V) = colnames(B)
colnames(V) = colnames(countList[[islice]])
initVList[[islice]] <- t(V)
}
}else if(version == "IRISfree"){
numK = length(IRIS_object@internal_info$markerList)
###### Initialize B matrix
for(i in 1:length(countList)){
  colnames(countList[[i]]) = paste0("Slice",i,"_",colnames(countList[[i]]) )
}
ifnb_liger <- createLiger(countList,remove.missing = F,take.gene.union = T, verbose = F)
ifnb_liger <- normalize(ifnb_liger,remove.missing = F, verbose = F)
#ifnb_liger <- selectGenes(ifnb_liger,var.thresh = -1)
ifnb_liger@var.genes <- rownames(countList[[1]])
ifnb_liger <- scaleNotCenter(ifnb_liger, verbose = F)
ifnb_liger <- optimizeALS(ifnb_liger, k = numK, verbose = F)
initB = t(ifnb_liger@W)
###### remove the "_" in liger
for(i in 1:length(countList)){
  colnames(countList[[i]]) = gsub(paste0("Slice",i,"_"),"",colnames(countList[[i]]) )
}
rownames(initB) = rownames(countList[[1]])
colnames(initB) = paste0("CT",1:numK)
###### initialize V list
initVList = list()
for(islice in 1:length(ifnb_liger@H)){
set.seed(islice)
V = ifnb_liger@H[[islice]]
colnames(V) = paste0("CT",1:numK)
rownames(V) = colnames(countList[[islice]])
initVList[[islice]] <- V
}
}
###### create the adjacency matrix between spatial locations in each slice
cat(paste0("### Graph construction! ...\n"))
AList <- createA(locationList)

###### scale the data 
countList = lapply(countList, function(x) { x * 1e-01  / mean(x) }) ############ faster needed
if(version == "IRIS"){
B = B * 1e-01  / mean(B)
CommonCellTypes = colnames(B)
}else if(version == "IRISfree"){
  initB = initB * 1e-01  / mean(initB)
  CommonCellTypes = colnames(initB)
}

######  initialize the contatenated matrix
MatrixCombined = getConcated(initVList,locationList,sliceNames)
###### initialize the spatial domain labels for each tissue slice
kmeans = kmeansFunc_Initialize(MatrixCombined[,CommonCellTypes],K)
numCluster = length(unique(kmeans$kmeans))
MatrixCombined$kmeans_cluster = as.numeric(kmeans$kmeans)
MatrixCombined$kmeans_cluster = MatrixCombined$kmeans_cluster - 1
TotalClusters = unique(MatrixCombined$kmeans_cluster)[order(unique(MatrixCombined$kmeans_cluster))]

###### Initialize the Mu and STRlist
Mu_STR = getMu_and_STR(CommonCellTypes,TotalClusters,MatrixCombined,nSlices,sliceNames)
cat(paste0("### ",version," starts to run the domain detection! ...\n"))
###### Initialize IRIS
if(version == "IRIS"){
ResList= IRIS_ref_iter(countList,B,AList,initVList,Mu_STR$STRList,Mu_STR$Mu)
}else if(version == "IRISfree"){
  ResList= IRIS_Marker_iter(countList,initB,AList,initVList,Mu_STR$STRList,Mu_STR$Mu)
}
obj_old = ResList$Obj 
epsilon = 5e-04

for(iter in 1:1000){
###### add the column names and rownames to the V matrix of each tissue slice
for(islice in 1:length(ResList$VList)){
  rownames(ResList$VList[[islice]]) <- colnames(countList[[islice]])
  colnames(ResList$VList[[islice]]) <- CommonCellTypes
}
VList = ResList$VList
if(version == "IRISfree"){
  B = ResList$B
  colnames(B) = colnames(initB)
  rownames(B) = rownames(initB)
  B = B * 1e-01  / mean(B)
}
rm(ResList)
######  update the contatenated matrix
MatrixCombined = getConcated(VList,locationList,sliceNames)
if(iter == 1){
  centers_old = NULL
}
kmeans = kmeansFunc_Iter(MatrixCombined[,CommonCellTypes],K,centers_old)
numCluster = length(unique(kmeans$kmeans))
MatrixCombined$kmeans_cluster = as.numeric(kmeans$kmeans)
rm(kmeans)
MatrixCombined$kmeans_cluster = MatrixCombined$kmeans_cluster - 1
TotalClusters = unique(MatrixCombined$kmeans_cluster)[order(unique(MatrixCombined$kmeans_cluster))]
##### update Mu and STRLIst
Mu_STR = getMu_and_STR(CommonCellTypes,TotalClusters,MatrixCombined,nSlices,sliceNames)
centers_old = t(Mu_STR$Mu)
countList = lapply(countList, function(x) {as(x,"sparseMatrix")})
if(version == "IRIS"){
ResList= IRIS_ref_iter(countList,B,AList,initVList,Mu_STR$STRList,Mu_STR$Mu)
}else if(version == "IRISfree"){
  ResList= IRIS_Marker_iter(countList,B,AList,initVList,Mu_STR$STRList,Mu_STR$Mu)
}
obj = ResList$Obj 
logicalObjective = (obj_old - obj) * 2.0 / abs(obj + obj_old) < epsilon;
if(is.na(obj) || logicalObjective){
if(iter >= 30){  #### run at least 30 iterations
  break;
}else{
  obj_old = obj
}
}else if(obj >= obj_old){
  break;
  }else{
  obj_old = obj
}
}
cat(paste0("### ",version," Finished the domain detection! ...\n"))
###### store the current domain that makes the minimum objective function
spatialDomain = MatrixCombined[,c("Slice","spotName","x","y","kmeans_cluster")]
colnames(spatialDomain)[colnames(spatialDomain) == "kmeans_cluster"] <- paste0(version,"_domain")
if(version == "IRISfree"){
  IRIS_object@internal_info$B = B
}
IRIS_object@spatialDomain = spatialDomain
IRIS_object@IRIS_Prop = MatrixCombined[,c("Slice","spotName","x","y",CommonCellTypes)]
return(IRIS_object)
}


#' Construct the mean gene expression basis matrix (B), this is the faster version 
#'
#' @param countMat single-cell RNA-seq count data as the reference
#' @param metadata meta data containing the cell type annotations in the scRNA-seq reference data. The rows matched with the columns of countMat
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' @param sample.varname character,the name of the column in metaData that specifies the sample information. If NULL, we just use the whole as one sample.
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom wrMisc rowGrpMeans
#' @return Return a list of basis (B) matrix
#'
#'
createscRef <- function(countMat, metadata, ct.select = NULL, ct.varname, sample.varname = NULL){
  metadata <- metadata[metadata[, ct.varname] %in% ct.select,]
  ct.id <- droplevels(as.factor(metadata[, ct.varname]))
  sample.id <- as.character(metadata[, sample.varname])
  ct_sample.id <- paste(ct.id, sample.id, sep = "$*$")
  colSums_countMat <- colSums(countMat)
  colSums_countMat_Ct = aggregate(colSums_countMat ~ ct.id + sample.id, FUN = 'sum')
  colSums_countMat_Ct_wide = reshape(colSums_countMat_Ct, idvar = "sample.id", timevar = "ct.id", direction = "wide")
  colnames(colSums_countMat_Ct_wide) = gsub("colSums_countMat.","",colnames(colSums_countMat_Ct_wide))
  rownames(colSums_countMat_Ct_wide) = colSums_countMat_Ct_wide$sample.id
  colSums_countMat_Ct_wide$sample.id <- NULL
  tbl <- table(sample.id,ct.id)
  colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[,match(colnames(tbl),colnames(colSums_countMat_Ct_wide))]
  colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[match(rownames(tbl),rownames(colSums_countMat_Ct_wide)),]
  S_JK <- colSums_countMat_Ct_wide / tbl
  S_JK <- as.matrix(S_JK)
  S_JK[S_JK == 0] = NA
  S_JK[!is.finite(S_JK)] = NA
  S = colMeans(S_JK, na.rm = TRUE)
  rm(S_JK)
  rm(colSums_countMat_Ct_wide)
  gc()
  S = S[match(unique(ct.id),names(S))]
  Theta_S_rowSums <- rowsum(t(countMat),ct_sample.id)
  Theta_S_rowSums <- t(Theta_S_rowSums)
  Theta_S <- sweep(Theta_S_rowSums,2,colSums(Theta_S_rowSums),"/")
  grp <- sapply(strsplit(colnames(Theta_S),split="$*$",fixed = TRUE),"[",1)
  Theta = rowGrpMeans(Theta_S, grp = grp, na.rm = TRUE)
  Theta = Theta[,match(unique(ct.id),colnames(Theta))]
  S = S[match(colnames(Theta),names(S))]
  basis = sweep(Theta,2,S,"*")
  colnames(basis) = colnames(Theta)
  rownames(basis) = rownames(Theta)
  return(basis)
}

#' Select Informative Genes used in the deconvolution
#'
#' @param Basis Reference basis matrix.
#' @param sc_eset scRNAseq data along with meta data stored in the S4 class format (SingleCellExperiment).
#' @param commonGene common genes between scRNAseq count data and spatial resolved transcriptomics data.
#' @param ct.select vector of cell type names that you are interested in to deconvolute, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param ct.varname character, the name of the column in metaData that specifies the cell type annotation information
#' 
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats quantile
#' @return a vector of informative genes selected
#'
selectInfo <- function(Basis,sc_eset,commonGene,ct.select,ct.varname){
#### log2 mean fold change >0.5
gene1 = lapply(ct.select,function(ict){
rest = rowMeans(Basis[,colnames(Basis) != ict])
FC = log((Basis[,ict] + 1e-06)) - log((rest + 1e-06))
rownames(Basis)[FC > 1.25 & Basis[,ict] > 0]
})
gene1 = unique(unlist(gene1))
gene1 = intersect(gene1,commonGene)
counts = assays(sc_eset)$counts
counts = counts[rownames(counts) %in% gene1,]
sd_within = sapply(ct.select,function(ict){
  temp = counts[,colData(sc_eset)[,ct.varname] == ict]
  apply(temp,1,var) / apply(temp,1,mean)
  })
##### remove the outliers that have high dispersion across cell types
gene2 = rownames(sd_within)[apply(sd_within,1,mean,na.rm = T) < quantile(apply(sd_within,1,mean,na.rm = T),prob = 0.99,na.rm = T)]
return(gene2)
}

#' Construct the adjacency matrix for each tissue slice
#'
#' @param locationList a list of spatial location matrix
#' @import Matrix
#' @importFrom RANN nn2
#' @return Return a list of Ajacency matrix
#'
#'
createA <- function(locationList){
  nSlices = length(locationList)
  AList = list()
  for(islice in 1:nSlices){
    location = as.data.frame(locationList[[islice]])
    norm_cords = location[ ,c("x","y")]
    norm_cords$x = norm_cords$x - min(norm_cords$x)
    norm_cords$y = norm_cords$y - min(norm_cords$y)
    scaleFactor = max(norm_cords$x,norm_cords$y)
    norm_cords$x = norm_cords$x / scaleFactor
    norm_cords$y = norm_cords$y / scaleFactor
    rownames(norm_cords) <- rownames(location)
    ##### find 10 nearest neighbors
    ineibor = 11
    ##### find neiibors
    near_data = nn2(norm_cords[ ,1:2],k = ineibor)
    neibors = near_data$nn.idx
    neibors = neibors[,-1] #### delete itself
    Nmat = Matrix(0,nrow = nrow(neibors),ncol = nrow(neibors),sparse = TRUE)
    for(icol in 1:ncol(neibors)){
    edges = data.frame(i = 1:nrow(neibors), j = neibors[,icol])
    adjacency = sparseMatrix(i = as.integer(edges$i),
                              j = as.integer(edges$j),
                              x = 1,
                              dims = rep(nrow(neibors), 2),
                              use.last.ij = TRUE)
    Nmat = Nmat + adjacency
    }
    ##### only keep mutual neiibors
    Nmat= Nmat * t(Nmat)
    rownames(Nmat) = colnames(Nmat) = rownames(norm_cords)
    AList[[islice]] = Nmat
  }
  return(AList)

}

#' Construct the concatenated matrix from V across tissue_slices
#'
#' @param VList a list of V matrix across all tissue slices
#' @param locationList a list of location matrix matched with the spatial locations in the VList
#' @param sliceNames a vector of slice names corresponding to the input tissue slices
#' @return Return a list of Ajacency matrix
#'
#'
getConcated <- function(VList,locationList,sliceNames){
  MatrixCombined = NULL
  for(islice in 1:length(VList)){
   V_slice = VList[[islice]]
   V_slice = sweep(V_slice,1,rowSums(V_slice),"/")
   V_slice = as.data.frame(V_slice)
   V_slice$Slice = sliceNames[islice]
   V_slice$ID = paste0(V_slice$Slice,"_",rownames(V_slice))
   V_slice$spotName = rownames(V_slice)
   V_slice$x = locationList[[islice]]$x
   V_slice$y = locationList[[islice]]$y
   MatrixCombined = rbind(MatrixCombined,V_slice)
 }
  return(MatrixCombined)
}

#' get the Mu and spatial domain label list
#'
#' @param CommonCellTypes cell types in the reference matrix
#' @param TotalClusters number of total spatial odmains across all tissue slices
#' @param MatrixCombined a concatenated matrix of V matrix across tissue slices
#' @param nSlices number of total slices analyzed
#' @param sliceNames names of the tissue slices in the count list
#' @return Return a list of Mu and the spatial domain labels for each tissue slice
#'
#'
getMu_and_STR <- function(CommonCellTypes,TotalClusters,MatrixCombined,nSlices,sliceNames){
  STRList = list()
  Mu = Matrix(0,nrow = length(CommonCellTypes),ncol = length(TotalClusters))
  for(islice in 1:nSlices){
   Vslice = MatrixCombined[MatrixCombined$Slice == sliceNames[islice],]
   STRList[[islice]] = Vslice$kmeans_cluster
   Clusters = unique(Vslice$kmeans_cluster)[order(unique(Vslice$kmeans_cluster))]
   mu = Matrix(0,nrow = length(CommonCellTypes),ncol = length(TotalClusters))
   colnames(mu) = TotalClusters
   for(i in Clusters){
      mu[,colnames(mu) == i] <- colMeans(Vslice[Vslice$kmeans_cluster == i,CommonCellTypes])
   }
   Mu = Mu + mu
  }
  Mu = Mu / nSlices
  Mu = as.matrix(Mu)
  rownames(Mu) = CommonCellTypes
  return(list(Mu = Mu,STRList = STRList))
}

#' Initialize the spatial domain labels
#'
#' @param data Contatnated matrix from all tissue slices.
#' @param k number of spatial domains
#' 
#' @importFrom stats kmeans
#' @return a list of the kmeans clustering results
#'
kmeansFunc_Initialize <- function(data,k){
set.seed(12345678)
##### use more iterations to initialize the data
if(nrow(data)  < 300000){
  numStart = 100
}else{
    numStart = 1
}
cl <- suppressWarnings(try(kmeans(data, k,nstart = numStart,iter.max=100),silent = TRUE))

return(list(kmeans =cl$cluster))
}

#' Cluster the spatial domain labels in each iterartion
#'
#' @param data Contatnated matrix from all tissue slices.
#' @param k number of spatial domains
#' @param centers_old centers from last iteration
#' 
#' @importFrom stats kmeans
#' @return a list of the kmeans clustering results and the centers for this iteration
#'
kmeansFunc_Iter <- function(data,k,centers_old = NULL){
set.seed(12345678)
if(nrow(data)  < 300000){
  numStart = 10
}else{
    numStart = 1
}
if(is.null(centers_old)){
cl <- suppressWarnings(try(kmeans(data, k,nstart = 1,iter.max=numStart),silent = TRUE))
}else{
cl <- suppressWarnings(try(kmeans(data, centers_old,nstart = 2,iter.max=numStart),silent = TRUE))
}
if(class(cl) == "try-error"){
  cl <- suppressWarnings(try(kmeans(data, k,nstart = 1,iter.max=numStart),silent = TRUE))
}
resList = list(kmeans =cl$cluster,centers = cl$centers)
return(resList)
}

