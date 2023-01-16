#include <iostream>
#include <fstream>
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// calculate objective function
double calculateObject(List &XList,List &VList,List &AList,List &StrList,arma::mat &B,arma::mat &Mu,double beta,double lambda){
    double obj = 0;
    int nSlice = XList.size();
    for(int islice = 0; islice < nSlice; ++islice){
        double objt = 0;
        arma::sp_mat Xt = XList(islice);
        arma::mat Vt = VList(islice);
        arma::sp_mat At = AList(islice);
        arma::sp_mat colsum_At = arma::sum(At, 1);
        arma::sp_mat Dt((int)Vt.n_rows,(int)Vt.n_rows);
        Dt.diag() = colsum_At;
        arma::sp_mat Lt = Dt - At;
        double Obj_NMF = norm(Xt - B * Vt.t(),"fro");
        arma::vec vecStr_t = StrList(islice);
        arma::vec UniqueStr_t = arma::unique(vecStr_t);
        double Obj_V = 0.0;
        for(int istr : UniqueStr_t){
            uvec indexStr = find(vecStr_t == istr);
            int NspotsInStr = indexStr.size();
            arma::vec vecOnetnr = ones<vec>(NspotsInStr);
            arma::vec VconstraintStr = Vt.rows(indexStr).t() * vecOnetnr / NspotsInStr - Mu.col(istr);
            double Obj_VStr = norm(VconstraintStr,2); // norm is the sqrt norm
            Obj_V = Obj_V + Obj_VStr * Obj_VStr;

        }
        //double Obj_NMF = norm(mainNMF,"fro");
        objt = Obj_NMF * Obj_NMF + lambda * trace(Vt.t() * Lt * Vt);;
        obj = obj + objt;
    }
    return obj;

}


arma::mat calc_VtDt(arma::sp_mat &colsum_At, arma::mat &Vt, arma::uvec &indexStr){
    arma::vec colsum_At_dense = vec(colsum_At.col(0));
    arma::mat res = Vt.rows(indexStr);
    res.each_col() %= colsum_At_dense.elem(indexStr);
    return(res);
}


//*******************************************************************//
// IRIS for Accurate and Scalable Spatial Domain Detection via Integrated Reference-Informed Segmentation for Spatial transcriptomics //
//*******************************************************************//
//' IRIS_ref function 
//' @param XList The input list of normalized spatial data
//' @param BIn The input list of cell type specific basis matrix B
//' @param AList The constructed Ajacency matrix for each tissue slice
//' @param VList List of initial matrix of cell type compositions V
//' @param StrList List of initial spatial domains for each tissue slice
//' @param MuIn Initial matrix of mean cell type proportion across slices
//'
//' @return A list
//'
//' @export
// [[Rcpp::export]]
SEXP IRIS_ref_iter(Rcpp::List &XList, SEXP BIn,Rcpp::List &AList,Rcpp::List VList, Rcpp::List StrList,SEXP MuIn)
{    
    try {
        arma::mat B = as<mat>(BIn); 
        double beta = 1000;
        double lambda = 2000;
        arma::mat Mu = as<mat>(MuIn);// with each row represents the cell type and each column represents the structure
        int nSlice = XList.size();
        // initialize the objective function
        double obj = 0;         
        // iterate through the slices
        arma::mat Mu_Slice = zeros<mat>((int)Mu.n_rows,(int)Mu.n_cols);
        for(int islice = 0; islice < nSlice; ++islice){
            arma::sp_mat Xt = XList(islice);
            arma::mat Vt = VList(islice);
            int Nt = (int)Xt.n_cols; // number of spatial sample points for the t-th slice
            arma::vec vecOnet = ones<vec>(Nt);
            arma::sp_mat At = AList(islice);
            arma::sp_mat colsum_At = arma::sum(At, 1);
            arma::sp_mat Dt((int)Vt.n_rows,(int)Vt.n_rows);
            Dt.diag() = colsum_At;
            // update Vt
            arma::vec vecStr_t = StrList(islice);
            arma::vec UniqueStr_t = unique(vecStr_t);
            arma::mat Mu_Str = zeros<mat>((int)Mu.n_rows,(int)Mu.n_cols);
            for(int istr : UniqueStr_t){
                uvec indexStr = find(vecStr_t == istr);
                uvec notIndexStr = find(vecStr_t != istr);
                int NspotsInStr = indexStr.size();
                arma::vec vecOnetnr = ones<vec>(NspotsInStr);
                arma::mat vecOnetnrMat = ones<mat>(NspotsInStr,NspotsInStr);
                arma::mat VtDt = calc_VtDt(colsum_At,Vt,indexStr);
                arma::mat VtAt = At.cols(indexStr).t() * Vt;
                arma::mat nomUpdateVtStr = (Xt.cols(indexStr).t() * B + beta / NspotsInStr * (vecOnetnr * Mu.col(istr).t()) + lambda * VtAt);
                arma::mat denomUpdateVtStr = (Vt.rows(indexStr) * B.t() * B + beta / NspotsInStr / NspotsInStr * vecOnetnrMat * Vt.rows(indexStr) + lambda * VtDt);
                arma::mat updateVtStr = nomUpdateVtStr / denomUpdateVtStr;
                Vt.rows(indexStr) = Vt.rows(indexStr) % updateVtStr;
                arma::vec mean_VtStr = mean(Vt.rows(indexStr), 0).t();
                Mu_Str.col(istr) = mean_VtStr; 
            }
            VList(islice) = Vt;
            Mu_Slice = Mu_Slice + Mu_Str;
        }
        // update mu
        Mu = Mu_Slice / nSlice;
        // calculate objective function
        obj = calculateObject(XList,VList,AList,StrList,B,Mu,beta,lambda);
       return List::create(
                           Named("VList") = VList,
                           Named("Mu") = Mu,
                           Named("Obj") = obj);
        }//end try 
        catch (std::exception &ex)
        {
            forward_exception_to_r(ex);
        }
        catch (...)
        {
            ::Rf_error("C++ exception (unknown reason)...");
        }
        return R_NilValue;
} // end funcs

