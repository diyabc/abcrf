#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix findweights_oob(NumericMatrix trainingNodeID, NumericMatrix testingNodeID, IntegerMatrix inbag, int ntrain, int nnew, int ntree){
  NumericMatrix result(ntrain,nnew);
  IntegerVector counti(ntrain);
  double sumweights;
  for(int i=0; i<nnew; i++){
    sumweights = 0;
    for(int j=0; j<ntrain; j++){
      for(int k=0; k<ntree; k++){
        if( ( trainingNodeID(j,k) == testingNodeID(i,k) ) && (inbag(j,k) == 0) ){
          result(j,i) = result(j,i) + 1;
        }
      }
      sumweights = sumweights + result(j,i);
    }
    for(int j=0; j<ntrain; j++){
      result(j,i) = result(j,i) / sumweights;
    }
  }
  return(result);
}
