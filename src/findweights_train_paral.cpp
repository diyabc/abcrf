#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

struct FindWeightsTrain : public Worker
{
  // trainingNodeID matrix
  const RMatrix<double> trainingNodeID;
  
  // inbag matrix
  const RMatrix<double> inbag;
  
  // trees and train
  const int ntrain;
  const int ntree;
  
  // destination matrix
  RMatrix<double> result;
  
  // initialize with source and destination
  FindWeightsTrain(const NumericMatrix trainingNodeID, const NumericMatrix inbag,
                   const int ntrain, const int ntree, NumericMatrix result)
    : trainingNodeID(trainingNodeID), inbag(inbag), ntrain(ntrain), ntree(ntree), result(result) {}
  
  // compute weight for one training
  void operator()(std::size_t begin, std::size_t end) {
    
    for (std::size_t trainIdx = begin; trainIdx < end; trainIdx++) {
      
      double Boob = 0;                               //the denominator for the weight computation
      
      for(int k=0; k<ntree; k++){
        
        double LbOob=0;
        
        std::vector<double> colK(ntrain);
        for(int i = 0; i < ntrain; i++) {
          colK[i] = inbag(i,k);  // store the inbag identifier
        }
        
        std::vector<double> nbOob(ntrain);
        if( inbag(trainIdx,k)==0 ){                // if the i-th training data is out-of-bag for tree k
          for(int i = 0; i < ntrain; i++) {
            nbOob[i] = 0;
            if(trainingNodeID(i,k) == trainingNodeID(trainIdx,k)) {   //identify what training data falls in the same leaf than the i-th training
              nbOob[i] = colK[i];
              LbOob = LbOob + nbOob[i];
            }
          }
          if(LbOob != 0){
            for(int i = 0; i < ntrain; i++) {
              nbOob[i] = nbOob[i] / LbOob;
            }
            for(int i = 0; i < ntrain; i++) {
              result(i,trainIdx) = result(i,trainIdx) + nbOob[i];                    //update weight vector
            }
          }
          Boob++;                                               //update the number of trees where oob
        }
      }
      
      if(Boob > 0){
        for(int i = 0; i < ntrain; i++) {
          result(i,trainIdx) = result(i,trainIdx) / Boob;
        }
      }
    }
    
  }
};

// [[Rcpp::export]]
NumericMatrix findweights_train_paral(NumericMatrix trainingNodeID,
                                      NumericMatrix inbag, 
                                      int ntrain, int ntree) {
  
  // allocate the output matrix
  NumericMatrix output(ntrain, ntrain);
  
  FindWeightsTrain findWeightsTrain(trainingNodeID, inbag, ntrain, ntree, output);
  
  // call parallelFor to do the work
  parallelFor(0, ntrain, findWeightsTrain);
  
  // return the output matrix
  return output;
}

/*** R
findweights_train_all_parallel <- function(trainingNodeID, inbag, ntrain, ntree,
                                           paral = FALSE, ncores = if(paral) max(detectCores()-1,1) else 1) {
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  result <- foreach(idxLoop=0:(ntrain-1), .combine='cbind') %dopar% {
    abcrf:::findweights_train(trainingNodeID = trainingNodeID, inbag = inbag, idxLoop,
                              ntrain = ntrain, ntree = ntree)
  }
  
  stopCluster(cl)
  
  return(result)
}

data(snp)
modindex <- snp$modindex
sumsta <- snp$sumsta[modindex == "3",]
r <- snp$param$r[modindex == "3"]
r <- r[1:5]
sumsta <- sumsta[1:5,]
data2 <- data.frame(r, sumsta)
model.rf.r <- regAbcrf(r~., data2, ntree=10)

ncores <- 3

obj <- model.rf.r$model.rf
nodeIDTrain <- predict(obj, data2, predict.all=TRUE, num.threads=ncores, type="terminalNodes")$predictions
nodeIDObs <- predict(obj, snp.obs, predict.all=TRUE, num.threads=ncores, type="terminalNodes")$predictions
inbag <- simplify2array(obj$inbag.counts)
ntrain <- obj$num.samples
ntree <- obj$num.trees

RcppParallel::setThreadOptions(numThreads = ncores)

weights_train_1 <- findweights_train_paral(nodeIDTrain, inbag, ntrain, ntree)

weights_train_2 <- findweights_train_all_parallel(nodeIDTrain, inbag, ntrain, ntree,
                                                  paral = TRUE, ncores = ncores)

all.equal(weights_train_1, weights_train_2)

microbenchmark::microbenchmark(findweights_train_paral(nodeIDTrain, inbag, ntrain, ntree),
                               findweights_train_all_parallel(nodeIDTrain, inbag, ntrain, ntree,
                                                              paral = TRUE, ncores = ncores),
                               times = 10)

*/

// void findweights_train_int(NumericMatrix trainingNodeID, NumericMatrix inbag, int ntrain, int trainIdx, int ntree,
//                            NumericMatrix result, int resIdx){
//   LogicalVector oobPerTree(ntrain);       //to store the number of out-of-bag per tree k
//   //LogicalVector nonOobPerTree(ntrain);
//   LogicalVector isEqual(ntrain);          //to identify the positions
//   NumericVector nbOob(ntrain);
//   double Boob = 0;                               //the denominator for the weight computation
//   NumericVector colK(ntrain);
//   double LbOob=NA_REAL;
//   for(int k=0; k<ntree; k++){
//     //oobPerTree = (inbag(_,k) == 0);       //identify what training data is out-of-bag in tree k
//     colK = inbag(_,k);                    // store the inbag identifier
//     //nonOobPerTree = (inbag(_,k) != 0);
//     if( inbag(trainIdx,k)==0 ){                // if the i-th training data is out-of-bag for tree k
//       isEqual = trainingNodeID(_,k) == trainingNodeID(trainIdx,k);   //identify what training data falls in the same leaf than the i-th training
//       nbOob = colK * as<NumericVector>(isEqual);              //recover a vector with the number of times an in-bag data fell with the i-th training
//       if(is_true(any(nbOob!=0))){                             //if there is at least one training inbag that falls that falls with the i-th training
//         LbOob = sum(nbOob);                                   //count their number (L_b.oob)
//         nbOob = nbOob/LbOob;
//         result(_,resIdx) = result(_,resIdx) + nbOob;                    //update weight vector
//         Boob++;                                               //update the number of trees where oob
//       }
//     }
//   }
//   
//   if(Boob > 0){
//     result(_,resIdx) = result(_,resIdx)/Boob;
//   }
// }

// NumericMatrix findweights_train(NumericMatrix trainingNodeID, NumericMatrix inbag, int ntrain, int trainIdx, int ntree){
//   NumericMatrix result(ntrain, 1);     //to store the result, a column vector with the oob weights for 1 training data
//   findweights_train_int(trainingNodeID,inbag, ntrain, trainIdx, ntree, result, 0);
//   return(result);
//   
// }

// NumericMatrix findweights_train_all(NumericMatrix trainingNodeID, NumericMatrix inbag, int ntrain, int ntree){
//   NumericMatrix result(ntrain,ntrain);
//   IntegerVector counti(ntrain);
//   double meancount;
//   for(int i=0; i<ntrain; i++){
//     for(int k=0; k<ntree; k++){
//       if( inbag(i,k) == 0 ){
//         meancount = 0;
//         for(int j=0; j<ntrain; j++){
//           if( ( trainingNodeID(j,k) == trainingNodeID(i,k) ) && (inbag(j,k) != 0) ){
//             counti[j] = inbag(j,k);
//             meancount = meancount + inbag(j,k);
//           } else{
//             counti[j] = 0;
//           }
//         }
//         if( meancount >=1 ){
//           for(int j=0; j<ntrain; j++){
//             result(j,i) = result(j,i) + counti[j]/ meancount;
//           }
//         }
//       }
//     }
//   }
//   return(result);
// }

// NumericMatrix findweights_train_all(NumericMatrix trainingNodeID, NumericMatrix inbag, int ntrain, int ntree){
//   NumericMatrix result(ntrain, ntrain);
//   for (int trainIdx = 0; trainIdx < ntrain; trainIdx++) {
//     result(_,trainIdx) = findweights_train(trainingNodeID, inbag, ntrain, trainIdx, ntree)(_,0);
//   }
//
//   return(result);
//
// }

// void operator()(std::size_t begin, std::size_t end) {
//   for (std::size_t trainIdx = begin; trainIdx < end; trainIdx++) {
//     
//     NumericVector nbOob(ntrain);
//     double Boob = 0;                               //the denominator for the weight computation
//     NumericVector colK(ntrain);
//     double LbOob=NA_REAL;
//     for(int k=0; k<ntree; k++){
//       for(int i = 0; i < ntrain; i++) {
//         colK[i] = inbag(i,k);  // store the inbag identifier
//       }                   
//       if( inbag(trainIdx,k)==0 ){                // if the i-th training data is out-of-bag for tree k
//         for(int i = 0; i < ntrain; i++) {
//           if(trainingNodeID(i,k) == trainingNodeID(trainIdx,k)) {   //identify what training data falls in the same leaf than the i-th training
//             nbOob[i] = colK[i];
//           } else {
//             nbOob[i] = 0;
//           }
//         }  
//         if(is_true(any(nbOob!=0))){                             //if there is at least one training inbag that falls that falls with the i-th training
//           LbOob = sum(nbOob);                                   //count their number (L_b.oob)
//           nbOob = nbOob/LbOob;
//           for(int i = 0; i < ntrain; i++) {
//             result(i,trainIdx) = result(i,trainIdx) + nbOob[i];                    //update weight vector
//           }  
//           Boob++;                                               //update the number of trees where oob
//         }
//       }
//     }
//     
//     if(Boob > 0){
//       for(int i = 0; i < ntrain; i++) {
//         result(i,trainIdx) = result(i,trainIdx) / Boob;
//       }  
//     }
//   }
//   
// }
