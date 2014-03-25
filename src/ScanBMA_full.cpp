// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::depends(RcppEigen)]] 

#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;

using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Function for getting the r^2 value for a linear regression
// y: the response vector
// x: the prediction matrix
// xIndices: the columns of x to use
//   automatically adds a column of ones
const double GetR2( NumericVector yy, NumericMatrix x, std::set<int> xIndices ) {
  NumericMatrix xsub( x.nrow(), 1 + (int)xIndices.size() );
  int i;
  for ( i = 0; i < x.nrow(); i++ ) {
    xsub(i,0) = 1;
  }
  i = 1;
  for ( std::set<int>::iterator it = xIndices.begin(); it != xIndices.end(); it++ ) {
    for ( int j = 0; j < x.nrow(); j++ ) {
      xsub(j,i) = x(j, (*it));
    }
    i++;
  }
  const Map<MatrixXd> X(as<Map<MatrixXd> >(xsub));
  const Map<VectorXd> y (as<Map<VectorXd> >(yy)) ;
  const int n(X.rows());
  const int p(X.cols());
  const LLT<MatrixXd> llt(MatrixXd(p,p).setZero().selfadjointView<Lower>().rankUpdate(X.adjoint()));
  const VectorXd betahat(llt.solve(X.adjoint() * y));
  const VectorXd fitted(X * betahat);
  const VectorXd resid(y - fitted);
  const double sse(resid.dot(resid));
  const double ybar(y.sum() / n);
  VectorXd ymybar(n);
  for ( int i = 0; i < n; i++ ) {
    ymybar[i] = y[i] - ybar;
  }
  const double sst(ymybar.dot(ymybar));

  return 1 - sse/sst;
};

class Model {
 public:
  std::set<int> modelIndices;
  double r2;
  double bic;
  Model(std::set<int>, double, double);

  const List Output() const {
    List ret;
    ret["modelLength"] = (int)modelIndices.size();
    ret["modelIndices"] = modelIndices;
    ret["r2"] = r2;
    ret["bic"] = bic;
    return ret;
  }

  bool operator<(const Model& b) const { return bic < b.bic; }
  bool operator>(const Model& b) const { return bic > b.bic; }
  bool operator>=(const Model& b) const { return !(*this < b); }
  bool operator<=(const Model& b) const { return !(*this > b); }
  bool operator==(const Model& b) const { return false; }
  bool operator!=(const Model& b) const { return !(*this == b); }
};


Model::Model( std::set<int> mInd, double mr2, double mbic ) {
  modelIndices = mInd;
  r2 = mr2;
  bic = mbic;
};

// Function for generating a model string from a set of indices
const std::string ModelString( std::set<int> indices ) {
  std::stringstream tmp;
  std::set<int>::iterator it;
  for ( it=indices.begin(); it != indices.end(); it++ ) {
    tmp << (*it) << ".";
  }
  return tmp.str();
};

// Function for returning an empty list (if no models were found)
const List EmptyReturn( int nModelsChecked ) {
  List ret;
  return ret;
};

// ScanBMA function using g-prior
// [[Rcpp::export]]
const List ScanBMA_BIC( NumericVector y,
                        NumericMatrix x,
                        NumericVector priorProbs_,
                        double oddsRatio ) {

  arma::vec priorProbs(priorProbs_);

  const double logOR = 2 * log(oddsRatio);
  const int n = y.size();
  const int H = x.ncol(); // # of predictors

  // diagnostics
  int curpass = 0;

  // Create a vector of indices sorted by the prior
  // then use that to get the correct indices
  // (decreasing)
  arma::uvec pord = sort_index(priorProbs, 1);

  // Compute log(prior/(1-prior))
  NumericVector logprior(H);
  for ( int i = 0; i < H; i++ ) {
    logprior[i] = log(priorProbs[pord[i]]/(1-priorProbs[pord[i]]));
  }

  // Active set (only the null model at first)
  // A sorted linked list
  // Attributes: set of model indices, r2, bic
  std::set<Model> keepModels;
  std::set<Model> activeModels;
  std::set<Model> nextModels;

  std::set<int> modelIndices;
  activeModels.insert(Model(modelIndices, 0, 0));

  // Checked set (empty at first)
  // A quickly searchable table
  // only want to add model strings and
  // check if a model string has been added
  // hashtable or map
  // boost library?
  std::tr1::unordered_set<std::string> checkedModels;
  checkedModels.insert("");

  // variables for use in the loops
  double minBIC = 0;
  double candidateR2 = 0;
  double candidateBIC = 0;
  double cutoffBIC = logOR;
  std::set<Model>::iterator it;
  std::set<int>::iterator intit;
  int h, i;
  std::string mString;

  // Loop through while we have active models
  // to search around
  while ( ((int)activeModels.size()) > 0 ) {
    curpass++;
    for ( i = 0; i < H; i++ ) {
      // Get the ordered predictor index
      h = pord[i];
      // Go through all active models and add or remove predictor h
      // then check the model if it has not been checked and add if
      // appropriate
      for ( it = activeModels.begin(); it != activeModels.end(); it++ ) {
        modelIndices = (*it).modelIndices;
        if (((int)modelIndices.count(h)) > 0) {
          modelIndices.erase(h);
        }
        else {
          modelIndices.insert(h);
        }
        mString = ModelString(modelIndices);
        // If we have not checked the model yet, check and add appropriately
        if ( ((int)checkedModels.count(mString)) < 1 ) {
          // get r2 value
          candidateR2 = GetR2( y, x, modelIndices );
	  candidateBIC = n * log(1-candidateR2) + ((int)modelIndices.size()) * log(n);
  //        printf( mString.c_str() );
  //        printf( ", candidateBIC = %f\\n", candidateBIC );
          for ( intit = modelIndices.begin(); intit != modelIndices.end(); intit++ ) {
            candidateBIC -= 2*logprior[(*intit)];
          }
          checkedModels.insert(mString);
          if ( candidateBIC - minBIC < logOR ) {
            // Add model
            nextModels.insert(Model( modelIndices, candidateR2, candidateBIC ));
            minBIC = ((minBIC < candidateBIC) ? minBIC : candidateBIC);
          }
        }
      }
    }
    cutoffBIC = minBIC + logOR;

    // Remove models with bic greater than the cutoff
    it = keepModels.begin();
    while( it != keepModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    keepModels.erase(it, keepModels.end());

    // Add good models from the active models to keep
    // and clear the active models
    it = activeModels.begin();
    while( it != activeModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    keepModels.insert(activeModels.begin(), it);
    activeModels.clear();

    // Add good models from the next models to the active models
    // and clear the next models
    it = nextModels.begin();
    while( it != nextModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    activeModels.insert(nextModels.begin(), it);
    nextModels.clear();
  }

  // For now, return raw results for processing in R (as a list?)
  int nModels = (int)keepModels.size();
  IntegerMatrix indicesMatrix(nModels, H);
  NumericVector r2s(nModels);
  NumericVector bics(nModels);
  IntegerVector sizes(nModels);

  std::set<int>::iterator itind;
  i = 0;
  for ( it = keepModels.begin(); it != keepModels.end(); it++ ) {
    for ( itind = (*it).modelIndices.begin();
          itind != (*it).modelIndices.end();
          itind++ ) {
      indicesMatrix(i, (*itind)) = 1;
    }
    r2s[i] = (*it).r2 * 100;
    bics[i] = (*it).bic;
    sizes[i] = (int)(*it).modelIndices.size();
    i++;
  }

  // Calculate posterior probabilities for the models
  // the maximum BIC is last (already sorted)
  const double maxBIC = bics[((int)bics.size()) - 1];
  NumericVector mbics = bics - maxBIC;
  NumericVector postprobs = exp(mbics * (-0.5));
  // Difference: R code includes check for is.na(postprob) and sets to 1
  double postprobnorm = 0;
  for ( i = 0; i < (int)postprobs.size(); i++ ) {
    postprobnorm += postprobs[i];
  }
  for ( i = 0; i < (int)postprobs.size(); i++ ) {
    postprobs[i] = postprobs[i]/ postprobnorm;
  }

  // Calculate posterior probabilities for the predictors
  NumericVector probne0(H);
  i = 0;
  for ( it = keepModels.begin(); it != keepModels.end(); it++ ) {
    for ( itind = (*it).modelIndices.begin();
          itind != (*it).modelIndices.end();
          itind++ ) {
      probne0(*itind) += postprobs[i] * 100;
    }
    i++;
  }

  List ret;
  ret["postprob"] = postprobs;
  ret["r2"] = r2s;
  ret["bic"] = bics;
  ret["size"] = sizes;
  ret["which"] = indicesMatrix;
  ret["probne0"] = probne0;
  ret["reduced"] = false;
  ret["n.models"] = nModels;
  ret["nmodelschecked"] = (int)checkedModels.size();

  keepModels.clear();
  checkedModels.clear();

  return ret;
}


// Differentiation function using g-prior
// [[Rcpp::export]]
const List BMA_Diff_BIC( NumericVector y,
			   NumericMatrix x,
			   NumericVector priorProbs_,
			   IntegerVector bestModel,
			   bool diff100,
			   bool diff0,
			   IntegerVector pred100,
			   IntegerVector pred0,
			   double minprob,
			   double epsilon ) {

arma::vec priorProbs(priorProbs_);

const int n = y.size();
const int H = x.ncol(); // # of predictors

int i;

// Create a vector of indices sorted by the prior
// then use that to get the correct indices
// (decreasing)
arma::uvec pord = sort_index(priorProbs, 1);
// Compute log(prior/(1-prior))
NumericVector logprior(H);
for ( int i = 0; i < H; i++ ) {
  logprior[i] = log(priorProbs[pord[i]]/(1-priorProbs[pord[i]]));
}

// Set up the best model
std::set<int> modelIndices;
for ( i = 0; i < bestModel.length(); i++ ) {
  modelIndices.insert(bestModel[i]);
}

// Compute the best model score (R^2)
double bestR2 = GetR2( y, x, modelIndices );
double otherR2, Ak;

NumericVector adjProb100(pred100.length());
NumericVector ak100(pred100.length());
NumericVector adjProb0(pred0.length());
NumericVector ak0(pred0.length());

 if ( diff100 ) {
   // for each predictor in pred100
   //   remove from best model
   //   compute score for submodel (Ak)
   //   compute adjusted probability (exp(Ak)/(1 + exp(Ak)))
   for ( i = 0; i < pred100.length(); i++ ) {
     modelIndices.erase(pred100[i]);
     otherR2 = GetR2( y, x, modelIndices );
     Ak = logprior[pred100[i]] - log(n)/2 - (n/2) * log((1-bestR2)/(1-otherR2));
     ak100[i] = Ak;
     adjProb100[i] = exp(Ak)/(1 + exp(Ak));
     modelIndices.insert(pred100[i]);
   }
 }

 if ( diff0 ) {
   double maxAdjProb0 = 0;

   // for each predictor in pred0
   //   add tobest model
   //   compute score for supermodel (Ak)
   //   compute adjusted probability (1/(1 + exp(Ak)))
   //   keep track of maximum adjusted probability
   for ( i = 0; i < pred0.length(); i++ ) {
     modelIndices.insert(pred0[i]);
     otherR2 = GetR2( y, x, modelIndices );
     Ak = -logprior[pred0[i]] + log(n)/2 - (n/2) * log((1-bestR2)/(1-otherR2));
     ak0[i] = Ak;
     adjProb0[i] = 1/(1 + exp(Ak));
     if ( adjProb0[i] > maxAdjProb0 ) {
       maxAdjProb0 = adjProb0[i];
     }
     modelIndices.erase(pred0[i]);
   }

   // if maxAdjProb0 > minprob
   //   adjProb0 = adjProb0 * minprob / (maxAdjProb0 + epsilon)
   if ( maxAdjProb0 > minprob ) {
     for ( i = 0; i < adjProb0.length(); i++ ) {
       adjProb0[i] = adjProb0[i] * minprob / maxAdjProb0 - epsilon;
     }
   }
 }

 // return list of adjusted probabilities
 List ret;
 ret["adjProb100"] = adjProb100;
 ret["adjProb0"] = adjProb0;
 ret["Ak100"] = ak100;
 ret["Ak0"] = ak0;
 return ret;
}

// ScanBMA function using g-prior
// [[Rcpp::export]]
const List ScanBMA_g( NumericVector y,
                        NumericMatrix x,
                        NumericVector priorProbs_,
                        double oddsRatio,
                        double g ) {

  arma::vec priorProbs(priorProbs_);

  const double logOR = 2 * log(oddsRatio);
  const int n = y.size();
  const int H = x.ncol(); // # of predictors

  // diagnostics
  int curpass = 0;

  // Create a vector of indices sorted by the prior
  // then use that to get the correct indices
  // (decreasing)
  arma::uvec pord = sort_index(priorProbs, 1);

  // Compute log(prior/(1-prior))
  NumericVector logprior(H);
  for ( int i = 0; i < H; i++ ) {
    logprior[i] = log(priorProbs[pord[i]]/(1-priorProbs[pord[i]]));
  }

  // Active set (only the null model at first)
  // A sorted linked list
  // Attributes: set of model indices, r2, bic
  std::set<Model> keepModels;
  std::set<Model> activeModels;
  std::set<Model> nextModels;

  std::set<int> modelIndices;
  activeModels.insert(Model(modelIndices, 0, 0));

  // Checked set (empty at first)
  // A quickly searchable table
  // only want to add model strings and
  // check if a model string has been added
  // hashtable or map
  // boost library?
  std::tr1::unordered_set<std::string> checkedModels;
  checkedModels.insert("");

  // variables for use in the loops
  double minBIC = 0;
  double candidateR2 = 0;
  double candidateBIC = 0;
  double cutoffBIC = logOR;
  std::set<Model>::iterator it;
  std::set<int>::iterator intit;
  int h, i;
  std::string mString;

  // Loop through while we have active models
  // to search around
  while ( ((int)activeModels.size()) > 0 ) {
    curpass++;
    for ( i = 0; i < H; i++ ) {
      // Get the ordered predictor index
      h = pord[i];
      // Go through all active models and add or remove predictor h
      // then check the model if it has not been checked and add if
      // appropriate
      for ( it = activeModels.begin(); it != activeModels.end(); it++ ) {
        modelIndices = (*it).modelIndices;
        if (((int)modelIndices.count(h)) > 0) {
          modelIndices.erase(h);
        }
        else {
          modelIndices.insert(h);
        }
        mString = ModelString(modelIndices);
        // If we have not checked the model yet, check and add appropriately
        if ( ((int)checkedModels.count(mString)) < 1 ) {
          // get r2 value
          candidateR2 = GetR2( y, x, modelIndices );
  //        candidateBIC = n * log(1-candidateR2) + ((int)modelIndices.size()) * log(n);
          candidateBIC = (n-1) * log(1 + g*(1-candidateR2)) + (1 + (int)modelIndices.size() - n) * log(1 + g);
  //        printf( mString.c_str() );
  //        printf( ", candidateBIC = %f\\n", candidateBIC );
          for ( intit = modelIndices.begin(); intit != modelIndices.end(); intit++ ) {
            candidateBIC -= 2*logprior[(*intit)];
          }
          checkedModels.insert(mString);
          if ( candidateBIC - minBIC < logOR ) {
            // Add model
            nextModels.insert(Model( modelIndices, candidateR2, candidateBIC ));
            minBIC = ((minBIC < candidateBIC) ? minBIC : candidateBIC);
          }
        }
      }
    }
    cutoffBIC = minBIC + logOR;

    // Remove models with bic greater than the cutoff
    it = keepModels.begin();
    while( it != keepModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    keepModels.erase(it, keepModels.end());

    // Add good models from the active models to keep
    // and clear the active models
    it = activeModels.begin();
    while( it != activeModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    keepModels.insert(activeModels.begin(), it);
    activeModels.clear();

    // Add good models from the next models to the active models
    // and clear the next models
    it = nextModels.begin();
    while( it != nextModels.end() && (*it).bic <= cutoffBIC ) {
      it++;
    }
    activeModels.insert(nextModels.begin(), it);
    nextModels.clear();
  }

  // For now, return raw results for processing in R (as a list?)
  int nModels = (int)keepModels.size();
  IntegerMatrix indicesMatrix(nModels, H);
  NumericVector r2s(nModels);
  NumericVector bics(nModels);
  IntegerVector sizes(nModels);

  std::set<int>::iterator itind;
  i = 0;
  for ( it = keepModels.begin(); it != keepModels.end(); it++ ) {
    for ( itind = (*it).modelIndices.begin();
          itind != (*it).modelIndices.end();
          itind++ ) {
      indicesMatrix(i, (*itind)) = 1;
    }
    r2s[i] = (*it).r2 * 100;
    bics[i] = (*it).bic;
    sizes[i] = (int)(*it).modelIndices.size();
    i++;
  }

  // Calculate posterior probabilities for the models
  // the maximum BIC is last (already sorted)
  const double maxBIC = bics[((int)bics.size()) - 1];
  NumericVector mbics = bics - maxBIC;
  NumericVector postprobs = exp(mbics * (-0.5));
  // Difference: R code includes check for is.na(postprob) and sets to 1
  double postprobnorm = 0;
  for ( i = 0; i < (int)postprobs.size(); i++ ) {
    postprobnorm += postprobs[i];
  }
  for ( i = 0; i < (int)postprobs.size(); i++ ) {
    postprobs[i] = postprobs[i]/ postprobnorm;
  }

  // Calculate posterior probabilities for the predictors
  NumericVector probne0(H);
  i = 0;
  for ( it = keepModels.begin(); it != keepModels.end(); it++ ) {
    for ( itind = (*it).modelIndices.begin();
          itind != (*it).modelIndices.end();
          itind++ ) {
      probne0(*itind) += postprobs[i] * 100;
    }
    i++;
  }

  List ret;
  ret["postprob"] = postprobs;
  ret["r2"] = r2s;
  ret["bic"] = bics;
  ret["size"] = sizes;
  ret["which"] = indicesMatrix;
  ret["probne0"] = probne0;
  ret["reduced"] = false;
  ret["n.models"] = nModels;
  ret["nmodelschecked"] = (int)checkedModels.size();

  keepModels.clear();
  checkedModels.clear();

  return ret;
}


// Differentiation function using g-prior
// [[Rcpp::export]]
const List BMA_Diff_g( NumericVector y,
			   NumericMatrix x,
			   NumericVector priorProbs_,
			   double g,
			   IntegerVector bestModel,
			   bool diff100,
			   bool diff0,
			   IntegerVector pred100,
			   IntegerVector pred0,
			   double minprob,
			   double epsilon ) {

arma::vec priorProbs(priorProbs_);

const int n = y.size();
const int H = x.ncol(); // # of predictors

int i;

// Create a vector of indices sorted by the prior
// then use that to get the correct indices
// (decreasing)
arma::uvec pord = sort_index(priorProbs, 1);
// Compute log(prior/(1-prior))
NumericVector logprior(H);
for ( int i = 0; i < H; i++ ) {
  logprior[i] = log(priorProbs[pord[i]]/(1-priorProbs[pord[i]]));
}

// Set up the best model
std::set<int> modelIndices;
for ( i = 0; i < bestModel.length(); i++ ) {
  modelIndices.insert(bestModel[i]);
}

// Compute the best model score (R^2)
double bestR2 = GetR2( y, x, modelIndices );
double otherR2, Ak;

NumericVector adjProb100(pred100.length());
NumericVector ak100(pred100.length());
NumericVector adjProb0(pred0.length());
NumericVector ak0(pred0.length());

 if ( diff100 ) {
   // for each predictor in pred100
   //   remove from best model
   //   compute score for submodel (Ak)
   //   compute adjusted probability (exp(Ak)/(1 + exp(Ak)))
   for ( i = 0; i < pred100.length(); i++ ) {
     modelIndices.erase(pred100[i]);
     otherR2 = GetR2( y, x, modelIndices );
     Ak = logprior[pred100[i]] - log(1+g)/2 - ((n-1)/2) * log((1 + g*(1-bestR2))/(1 + g*(1-otherR2)));
     ak100[i] = Ak;
     adjProb100[i] = exp(Ak)/(1 + exp(Ak));
     modelIndices.insert(pred100[i]);
   }
 }

 if ( diff0 ) {
   double maxAdjProb0 = 0;

   // for each predictor in pred0
   //   add tobest model
   //   compute score for supermodel (Ak)
   //   compute adjusted probability (1/(1 + exp(Ak)))
   //   keep track of maximum adjusted probability
   for ( i = 0; i < pred0.length(); i++ ) {
     modelIndices.insert(pred0[i]);
     otherR2 = GetR2( y, x, modelIndices );
     Ak = -logprior[pred0[i]] + log(1+g)/2 - ((n-1)/2) * log((1 + g*(1-bestR2))/(1 + g*(1-otherR2)));
     ak0[i] = Ak;
     adjProb0[i] = 1/(1 + exp(Ak));
     if ( adjProb0[i] > maxAdjProb0 ) {
       maxAdjProb0 = adjProb0[i];
     }
     modelIndices.erase(pred0[i]);
   }

   // if maxAdjProb0 > minprob
   //   adjProb0 = adjProb0 * minprob / (maxAdjProb0 + epsilon)
   if ( maxAdjProb0 > minprob ) {
     for ( i = 0; i < adjProb0.length(); i++ ) {
       adjProb0[i] = adjProb0[i] * minprob / maxAdjProb0 - epsilon;
     }
   }
 }

 // return list of adjusted probabilities
 List ret;
 ret["adjProb100"] = adjProb100;
 ret["adjProb0"] = adjProb0;
 ret["Ak100"] = ak100;
 ret["Ak0"] = ak0;
 return ret;
}
