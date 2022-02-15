#' @title get_test_data
#' @description This function generates a bunch of testing data given the tree and true shifts
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param true_variable The true shifts
#' @param beta The coefficients of shifts
#' @param alpha The selective rate
#' @param n_test The number of testing datasets
#' @param xtype the type of design matrix
#' @return A matrix, each row is a testing dataset for the true model
#' @export
#' @import phylolm

get_test_data = function(tree, true_variable, beta, alpha, n_test,xtype=c('simpX','orgX')){
  sigma = rTrait(n=1,tree,model='OU',parameters = list(alpha=alpha,sigma2=2*alpha))
  if(xtype == 'orgX'){
      X = generate_design_matrix(tree,type='orgX',alpha = alpha)
  }else if(xtype == 'simpX'){
    X = generate_design_matrix(tree,type='simpX')
  }
  ret = if(length(true_variable)>1) as.data.frame(matrix(X[,true_variable] %*% beta + sigma,nrow=1)) else as.data.frame(matrix(X[,true_variable] * beta + sigma,nrow=1))
  for(i in 2:n_test){
      sigma = rTrait(n=1,tree,model='OU',parameters = list(alpha=alpha))
      Y = if(length(true_variable)>1) X[,true_variable] %*% beta + sigma else X[,true_variable] * beta + sigma
      ret[nrow(ret)+1,] = Y
  }
  return(as.matrix(ret))

}

#' @title get_prediction_likelihood
#' @description This function calculates the prediction log likelihood to evaluate the accuracy of estimates
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param selected_variable The variables selected
#' @param alpha The selective rate
#' @param test_data A bunch of testing data given by the true model
#' @param xtype the type of design matrix
#' @return The prediction log likelihood value
#' @export
#' @import phylolm

get_prediction_likelihood = function(tree,Y,selected_variable,alpha,test_data,xtype=c('simpX','orgX')){
    sigma2 = 2*alpha
    xtype = match.arg(xtype)
    if(xtype == 'orgX'){
      X = generate_design_matrix(tree,type='orgX',alpha = alpha)
  }else if(xtype == 'simpX'){
    X = generate_design_matrix(tree,type='simpX')
  }
    
    lmod = if(length(selected_variable)>0) phylolm(Y~X[,selected_variable],phy=tree,model='OUfixedRoot') else phylolm(Y~1,phy=tree,model='OUfixedRoot')

    Y_pred = lmod$fitted.values
    
    sqrtInvSigma = t(sqrt_OU_covariance(tree, alpha=alpha,root.model ="OUfixedRoot")$sqrtInvSigma)
    sqrtSigma = sqrt_OU_covariance(tree, alpha=alpha,root.model ="OUfixedRoot")$sqrtSigma

    cmp_likelihood = function(target,prediction,sqrtInvSigma,sqrtSigma,sigma2){
        n = length(target)
        return(-n*log(2*pi)/2-n/2*log(sigma2)-(t(prediction-target) %*% t(sqrtInvSigma) %*% sqrtInvSigma %*% (prediction-target))/(2*sigma2)-1/2*log(det(sqrtSigma %*% t(sqrtSigma))))
    }

    loglik = apply(test_data,1,cmp_likelihood,prediction=Y_pred,sqrtInvSigma = sqrtInvSigma,sqrtSigma = sqrtSigma,sigma2 = sigma2)

    return(mean(loglik))

}
