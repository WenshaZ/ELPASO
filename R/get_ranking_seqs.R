#' @title get_ranking_seqs 
#' @description  This function generates subsamples of a given combination of tree and trait values, applys LASSO/SCAD on each subsample, and produces ranking sequence of each subsample
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param nsamples The number subsamples are generated. Integer.
#' @param size The size of each subsample. Integer.
#' @param replace if replace while subsampling, TRUE/FALSE(default)
#' @param alpha the adaption rate (known or estimated) 
#' @param sigma2 the variance of OU model (known or estimated)
#' @param sigma2_error the variance of gaussian measurement error
#' @param xtype type of design matrix, 'simpX'(simple design matrix, default) or 'orgX'(original design matrix)
#' @param penalty Penalty function to be applied. Either "LASSO"(default) or "SCAD"
#' @return A matrix, each row is a ranking sequence of variables, each column represents the rank sets of each variable
#' @export 
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg setupLambda
#' @import l1ou
#' @importFrom PIGShift OU.vcv

get_ranking_seqs = function(tree,Y,nsamples,size,replace=FALSE,alpha,sigma2,sigma2_error=0,xtype=c('simpX','orgX'), penalty = c('LASSO','SCAD')){
  #xtype = match.arg(xtype)
  #penalty = match.arg(penalty)

  if(alpha==0|xtype=='simpX'){
    X = generate_design_matrix(tree,type='simpX')
  }else{
    X = generate_design_matrix(tree,type='orgX',alpha)
  }

  seqs = X[0,]

  
    if(sigma2_error == 0){
      t = t(sqrt_OU_covariance(tree, alpha=alpha,root.model = "OUfixedRoot")$sqrtInvSigma)

    }else{
      Vc = OU.vcv(tree,alpha)*sigma2 + diag(sigma2_error,length(tree$tip.label))
      cp = svd(Vc)
      t = cp$u %*% diag(1/sqrt(cp$d)) %*% t(cp$v)

    }
    X_t = t %*% X
    Y_t = t %*% Y

    for(i in 1:nsamples){
      index = sample(1:nrow(X),size=size,replace=replace)
      X_sample = X_t[index,]
      Y_sample = Y_t[index,]

      lambda = setupLambda(X_sample, Y_sample, family='gaussian', alpha=1, lambda.min=0.001, nlambda=100, penalty.factor=rep(1,ncol(X)))
      lasso = glmnet(as.matrix(X_sample), as.vector(Y_sample),
                     family='gaussian',intercept = FALSE,standaraize=FALSE)

      scad = ncvreg(X_sample,Y_sample,
                    family='gaussian',penalty='SCAD',normalize=FALSE,lambda = lambda)

      if(penalty=='LASSO'){
        the_rank = rank(-rowSums(as.matrix(lasso$beta[1:(nrow(lasso$beta)),]!=0)))
      }else{
        the_rank = rank(-rowSums(scad$beta[2:(nrow(scad$beta)),]!=0))
      }
      
      seqs = rbind(seqs,the_rank)
    }
  

  return(seqs)
}

#' @title combine_ranking_seqs
#' @description ensemble the ranking sequnces with median/arith.mean/geom.mean/quantile
#' @param rank_seqs A ranking sequence of all the variables
#' @param method The ensemble method "median"/"arith.mean"/"geom.mean"/"quantile"
#' @param q Required only when method is "quantile", specify the quantile value
#' @return A vector (length is number of potential shift postions). The ensembled ranking sequence
#' @export
#' @importFrom psych geometric.mean

combine_ranking_seqs = function(rank_seqs,method = "quantile",q=0.25){
 if(method == "min"){
   cseq = apply(rank_seqs,2,min)
 } else if(method == "median"){
   cseq = apply(rank_seqs,2,median)
 } else if(method == "arith.mean"){
   cseq = apply(rank_seqs,2,mean)
 } else if(method == "geom.mean"){
   cseq = apply(rank_seqs,2,geometric.mean)
 } else if(method == "quantile"){
   cseq = apply(rank_seqs,2,function(x)quantile(x,probs=q))
 }
 return(cseq)
}

generate_design_matrix = l1ou:::generate_design_matrix
