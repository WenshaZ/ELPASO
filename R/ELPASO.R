#' @title forward_selection
#' @description This function provides forward selection given a variable ranking sequence
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param seq A vector gives that selection order of the shifts. (For example, the ensembled rank)
#' @param increasing If a smaller value in param seq means more important variable, increasing = TRUE, vice versa.
#' @param criterion the model selection criterion
#' @param stop TRUE/FALSE. If the selection process will be stopped if score stops decreasing
#' @param maxShifts the maximum number of shifts
#' @param opt options for l1ou functions
#' @return 
#' \item{s.c}{The position of shifts}
#' \item{score}{The criterion score}
#' @export

forward_selection = function(tree,Y,seq,increasing = TRUE,criterion = c('BIC','pBIC','pBICess','AICc','mBIC'),stop = TRUE,maxShifts = 20,opt){
  if(increasing == TRUE){
    v_rank = as.numeric(substr(names(sort(seq,FALSE)),2,5))
  } else{
    v_rank = as.numeric(substr(names(sort(seq,TRUE)),2,5))
  }
  if(stop==TRUE){
    s.c = NULL
    min.score = get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),NULL,opt)$score
    for(i in 1:maxShifts){
      model =  tryCatch(get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),v_rank[1:i],opt) ,error = function(e) NA)
      if(is.na(model)){
        break
      }
      else{
              score = model$score
            if(score <= min.score){
              min.score = score
              s.c = v_rank[1:i]
            }
            else{
              break
            }
          }
      }

  }
  
  if(stop==FALSE){
    s.c = NULL
    min.score = get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),NULL,opt)$score
    for(i in 1:maxShifts){
     model =  tryCatch(get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),c(s.c,v_rank[i]),opt) ,error = function(e) NA)
     if(is.na(model)){
      next
     }
     else{
      score = model$score
      if(score<=min.score){
        min.score = score
        s.c = c(s.c,v_rank[i])
      }
     }
  }
}
  
  s.c = if (is.null(s.c)) numeric(0) else s.c
  return(list('score'=min.score,'s.c'=s.c))
  
  
} 

#' @title backward_selection
#' @description This function provides forward selection given a variable ranking sequence
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param seq A vector gives that selection order of the shifts. (For example, the ensembled rank)
#' @param increasing If a smaller value in param seq means more important variable, increasing = TRUE, vice versa.
#' @param criterion the model selection criterion
#' @param stop TRUE/FALSE. If the selection process will be stopped if score stops decreasing
#' @param maxShifts the maximum number of shifts
#' @param opt options for l1ou functions
#' @return 
#' \item{s.c}{The position of shifts}
#' \item{score}{The criterion score}
#' @export

backward_selection = function(tree,Y,seq,increasing = TRUE,criterion = c('BIC','pBIC','pBICess','AICc','mBIC'),stop=TRUE,maxShifts = 10,opt){
  if(increasing == TRUE){
    v_rank = as.numeric(substr(names(sort(seq,FALSE)),2,5))[1:min(maxShifts,length(seq))]
  } else{
    v_rank = as.numeric(substr(names(sort(seq,TRUE)),2,5))[1:min(maxShifts,length(seq))]
  }
  if(stop==TRUE){
      s.c = v_rank
      m = maxShifts
      min.score = Inf
      for(i in maxShifts:0){
         the.sc = if (i>0) v_rank[1:i] else  NULL
          model =  tryCatch(get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),the.sc,opt),error = function(e) NA)
          if(is.na(model)){
            next
          } else{
            score = model$score
            if(score<=min.score){
              min.score = score
              s.c = the.sc
            }
            else{
              break
            }
          }

      }
        

  }
  if(stop==FALSE){
    s.c = v_rank
    min.score = tryCatch(get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),s.c,opt)$score,error = function(e) Inf)
      for(i in maxShifts:1){
          the.sc = setdiff(s.c,v_rank[i])
          model = tryCatch(get(paste('cmp_',criterion,sep = ''))(tree,as.matrix(Y),the.sc,opt),error = function(e) NA)
          if(is.na(model)){
            
            s.c = the.sc
            next

          } else{
            score = model$score
            if(score<=min.score){
              min.score = score
              s.c = the.sc
            }
          }
      }
  }
  
  
  s.c = if (is.null(s.c)) numeric(0) else s.c
  return(list('score'=min.score,'s.c'=s.c))
  
}

#' @title forward_backward_selection
#' @description This function provides forward+bacward selection given a variable ranking sequence
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param seq A vector gives that selection order of the shifts. (For example, the ensembled rank)
#' @param increasing If a smaller value in param seq means more important variable, increasing = TRUE, vice versa.
#' @param criterions A two-element vector. The model selection criterions for forward selection and backward selection respectively
#' @param stops A two-element vector. Each element is TRUE/FALSE. If the selection process will be stopped if score stops decreasing 
#' @param maxShifts the maximum number of shifts
#' @param direction 'same'(default)/'opposite'. If the forward selection and the backward selection are produced in the same direction.
#' @param opt options for l1ou functions
#' @return 
#' \item{s.c}{The position of shifts}
#' \item{score}{The criterion score}
#' @export


forward_backward_selection = function(tree,Y,seq,increasing = TRUE,criterions,stops,maxShifts,direction = c('same','opposite'),opt){
  f = forward_selection(tree,Y,seq,increasing = increasing, 
  criterion = criterions[1],stop = stops[1], maxShifts = maxShifts,opt=opt)
  f.sc = f$s.c
  if(length(f.sc)==0){
    return(f)
  }
  seq2 = 1:length(f.sc)
  names(seq2) = paste('X',f.sc,sep="")
  if(direction =='same'){
    bresult = backward_selection(tree,Y,seq2,increasing = FALSE, criterion = criterions[2],stop= stops[2],maxShifts=min(length(seq2),maxShifts),opt = opt)
    return(list('score'=bresult$score,'s.c'=bresult$s.c[length(bresult$s.c):1]))    
  } else if(direction == 'opposite'){
    bresult = backward_selection(tree,Y,seq2,increasing = TRUE, criterion = criterions[2],stop= stops[2],maxShifts=min(length(seq2),maxShifts),opt = opt)
    
    s.c = if (is.null(s.c)) numeric(0) else s.c
    return(list('score'=bresult$score,'s.c'=bresult$s.c))
  }
}


#' @title ELPASO 
#' @description  This function is the main function of the package, the inputs of the function are the trait vector and phylogenetic tree, 
#' the outputs are the estimated positions and change values of the shifts. And the results are based on ensemble variable selection.
#' @param tree ultrametric tree of class phylo with branch lengths, and edges in postorder
#' @param Y Y trait vector without missing entries
#' @param criterion the criterion for model selection, "BIC"(default),"pBIC"
#' @param maxShifts The max number of shifts. Integer.
#' @param nsamples The number subsamples are generated. Integer.
#' @param xtype type of design matrix, 'simpX'(simple design matrix, default) or 'orgX'(original design matrix)
#' @param penalty Penalty function to be applied. Either "LASSO"(default) or "SCAD"
#' @param ensemble_method, the way to ensemble the ranking seqs. 'median','quantile'(default),'arith.mean','geom.mean'
#' @param q, if using "quantile" as the ensemble method, specify the quantile value in this parameter. 
#' The default value is 0.25, which is a good choice given by the simulation results
#' @return 
#' \item{Y}{input trait vector/matrix.}
#' \item{tree}{input tree.}
#' \item{ensemble_rank}{The final ensemble rank of shift positions}
#' \item{shifts}{estimated shift positions, i.e. vector of indices of edges where the estimated shifts occur.}
#' \item{shift.means}{estimates change of the expectation of the shift values}
#' \item{nShifts}{estimated number of shifts.}
#' \item{alpha}{maximum likelihood estimate of the adaptation rate \eqn{\alpha}{alpha}}
#' \item{sigma2}{maximum likelihood estimate of the variance rate \eqn{\sigma^2}{sigma^2}}
#' \item{fitted.values}{fitted values, i.e. estimated trait means.}
#' \item{residuals}{residuals.}
#' \item{logLik}{log likelihood of given model}
#' \item{criterion}{The criterion for model selection}
#' \item{score}{information criterion value of the estimated shift configuration.}
#' \item{penalty}{Penalty function to be applied}
#' @export 
#' @import phylolm
#' @import l1ou
#' @examples
#' require(l1ou)
#' data('lizard.tree')
#' data('lizard.traits')
#' tree = lizard.tree
#' Y = as.vector(lizard.traits[,1])
#' ELMODEL = ELPASO(tree,Y)


ELPASO = function(tree, Y, criterion = c('BIC','pBIC'), maxShifts = 20,nsamples = 200, xtype=c('simpX','orgX'), penalty = c('LASSO','SCAD'),ensemble_method = "quantile", q = 0.25 ){
  criterion =  match.arg(criterion)
  xtype = match.arg(xtype)
  penalty = match.arg(penalty)


  sc_list = c()
  score_list = c()
  alpha_list = c()
  sigam2_list = c()

  bmod = phylolm(Y~1,phy=tree,model='BM')
  rank_seqs = get_ranking_seqs(tree,Y,nsamples = nsamples,size = 0.8*length(Y),replace = FALSE,alpha = 0,sigma2 = bmod$sigma2,sigma2_error = 0,xtype = xtype,penalty = penalty)

  ensemble_rank = combine_ranking_seqs(rank_seqs,method = ensemble_method,q=q)
  rank_all = as.matrix(ensemble_rank)

  opt = list()
  opt$multivariate.missing = FALSE
  opt$Z = l1ou:::generate_design_matrix(tree,type='simpX')
  opt$alpha.upper.bound = l1ou:::alpha_upper_bound(tree)
  opt$quietly = FALSE
  opt$root.model = "OUfixedRoot"
  opt$alpha.lower.bound = NA
  opt$alpha.starting.value = NA


  if(criterion == "pBIC"){
    selection_res = backward_selection(tree,Y,ensemble_rank,TRUE,'pBIC',FALSE,maxShifts,opt)
  }else{
    selection_res = forward_backward_selection(tree,Y,ensemble_rank,TRUE,c(criterion, criterion),c(FALSE,FALSE),maxShifts=50,direction='same',opt)

  }

  sv = selection_res$s.c

  sc_list = c(sc_list, if(length(sv)>0) paste(sort(sv),collapse = ',') else '')
  score_list = c(score_list, selection_res$score)

  lmod = if(length(sv)>0) phylolm(Y~X[,sv],phy=tree,model='OUfixedRoot') else phylolm(Y~1,phy=tree,model='OUfixedRoot')
  alpha_list = c(alpha_list, lmod$optpar)
  sigam2_list = c(sigam2_list, lmod$sigma2)

  for(i in 2:3){
    
    rank_seqs = get_ranking_seqs(tree,Y,nsamples = nsamples,size = 0.8*length(Y),replace = FALSE,alpha = lmod$optpar,sigma2 = lmod$sigma2,sigma2_error = 0,xtype = xtype,penalty = penalty)
    ensemble_rank = combine_ranking_seqs(rank_seqs,method = ensemble_method,q=q)
    rank_all = cbind(rank_all, as.matrix(ensemble_rank))
    if(criterion == "pBIC"){
      selection_res = backward_selection(tree,Y,ensemble_rank,TRUE,'pBIC',FALSE,maxShifts,opt)
    }else{
      selection_res = forward_backward_selection(tree,Y,ensemble_rank,TRUE,c(criterion, criterion),c(FALSE,FALSE),maxShifts=50,direction='same',opt)

    }
    sv = selection_res$s.c

    sc_list = c(sc_list, if(length(sv)>0) paste(sort(sv),collapse = ',') else '')
    score_list = c(score_list, selection_res$score)

    lmod = if(length(sv)>0) phylolm(Y~X[,sv],phy=tree,model='OUfixedRoot') else phylolm(Y~1,phy=tree,model='OUfixedRoot')
    alpha_list = c(alpha_list, lmod$optpar)
    sigam2_list = c(sigam2_list, lmod$sigma2)
  }

  ssv = as.numeric(strsplit(sc_list[which.min(score_list)],split=',')[[1]])
  alpha = alpha_list[which.min(score_list)]
  sigma2 = sigam2_list[which.min(sigam2_list)]
  lmod = if(length(ssv)>0) phylolm(Y~X[,ssv],phy=tree,model='OUfixedRoot') else phylolm(Y~1,phy=tree,model='OUfixedRoot')
  shift.means = as.matrix(lmod$coefficients[-1])
  final_ensemble_rank = rank_all[,which.min(score_list)]
  rownames(shift.means)=NULL

return(list(
  'Y' = Y,
  'tree' = tree,
  'shifts' = ssv,
  'shift.means' = shift.means,
  'nShifts' = length(ssv),
  'alpha' = alpha,
  'sigma2' = sigma2,
  'fitted.values' = lmod$fitted.values,
  'residuals' = lmod$residuals,
  'logLik' = lmod$logLik,
  'criterion' = criterion,
  'score' = min(score_list),
  'penalty' = penalty,
  'ensemble_rank' = final_ensemble_rank


  ))


}

cmp_BIC = l1ou:::cmp_BIC
cmp_pBIC = l1ou:::cmp_pBIC
cmp_pBICess = l1ou:::cmp_pBICess
cmp_AICc = l1ou:::cmp_AICc
cmp_mBIC = l1ou:::cmp_mBIC