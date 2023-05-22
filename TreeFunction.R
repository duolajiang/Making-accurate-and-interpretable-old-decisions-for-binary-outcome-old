# xvars  <- c('age_at_diagnosis','lymph_assessed','max_LoS','BMI','pT_num',
#               'Stage','sex','MSI','BRAF','RAS','ASAclass','pN_num','cM_num',
#               "colon_perforation",'lym_agioinvasie',"Grade",'Morphology','sub_loc','PerformanceStatus'
#   )
# dat <- data[,xvars]
# dat$trt <- data$combined_chemo
# dat$y <- ate.test.mean
# # COLUMNS OF COVARIATES
# split.var <- 1:19;ctg=6:19; mtry=length(split.var) 
# #split.var <- 1:16;ctg=7:16; mtry=length(split.var) 
# test=NULL; 
# # FLEXIBILTIY OF TREE STRUCTURE 
# min.ndsz=min.ndsz # mininum number of observations in a node
# pt=pt # minimun number of treatment/control in a node 
# max.depth=max.depth # max depth of tree
# alpha=alpha # only when significance level larger than 0.05, then the split is taken 
#   
# tree <- grow.INT(dat, test=NULL, 
#                  split.var,ctg,
#                  min.ndsz, pt=pt, max.depth,alpha,
#                  mtry=length(split.var))
#   
# varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
# cate <- NULL
# for (var in varname){
#   if(class(dat[,var])=='factor'){
#     cate <- cbind(cate,var)
#   }
# }



# =============================================================
# for printing
# =============================================================
GetNodeName <- function(tree,dat,cate,nodeid){
  #browser()
  isleft <- (substr(nodeid,start = nchar(nodeid),stop = nchar(nodeid))==1)
  node.parent.id <- substr(nodeid,start = 1,stop = nchar(nodeid)-1)
  node.parent <- tree[tree$node==node.parent.id,]
  var.name <- as.character(node.parent[,'vname'])
  if(is.element(var.name,cate)){
    cut <- as.character(unlist(strsplit(as.character(node.parent[,'cut']),split=" ")))
    cut.index <- is.element(levels(dat[,var.name]),cut)
    cut.val.left <- levels(dat[,var.name])[cut.index];   cut.val.left <- CombineLevelsString(cut.val.left)
    cut.val.right <- levels(dat[,var.name])[!cut.index]; cut.val.right <- CombineLevelsString(cut.val.right)
    target <- ifelse(isleft==TRUE,paste(var.name,'=',cut.val.left,sep = ''),paste(var.name,'=',cut.val.right,sep = ''))
  } else{
    cut.val <- as.character(node.parent[,'cut'])
    target <- ifelse(isleft==TRUE,paste(var.name,'<=',cut.val,sep = ''),paste(var.name,'>',cut.val,sep = ''))
  }
  return(target)
}
CombineLevelsString <- function(stringarray){
  single <- NULL
  for(i in seq(length(stringarray))){
    single <- paste(single,as.character(stringarray[i]),'')
  }
  return(single)
}


# =================================================
# THE grow.INT() FUNCTION CONSTRUCTS A LARGE TREE 
# =================================================
grow.INT <- function(dat, test=NULL, 
                     split.var,ctg=NULL,
                     min.ndsz=20, pt=pt, max.depth=15,alpha,
                     mtry=length(split.var))
{
  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(dat); 
  if (!is.null(test)) list.test <- list(test)
  name <- 1
  #i <- 1
  while (length(list.nd)!=0) {    
    for (i in 1:length(list.nd)){
      print(i)
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        split <- partition.INT(list.nd[[i]], test0, name[i], min.ndsz=min.ndsz, 
                               pt=pt, split.var=split.var, ctg=ctg, max.depth=max.depth, alpha=alpha, mtry=mtry)
        out <- rbind(out, split$info)
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          
        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }   
  out$node <- as.character(out$node)
  out <- out[order(out$node),] 
  return(out)
}

grow.INT.ILL <- function(dat, 
                         split.var,ctg=NULL,
                         min.ndsz, pt, max.depth,
                         mtry,
                         alpha,loc)
{
  out <- list.nd <- temp.list <- temp.name <- NULL
  list.nd <- list(dat); 
  name <- 1
  #i <- 1
  while (length(list.nd)!=0) {    
    for (i in 1:length(list.nd)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        split <- partition.INT.ILL(list.nd[[i]], name[i], min.ndsz=min.ndsz, 
                                   pt=pt, split.var=split.var, ctg=ctg, max.depth=max.depth,mtry=mtry,alpha,loc)
        out <- rbind(out, split$info)
        if (!is.null(split$left)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } 
      }
    }
    list.nd <- temp.list; name <- temp.name
    temp.list <-  temp.name <- NULL
  }   
  out$node <- as.character(out$node)
  out <- out[order(out$node),]
  return(out)
}

# =================================================
# THE PARTITION INTINIAL LARGE TREE 
# =================================================
partition.INT <- function(dat, test=NULL, name="1", min.ndsz=20, pt=pt, 
                          split.var, ctg=NULL, 
                          max.depth=15, alpha = 0.05, mtry=length(split.var))
{   
  # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES. 
  #browser()
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat); 
  if (!is.null(test)) {n.test <- NA; score.test <- NA;}  ########## TEST SAMPLE ########  
  var <- vname <- NA; cut <- NA; max.score <- threhold.score <- qt(1-alpha/2,df=n-2);   
  trt <- dat$trt; y <- dat$y; vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  #browser()
  trt.effect <- NA; n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
  if (min(n.1, n.0) >= n0) {trt.effect <- mean(y)}
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name) 
  if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)  
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));  
      if(length(temp) > 1) { 
        if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
        else zcut <- temp[-length(temp)]
        # print(i); print(temp); print(zcut)
        for(j in zcut) {
          score <- NA
          if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
          else  {grp <- sign(x <= j); cut1 <- as.character(j)}   
          #browser()
          score <- ttest(dat, z=grp, min.ndsz, pt)
          # print(cbind(var=i, cut=j, score=score))
          if (!is.na(score) && score > max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j} 
        }}}}
  if (!is.null(test)) { 
    n.test <- nrow(test); score.test <- NA;
    if (!is.na(var)) {
      if (is.element(var,ctg)) grp.test <- sign(is.element(test[,var], best.cut))                              ############################
      else  grp.test <- sign(test[,var] <= best.cut)   
      score.test <- ttest(test, z=grp.test, n0=(n0/2))
      if (!is.na(score.test)){
        out$name.l <- name.l; out$name.r <- name.r
        out$left.test <- test[grp.test==1,  ]
        out$right.test <- test[grp.test==0,  ]
        if (is.element(var,ctg)) {                                                                               ############################                      
          out$left  <- dat[is.element(dat[,var], best.cut),]      
          out$right <- dat[!is.element(dat[,var], best.cut), ]}
        else {
          out$left  <- dat[dat[,var]<= best.cut,]
          out$right <- dat[dat[,var]> best.cut, ]
        }
      }
      else {var <- NA; vname <- NA; cut <- NA;  max.score <- NA}
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                             var = var, vname=vname, cut= cut, score=ifelse(max.score==threhold.score, NA, max.score), 
                             score.test, size.test=n.test)
    }
    else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,
                             var = NA, vname=NA, cut= NA, score=NA, 
                             score.test=NA, size.test=n.test)
    }
  }	else {
    if (!is.na(var)) {
      out$name.l <- name.l; out$name.r <- name.r
      if (is.element(var,ctg)) {                                                                               ############################                      
        out$left  <- dat[is.element(dat[,var], best.cut),]      
        out$right <- dat[!is.element(dat[,var], best.cut), ]}
      else {
        out$left  <- dat[dat[,var]<= best.cut,]
        out$right <- dat[dat[,var]> best.cut, ]
      }
      out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=dim(out$left)[1],n.r=dim(out$right)[1],
                             trt.effect=trt.effect, lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]], var = var, vname=vname, cut= cut, 
                             score=ifelse(max.score==threhold.score, NA, max.score))
    }
    else{
      out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect,
                             lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                             var=NA, vname=NA, cut=NA, score=NA)
    }
  }
  out 
}



# =================================================
# THE PARTITION INTINIAL LARGE TREE: use increment of log-loglikelihood as the score of splits
# =================================================
partition.INT.ILL <- function(dat, name="1", min.ndsz, pt, 
                              split.var, ctg=NULL, 
                              max.depth, mtry=length(split.var),alpha,loc)
{   
  # NOTE THAT CTG INDICATES THE COLUMNS FOR CATEGORICAL VARIABLES. 
  #browser()
  call <- match.call(); out <- match.call(expand = F)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  name.l <- paste(name, 1, sep=""); name.r <- paste(name, 2, sep="")
  n <- nrow(dat); 
  var <- vname <- NA; cut <- NA; max.score <- 0;   
  trt <- dat$trt; y <- dat$y; vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  #browser()
  trt.effect <- ifelse(loc,sum(dat$weight*y)/sum(dat$weight),mean(y))
  n.1 <- sum(trt==1); n.0 <- n-n.1 ; n0 <- n*pt
  # CONTROL THE MAX TREE DEPTH
  depth <- nchar(name) 
  if (depth <= max.depth && n >= min.ndsz && min(n.1, n.0) >= n0) {
    print(c(n,min.ndsz))
    m.try <- ifelse(is.null(mtry), length(split.var), mtry)  
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(x));  
      if(length(temp) > 1) { 
        if (is.element(i,ctg)) zcut <- power.set(temp)                                                  ############################ CLASS VARIABLE
        else zcut <- temp[-length(temp)]
        # print(i); print(temp); print(zcut)
        for(j in zcut) {
          score <- NA
          if (is.element(i,ctg)) {grp <- sign(is.element(x, j)); cut1 <- paste(j, collapse=" ")}      ############################ CLASS VARIABLE
          else  {grp <- sign(x <= j); cut1 <- as.character(j)}   
          #browser()
          if (loc) { score <- WeightedIncrementLL(dat, grp, min.ndsz ,pt) - alpha }
          else { score <- IncrementLL(dat, grp, min.ndsz ,pt) - alpha }
          # print(cbind(var=i, cut=j, score=score))
          if (!is.na(score) && score > max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j} 
        }}}}
  if (!is.na(var)) {
    out$name.l <- name.l; out$name.r <- name.r
    if (is.element(var,ctg)) {                                                                               ############################                      
      out$left  <- dat[is.element(dat[,var], best.cut),]      
      out$right <- dat[!is.element(dat[,var], best.cut), ]}
    else {
      out$left  <- dat[dat[,var]<= best.cut,]
      out$right <- dat[dat[,var]> best.cut, ]
    }
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=dim(out$left)[1],n.r=dim(out$right)[1],
                           trt.effect=trt.effect, lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]], var = var, vname=vname, cut= cut, 
                           score=max.score)
  }
  else{
    out$info <- data.frame(node=name, size = n, n.t=n.1, n.c=n.0, n.l=NA, n.r=NA, trt.effect=trt.effect,
                           lower= quantile(y,0.025)[[1]], upper=quantile(y,0.975)[[1]],
                           var=NA, vname=NA, cut=NA, score=NA)
  }
  return(out)
}



# ==============================================
# PREDICT TREATMENT OR CONTROL
# tree: tree structure
# test: data for prediction
# ctg: categorical covariate index
# ==============================================
predict_initree_ITE <- function(tree,test,ctg){
  #browser()
  y <- NULL
  for (i in seq(dim(test)[1])){
    obs <- test[i,]
    node.id <- 1
    node <- tree[tree$node==node.id,]
    while (!is.na(node$vname)){
      varname <- as.character(node$vname)
      if(is.element(node$var,ctg)){
        cut.value <- as.character(node$cut)
        if(is.element(obs[,varname],strsplit(cut.value,' '))){
          node.id <- paste(node.id,'1',sep = '')
        } else{
          node.id <- paste(node.id,'2',sep = '')
        }
      } else{
        # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
        cut.value <- as.numeric(levels(node$cut)[node$cut]) 
        if(obs[,varname]<=cut.value){
          node.id <- paste(node.id,'1',sep = '')
        }else{
          node.id <- paste(node.id,'2',sep = '')
        }
      }
      node <- tree[tree$node==node.id,]
    }
    y <- c(y,node$trt.effect)
  }
  return(y)
}


# ==============================================
# PREDICT node ID
# tree: tree structure
# test: data for prediction
# ctg: categorical covariate index
# ==============================================
predict_initree_NodeID <- function(tree,test,ctg){
  y <- c()
  for (i in seq(dim(test)[1])){
    obs <- test[i,]
    node.id <- 1
    node <- tree[tree$node==node.id,]
    while (!is.na(node$vname)){
      varname <- as.character(node$vname)
      if(is.element(node$var,ctg)){
        cut.value <- as.character(node$cut)
        if(is.element(obs[,varname],strsplit(cut.value,' '))){
          node.id <- paste(node.id,'1',sep = '')
        } else{
          node.id <- paste(node.id,'2',sep = '')
        }
      } else{
        # cant use as.numeric(node$cut); below is efficient thant as.numeric(as.character(node$cut))
        cut.value <- as.numeric(levels(node$cut)[node$cut]) 
        if(obs[,varname]<=cut.value){
          node.id <- paste(node.id,'1',sep = '')
        }else{
          node.id <- paste(node.id,'2',sep = '')
        }
      }
      node <- tree[tree$node==node.id,]
    }
    y <- c(y,node[,1])
  }
  return(y)
}

# ==============================================
# FUNCTION rdat() SIMULATES A SIMPLE DATA SET
# ==============================================

rdat <- function(n=100, K =50, 
                 beta0=2, beta1=2, beta2=2, beta3=2, beta4=2, beta5=2, 
                 sigma=1, cut1=.5, cut2=.5)
{
  trt <- sample(0:1, n, replace=T)
  #### Generate Covariates
  for (j in 1:4) assign(paste("x", j, sep=""),  sample(1:K, n, replace=T)/K)    
  ### 
  mean <- beta0 + beta1*trt + beta2*sign(x1<=cut1) + beta3*sign(x2<=cut2) + 
    beta4*sign(x1<=cut1)*trt + beta5*sign(x2<=cut2)*trt
  y <- mean + rnorm(n, mean=0, sd=sigma)            
  ##### Output
  data.frame(x1=x1, x2=x2, x3=x3, x4=x4, y=y, trt=trt)
}



# ===============================================================================
# FUNCTION ttest() COMPUTES THE t TEST FOR INTERACTION WITH CONTINUOUS RESPONSE
# ===============================================================================
IncrementLL <- function(dat, z, min.ndsz, pt)
{  
  trt <- dat$trt
  n <- nrow(dat)
  n.l <- sum(z==1)
  n.r <- sum(z==0)
  # n.l.t <- sum((z==1)&(trt==1))/sum((z==1))
  # n.l.c <- sum((z==1)&(trt==0))/sum((z==1))
  # n.r.t <- sum((z==0)&(trt==1))/sum((z==0))
  # n.r.c <- sum((z==0)&(trt==0))/sum((z==0))
  #if ((n.l<min.ndsz)|(n.r<min.ndsz)|(min(n.l.t,n.l.c,n.r.t,n.r.c) < pt)) {
  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  } 
  
  y.l <- dat[z==1,'y']
  y.r <- dat[z==0,'y']
  y <- dat$y
  
  mu.l <- mean(y.l)
  mu.r <- mean(y.r)
  mu   <- mean(y)
  
  predvar.l <- dat[z==1,'predictive_var']
  predvar.r <- dat[z==0,'predictive_var']
  predvar   <- dat$predictive_var
  
  sigma.l <- (sum(predvar.l) + sum((y.l-mu.l)^2))/n.l
  sigma.r <- (sum(predvar.r) + sum((y.r-mu.r)^2))/n.r
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n
  
  ll.l <- -1/2*n.l*log(sigma.l)
  ll.r <- -1/2*n.r*log(sigma.r)
  ll   <- -1/2*n*log(sigma)
  
  score <- ll.l+ll.r-ll
  
  return(score)
}



IncrementLLNoProject <- function(dat, z, min.ndsz, pt)
{  
  trt <- dat$trt
  n <- nrow(dat)
  n.l <- sum(z==1)
  n.r <- sum(z==0)

  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  } 
  
  y.l <- dat[z==1,'y']
  y.r <- dat[z==0,'y']
  y <- dat$y
  
  mu.l <- mean(y.l)
  mu.r <- mean(y.r)
  mu   <- mean(y)
  
  sigma.l <- (sum((y.l-mu.l)^2))/n.l
  sigma.r <- (sum((y.r-mu.r)^2))/n.r
  sigma   <- (sum((y-mu)^2))/n
  
  ll.l <- -1/2*n.l*log(sigma.l)
  ll.r <- -1/2*n.r*log(sigma.r)
  ll   <- -1/2*n*log(sigma)
  
  score <- ll.l+ll.r-ll
  
  return(score)
}



WeightedIncrementLL <- function(dat, z, min.ndsz, pt)
{  
  trt <- dat$trt
  # n <- nrow(dat)
  # n.l <- sum(z==1)
  # n.r <- sum(z==0)
  
  n <- sum(dat$weight)
  n.l <- sum(dat[z==1,'weight'])
  n.r <- sum(dat[z==0,'weight'])
  
  if ((n.l<min.ndsz)|(n.r<min.ndsz)) {
    score <- NA
    return(score)
  } 
  
  y.l <- dat[z==1,'y']*dat[z==1,'weight']
  y.r <- dat[z==0,'y']*dat[z==0,'weight']
  y <- dat$y*dat$weight
  
  mu.l <- sum(y.l)/n.l
  mu.r <- sum(y.r)/n.r
  mu   <- sum(y)/n
  
  predvar.l <- dat[z==1,'predictive_var']*(dat[z==1,'weight'])^2
  predvar.r <- dat[z==0,'predictive_var']*(dat[z==0,'weight'])^2
  predvar   <- dat$predictive_var*(dat$weight)^2
  
  sigma.l <- (sum(predvar.l) + sum((y.l-mu.l)^2))/n.l
  sigma.r <- (sum(predvar.r) + sum((y.r-mu.r)^2))/n.r
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n
  
  ll.l <- -1/2*n.l*log(sigma.l)
  ll.r <- -1/2*n.r*log(sigma.r)
  ll   <- -1/2*n*log(sigma)
  
  score <- ll.l+ll.r-ll
  
  return(score)
}


# =================================
# METHOD I: THE TEST SAMPLE METHOD
# =================================

# -----------------------------------------------------------------------------
# THE PRUNING ALGORITHM GOOD FOR THE TREE SIZE SELECTION VIA TEST SAMPLE METHOD
# -----------------------------------------------------------------------------
#tre <- tree
prune.size.testsample <- function(tre)
{
  #browser()
  out <- as.list(NULL)     
  ntest <- as.numeric(tre[1, ncol(tre)])
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  result <- NULL; n.tmnl <- sum(is.na(tre$var)); subtree <- 1            
  a <- cbind(Ga.2=20, Ga.3=30, Ga.4=40, Ga.BIC=log(ntest))
  max.Ga <- rep(-1e20, 4); size <- rep(0, 4); btree <-as.list(1:4) 
  while (n.tmnl > 1 ) {
    # print(tre)
    internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
    r.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- as.numeric(as.vector(branch$score))
      r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
    }
    alpha <- min(r.value)
    nod.rm <- internal[r.value == alpha]; 
    if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    G <- sum(as.numeric(as.vector(tre$score)), na.rm=T); 
    Ga <- G - a*l 
    for (k in 1:4){if (Ga[k] > max.Ga[k]) {max.Ga[k] <- Ga[k]; size[k] <- n.tmnl; btree[[k]] <- tre}}                        
    result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                  size.tmnl=nrow(tre)-l, alpha=alpha, G=G, Ga))
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score")] <- NA
    n.tmnl <- sum(is.na(tre$cut))
    if (n.tmnl ==1) {for (k in 1:4){if (0 > max.Ga[k]) {max.Ga[k] <- 0; size[k] <- 1; btree[[k]] <- tre}}}
    subtree <- subtree + 1          
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                size.tmnl=1, alpha=9999, G=0, Ga=cbind(Ga.2=0, Ga.3=0, Ga.4=0, Ga.BIC=0)))     
  result <- as.data.frame(result)
  out$result <- result; out$size <- size; out$btree <- btree
  out 
}



# ==========================================================
# FUNCTION de() FINDS ALL THE DESCENDENTS OF NODE x IN tree
# ==========================================================
de <- function(x, tree)
{
  #browser()
  if(length(x) != 1) stop("The length of x in function de must be 1.")    
  y <- tree$node;  de <- NA
  if(sum(match(x, y), na.rm = T) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1    
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1  #find the subtree's the least right leaf node, then the next node is another tree, therefore finding the first node whose index larger the root node and depth is the same as the root node, then one node before this node is the end of the subtree (the last node of the subtree)
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }}
  de
}



# ===========================================================================
# THE power.set() FUNCTION PROVIDES THE POWER SET FOR A CATEGORICAL VARIABLE
# ===========================================================================
power.set <- function(x) {
  if(length(x) == 0) return(vector(mode(x), 0))
  x <- sort(unique(x)); n <- length(x); K <- NULL
  for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out); i <- 1
  out[!sapply(out, length)>=ceiling(n/2+.5)]
}


# THIS FUNCTION as.numeric.factor() CONVERTS FACTOR INTO NUMERIC
as.numeric.factor <- function(x){as.numeric(levels(x))[x]}





# ========================================================================
# FUNCTION obtain.btree() OBTAINS THE BEST SUBTREE WITH KNOW SIZE bsize=
# ========================================================================
obtain.btree <- function(tre, bsize=6)
{
  btre <- NULL
  if (bsize==1) { btre <- tre[1,]; btre[, c("var", "cut", "score", "score.test")] <- NA}  
  else if (bsize <1) stop("THE BEST TREE SIZE bsize= MUST BE >=1!")
  else {
    n.tmnl <- sum(is.na(tre$cut)); indicator <- T  
    if (bsize > n.tmnl) stop("THE BEST TREE SIZE bsize PROVIDED IS LARGER THAN THE FULL TREE THAT YOU HAVE GROWN.")          
    while (n.tmnl >= bsize && indicator ==T) {
      # print(tre); print(cbind(n.tmnl, bsize))
      internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
      r.value <- 1:l
      for(i in 1:l) {
        branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
        score <- as.numeric(as.vector(branch$score))
        r.value[i] <- sum(score, na.rm=T) / sum(!is.na(score))
      }
      alpha <- min(r.value)
      nod.rm <- internal[r.value == alpha]; 
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score", "score.test")] <- NA
      n.tmnl <- sum(is.na(tre$cut))      
      # print(cbind(n.tmnl, bsize))  
      if (n.tmnl==bsize) {btre <- tre; print(btre); indicator <- F}  
    }
  }
  if (is.null(btre)) print(paste("The optimally-pruned subtree sequence does not have a subtree of bsize = ", bsize, sep="")) 
  return(btre)
}




# ================================================================
# TREE AUTOMATICALLY GENERATION; DATA IS CONVERTED INTO PARTY CLASS
# i: node (like 1,11,112)
# data: tree
# no return
# ================================================================
tree_search <- function(i,data){
  #browser()
  if(!is.element(i,data$node)){
    print(i)
    if (substr(i,nchar(i),nchar(i))==1){
      ii <- paste('split=sp_',substr(i,1,nchar(i)-1),sep='')
      ss <<- gsub(ii,'',ss,fixed=TRUE)
      return(NULL)
    } else{
      return(NULL)
    }
  }
  print(i)
  n <<- n+1
  left.kid <- paste(i,1,sep = '')
  right.kid <- paste(i,2,sep = '')
  
  if (substr(i,nchar(i),nchar(i))==1){
    if((!is.element(left.kid,data$node))&(!is.element(right.kid,data$node))){
      ss <<- paste(ss,'partynode(',n,'L,',paste('split=sp_',i,sep = ''),'info=',round(data[data$node==i,'trt.effect'],3),',kids=list{',sep = '')
    }else{
      ss <<- paste(ss,'partynode(',n,'L,',paste('split=sp_',i,sep = ''),',kids=list{',sep = '')
    }
  } else{
    if((!is.element(left.kid,data$node))&(!is.element(right.kid,data$node))){
      ss <<- paste(ss,',partynode(',n,'L,',paste('split=sp_',i,sep = ''),'info=',round(data[data$node==i,'trt.effect'],3),',kids=list{',sep = '')
    }else{
      ss <<- paste(ss,',partynode(',n,'L,',paste('split=sp_',i,sep = ''),',kids=list{',sep = '')
    }
  }
  
  tree_search(paste(i,1,sep = ''),data)
  tree_search(paste(i,2,sep = ''),data)
  ss <<- paste(ss,'})',sep = '')
}


# ================================================================
# GENERATE SP_NODE <- PARTYSPLIT(VAR,BREAKS(INDEX))
# dat: data for growing tree
# tree: tree
# ================================================================
generate_split <- function(dat,tree){
  #browser()
  cate <- NULL
  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }
  
  tree_split <- tree[!is.na(tree$var),]
  split_party <- NULL
  for (i in 1:dim(tree_split)[1]){
    name <- as.character(tree_split[i,'vname'])
    if (is.element(name,cate)){
      n <- length(levels(dat[,name]))
      index.cut <- 'c('
      cut.val <- unlist(strsplit(as.character(tree_split[i,'cut']),split=" "))
      for (j in 1:n){
        if (is.element(levels(dat[,name])[j],cut.val)) {index.cut <- ifelse(j==n,paste(index.cut,'1L)',sep = ''),paste(index.cut,'1L,',sep = ''))}
        else {index.cut <- ifelse(j==n,paste(index.cut,'2L)',sep = ''),paste(index.cut,'2L,',sep = ''))}
      }
      split_party <- paste(split_party,'sp_',tree_split[i,'node'],'<-','partysplit(',tree_split[i,'var'],'L,','index=',index.cut,');',sep = '')
    } else{
      split_party <- paste(split_party,'sp_',tree_split[i,'node'],'<-','partysplit(',tree_split[i,'var'],'L,','breaks=',tree_split[i,'cut'],');',sep = '')
    }
  }
  return(split_party)
}


# =============================================================
# PRUNE THE TREE BASED ON THE COST
# i: node, like 1, 11, 112
# data: tree 
# =============================================================
prune_search <- function(i,data){
  #browser()
  if(!is.element(i,data$node)){
    print(i)
    return(NULL)
  }
  print(i)
  
  left.subtree <- prune_search(paste(i,1,sep = ''),data)
  right.subtree <- prune_search(paste(i,2,sep = ''),data)
  tree <- rbind(data[data$node==i,],left.subtree,right.subtree)
  
  if((!is.null(left.subtree))&(!is.null(right.subtree))){
    n.l <- dim(left.subtree)[1]; n.r <- dim(right.subtree)[1]
    if(n.l==n.r){
      n.t <- sum(left.subtree$trt == right.subtree$trt)
      if (n.t==n.l){
        tree <- data[data$node==i,]
        tree[,c('n.l','n.r','var','vname','cut','score')] <- NA
        return(tree)
      } else{
        return(tree)
      }
    } else{
      return(tree)
    }
  } else{
    return(tree)
  } 
}

# ====================================
# EXTRACT AND PRINT RULES 
# ====================================
print_leaf_nodes <- function(i,path,data,dat,cate){
  #browser()
  if(!is.element(i,data$node)){
    return(NULL)
  }
  
  if (i!=1){
    parent.node <- substr(i,1,nchar(i)-1)
    vname <- as.character(data[data$node==parent.node,'vname'])
    if(substr(i,nchar(i),nchar(i))==1){
      if(is.element(vname,cate)){
        cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
        cut.index <- is.element(levels(dat[,vname]),cut)
        cut.val <- levels(dat[,vname])[cut.index]
        cut <- NULL
        for(j in 1:length(cut.val)){
          cut <- paste(cut,cut.val[j],sep = ' ')
        }
        vname.con <- paste(vname,'in', cut, sep = ' ')
        if (parent.node==1){
          path <- paste(path,'{',vname.con,sep = ' ')
        } else{
          path <- paste(path,vname.con,sep = ' AND ')
        }
      } else{
        cut <- data[data$node==parent.node,'cut']
        vname.con <- paste(vname,'<=', cut, sep = ' ')
        path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
      }
    } else{
      if(is.element(vname,cate)){
        cut <- as.character(unlist(strsplit(as.character(data[data$node==parent.node,'cut']),split=" ")))
        cut.index <- !is.element(levels(dat[,vname]),cut)
        cut.val <- levels(dat[,vname])[cut.index]
        cut <- NULL
        for(j in 1:length(cut.val)){
          cut <- paste(cut,cut.val[j],sep = ' ')
        }
        vname.con <- paste(vname,'in', cut, sep = ' ')
        if (parent.node==1){
          path <- paste(path,'{',vname.con,sep = ' ')
        } else{
          path <- paste(path,vname.con,sep = ' AND ')
        }
      } else{
        cut <- data[data$node==parent.node,'cut']
        vname.con <- paste(vname,'>', cut, sep = ' ')
        path <- ifelse(parent.node==1,paste(path,'{',vname.con,sep = ' '),paste(path,vname.con,sep = ' AND '))
      }
    }
  }
  
  if((!is.element(paste(i,1,sep = ''),data$node))&(!is.element(paste(i,2,sep = ''),data$node))){
    
    path <- paste(path,'}','ate=',round(data[data$node==i,'trt.effect'],3),'(',round(data[data$node==i,'lower'],3),',',round(data[data$node==i,'upper'],3),')','size=',data[data$node==i,'size'],'t%=', round(data[data$node==i,'n.t']/data[data$node==i,'size'],3)*100, sep = ' ')
    print(path)
  }
  
  print_leaf_nodes(paste(i,1,sep = ''),path,data,dat,cate)
  print_leaf_nodes(paste(i,2,sep = ''),path,data,dat,cate)
  
}

# ====================================
# compute the number of leaf node vs. sum of increment of loglikelihood
# tree is the initial large tree whose alpha=0
# alpha is a given pernalty term, 
# we aim to plot two curves: x=number of leaves, y = 1) loglikeloood of the tree; 2) tree score (LL-alpha*#leaeves) 
# ====================================
LLAndTreescoreGivenAlpha <- function(dat,tree,alpha){
  list.nd <- temp.list <- NULL
  list.nd <- tree[1,]
  
  n <- nrow(dat)
  y <- dat$y; mu<- mean(y)
  predvar   <- dat$predictive_var
  sigma   <- (sum(predvar) + sum((y-mu)^2))/n
  num.leaf <- 1
  LLTree  <- -1/2*n*log(sigma)
  Scoretree <- LLTree - num.leaf*alpha
  leaf_LL_Score <- c(num.leaf,LLTree,Scoretree)
  
  while (!is.null(list.nd)){    
    for (i in 1:dim(list.nd)[1]){
      node <- list.nd[i,]
      if (is.element(paste(node[,1],'1',sep = ''),tree$node) && is.element(paste(node[,1],'2',sep = ''), tree$node)) {
        left <- as.numeric(paste(node[,1],'1',sep = '')); right <- as.numeric(paste(node[,1],'2',sep = ''))
        temp.list <- rbind(temp.list, tree[tree$node==left,], tree[tree$node==right,])
        num.leaf <- (num.leaf-1) + 2
        LLTree <- LLTree + node[,'score']
        Scoretree <- LLTree - num.leaf*alpha
        leaf_LL_Score <- rbind(leaf_LL_Score,c(num.leaf,LLTree,Scoretree))
      }
    }
    list.nd <- temp.list; 
    temp.list <- NULL
  }
  leaf_LL_Score <- as.data.frame(leaf_LL_Score)
  colnames(leaf_LL_Score) <- c('treesize','LL','score')
  return(leaf_LL_Score)
} 

# score = c(rank(delta LL(namely, score in tree) in each internal node, increasing))
# for (i in score){
#   continue split, step i, check object function C(sum of leaf nodes sigma - alpha*leaf_nodes), if Ci < Ci-1, then stop growing,
#   compute sum of log-likelihood of all the data of test sample
#   repeat 5 times,
#}
# compute alpha from bottom to top 
# ======================================================
#tre <- tree
Get_Alpha_TreeSize <- function(tre) {
  n.tmnl <- sum(is.na(tre$cut)); 
  Alpha <- NULL
  TreeSize <- NULL
  while (n.tmnl > 1) {
    internal <- tre$node[!is.na(tre$cut)]; l <- length(internal); 
    g.value <- 1:l
    for(i in 1:l) {
      branch <- tre[is.element(tre$node,c(internal[i], de(internal[i], tree=tre))),]
      score <- branch$score
      g.value[i] <- sum(score, na.rm=T) / (sum(is.na(score))-1)
    }
    alpha <- min(g.value)
    Alpha <- c(Alpha,alpha)
    nod.rm <- internal[g.value == alpha]; 
    tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
    tre[match(nod.rm, tre$node), c("var", "vname", "cut", "score")] <- NA
    n.tmnl <- sum(is.na(tre$cut))   
    treesize <- sum(is.na(tre$score))
    TreeSize <- c(TreeSize,treesize)
    #print(treesize)
    #print(tre)
  }
  return(list(Alpha,TreeSize))
}


GetTrainSample <- function(samdat,n_sam){
  a <- sample(samdat,n_sam)
  return(a)
}


test.Rsquare <- function(test, tree,ctg){
  Outcome <- predict_initree_ITE(tree,test,ctg)
  SSres <- sum((test$y-Outcome)^2+test$predictive_var)
  SStot <- sum((dat$y-mean(dat$y))^2+dat$predictive_var)
  Rsquare <- 1 - SSres/SStot
  return(Rsquare)
}

test.logsigma <- function(test, tree,ctg){
  Outcome <- predict_initree_ITE(tree,test,ctg)
  Risk.test <- sum(log((test$y-Outcome)^2+test$predictive_var))
  return(Risk.test)
}




PlotAte <- function(ate){
  ate <- as.data.frame(ate)
  colnames(ate) <- c("ITE")
  p <- ggplot(data=ate,aes(ITE))+
    geom_density(color="lightblue", fill="lightblue") +
    ggtitle("Distribution of x_{i} individual treatment effect")
  return(p)
}


# input: 
# alpha: pernalty term of mode complexity(number of internal splits) 
# alpha=0: to grow a large inital tree without prunng
# alpha >0: to grow and prune a tree given a fixed alpha 
# dat(processed dat),xvars(for tree grow),min.ndsz(minimun number of obs in node),
# pt(minimun overlap between treatment and control groups), max.depth(max depth of tree)
TreeGrownAndPrune <- function(dat,xvars,min.ndsz,pt,max.depth,alpha,loc){
  # COLUMNS OF COVARIATES
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var) 
  # FLEXIBILTIY OF TREE STRUCTURE 
  min.ndsz=min.ndsz # mininum number of observations in a node
  pt=pt # minimun number of treatment/control in a node 
  max.depth=max.depth # max depth of tree
  
  tree <- grow.INT.ILL(dat, 
                       split.var,ctg,
                       min.ndsz, pt, max.depth,
                       mtry=length(split.var),alpha,loc)
  return(tree)
}

# input: processed data, tree
GetCategoryVarInTree <- function(dat,tree){
  varname <- as.character(unique(tree$vname[!is.na(tree$vname)]))
  cate <- NULL
  for (var in varname){
    if(class(dat[,var])=='factor'){
      cate <- cbind(cate,var)
    }
  }
  return(cate)
}

# input: processed data,xvars for growing tree, initial large tree
# output: alpha and error of pruned tree
TreeComplexity5foldCV <- function(dat,xvars,tree,min.ndsz,max.depth,pt,loc){
  ##### generate test index for k-fold cross validation
  k <- 5
  n_test <- floor(nrow(dat)/k)
  Test <- list()
  Test[[1]] <- GetTrainSample(rownames(dat),n_test)
  Test[[2]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),Test[[1]]),]),n_test)
  Test[[3]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]])),]),n_test)
  Test[[4]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]],Test[[3]])),]),n_test)
  Test[[5]] <- GetTrainSample(rownames(dat[!is.element(rownames(dat),c(Test[[1]],Test[[2]],Test[[3]],Test[[4]])),]),n_test)
  
  Alpha_TreeSize <- Get_Alpha_TreeSize(tree)
  Alpha <- Alpha_TreeSize[[1]]
  TreeSize <- Alpha_TreeSize[[2]]
  Risk <- NULL
  split.var <- 1:length(xvars);ctg=6:length(xvars); mtry=length(split.var) 
  
  CV <- function(i){
    test <- dat[is.element(rownames(dat),Test[[i]]),]
    train <- dat[!is.element(rownames(dat),Test[[i]]),]
    R.fold <- NULL
    
    for (j in 1:length(Alpha)){
      tree <- grow.INT.ILL(train, 
                           split.var,ctg,
                           min.ndsz=min.ndsz, pt=pt, max.depth,
                           mtry,
                           alpha = Alpha[j],loc)
      R.test <- test.logsigma(test,tree,ctg)
      R.fold <- c(R.fold,R.test)
    }
    return(R.fold)
  }
  
  Risk <- mclapply(seq(k),CV,mc.cores = 3)
  Risk <- t(as.data.frame(Risk))
  Risk <- apply(Risk, 2, sum)/dim(dat)[1]
  Risk <- rbind(Alpha,TreeSize,Risk)
  
  return(Risk)
}

plotLeafScore <- function(df){
  df <- as.data.frame(df)
  p1 <- ggplot(data=df, aes(x=treesize, y=score, color=Alpha)) +
    geom_line(linetype = "solid")+
    geom_point() +
    ggtitle(expression("Tree size versus G(T;"~alpha~")"))+
    xlab("Tree Size")+
    ylab(expression("G(T;"~alpha~")"))+
    #color=expression(alpha),size=6)+
    theme(legend.position = c(0.8, 0.3),
          legend.text = element_text(size=10),
          legend.background = element_rect(fill=alpha('white', 0.5)),
          plot.title = element_text(size=10)) +
    scale_colour_discrete(expression(alpha))
  return(p1)
}

TreeSize_Risk <- function(df,vadjust){
  df <- df[,order(df[2,])]
  df <- t(df)
  df <- as.data.frame(df)
  df[,1] <- round(df[,1],digits = 1)
  n <- dim(df)[1]
  alpha.yposition <- max(df[,'Risk'])
  if((floor(n)/2)<8){
    df[c(1,seq(3,floor(n)/2,1),seq(floor(n)/2+2,n)),1] <- NA 
  } else if((floor(n)/3)<8){
    df[c(1,3,4,5,6,seq(8,floor(n)/2,1),seq(floor(n)/2+2,n-1)),1] <- NA 
  } else{
    df[c(1,3,4,5,6,seq(8,floor(n)/3,1),seq(floor(n)/3+2,2*floor(n)/3,1),seq(2*floor(n)/3+2,n-1,1)),1] <- NA 
  }
  p1 <- ggplot(data=df, aes(x=TreeSize, y=Risk,label=Alpha)) +
    geom_point()+
    scale_x_continuous(breaks=seq(min(df[,2]),max(df[,2]),2)) +
    geom_line(linetype = "solid")+
    coord_cartesian(clip = 'off') +
    labs(#title = expression("Cross validated estimate of risk"),
      subtitle= expression(alpha),
      x="Tree Size",y="Risk")+
    theme(plot.margin = unit(c(1.5,1,1,1), "lines"),
          plot.subtitle = element_text(hjust = 0.5,vjust = 4.5,size=10),
          plot.title = element_text(vjust = 5.5,size=10)) +
    geom_text(y=max(alpha.yposition)+vadjust,size=3)
  return(p1)
}











PlotCostSensTree <- function(tree.pruned,dat,cate,cost){
  tree.pruned$trt <- 1*(tree.pruned$trt.effect > cost)
  tree.pruned <- prune_search(1,tree.pruned)
  path <- NULL
  print_leaf_nodes(1,path,tree.pruned,dat,cate)
  SankeyNetworkPlot(tree.pruned,dat,cate,link_group = TRUE)
}




Treesize_Score <- function(R.ave,dat,tree,alpha.costum){
  alpha.cv <- R.ave[1,which(R.ave[3,]==max(R.ave[3,]))]; alpha.cv <- round(alpha.cv,digits = 1)
  alpha.aic <- 2
  alpha.bic <- log(nrow(dat))
  leaf_LL_Score.BIC <- LLAndTreescoreGivenAlpha(dat,tree,alpha.bic)
  leaf_LL_Score.AIC <- LLAndTreescoreGivenAlpha(dat,tree,alpha.aic)
  leaf_LL_Score.CV  <-  LLAndTreescoreGivenAlpha(dat,tree,alpha.cv)
  leaf_LL_Score.CUS <-  LLAndTreescoreGivenAlpha(dat,tree,alpha.costum)
  
  leaf_LL_Score.BIC$Alpha <- "BIC:log(n) "
  leaf_LL_Score.AIC$Alpha <- "AIC:2"
  leaf_LL_Score.CV$Alpha  <- paste('CV:',alpha.cv,sep = '')
  leaf_LL_Score.CUS$Alpha <- round(alpha.costum,digits = 1)
  leaf_LL_Score <- rbind(leaf_LL_Score.BIC,leaf_LL_Score.AIC,leaf_LL_Score.CV,leaf_LL_Score.CUS)
  pGalpha <- plotLeafScore(leaf_LL_Score)
  return(pGalpha)
}


ChooseAlpha <- function(R.ave,threhold){
  n <- dim(R.ave)[2]
  Rsqure <- 0
  stop <- FALSE
  while((!stop)&(n>=1)){
    if(abs(R.ave[3,n]-Rsqure)>=threhold){
      Rsqure <- R.ave[3,n]
      alpha <- R.ave[1,n]
      n <- n-1
    } else{
      stop <- TRUE
    }
  }
  return(alpha)
}



SigmaToRisk <- function(R.ave){
  R.ave[3,] <- 1/2*log(2*pi)+R.ave[3,]/2+1/2
  return(R.ave)
}


