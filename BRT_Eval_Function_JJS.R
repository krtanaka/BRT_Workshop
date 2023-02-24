library(sqldf)
library(dplyr)
library(caret)
library(mlbench)
library(Metrics)
library(dismo)
library(gbm)

#the bulk of the following code was written with Heather Welch, Megan Cimino, and Justin Suca

#this is a balanced cross validation for PA, meaning it splits into equal prop occurrence
#between the training and test sets
#same concept as the ensemble random forest selection
eval_CV_BRT_Balanced<- function(dataInput, gbm.x, gbm.y, lr=lr, tc=tc, family=family,bag.fraction=bag.fraction, n.folds=n.folds){
  DataInput <- dataInput
  if(family=="bernoulli"){
    DataInput_bound_Positive  <- floor((nrow(DataInput[DataInput$gbm.y==1,])/4)*3) #define % of training and test set, by pres/abs
    DataInput_bound_Zero <- floor((nrow(DataInput[DataInput$gbm.y==0,])/4)*3) #define % of training and test set, by pres/abs
    DataInput_Positive<-DataInput[DataInput$gbm.y==1,]
    DataInput_Zero<-DataInput[DataInput$gbm.y==0,]
    DataInput_train_Positive<- DataInput_Positive[sample(nrow(DataInput_Positive),DataInput_bound_Positive),]
    DataInput_train_Zero<- DataInput[sample(nrow(DataInput_Zero),DataInput_bound_Zero),]
    DataInput_train<-rbind(DataInput_train_Positive, DataInput_train_Zero)
    DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')#this pulls the unsampled group
    #the function that really runs/generates the BRT 
    
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family=family, tree.complexity=tc,
                                 learning.rate = lr, bag.fraction=bag.fraction, n.folds=n.folds) 
    #predict the test set to evaluate performance
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    #extracting the observations
    d <- cbind(obs=DataInput_test[,c(gbm.y)], preds)
    
    pres <- d[d[,1]==1, 2]
    abs <- d[d[,1]==0, 2]
    #this functions gets TSS components and AUC
    e <- evaluate(p=pres, a=abs)
    e
    return(list(e, DataInput.kfolds))}
  
  else if (family=="poisson" | family=="gaussian"){
    
    #same thing but for abundance data, no need to rare occurrence splitting
    DataInput_bound <- floor((nrow(DataInput)/4)*3) #define % of training and test set
    DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
    DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family=family, tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = bag.fraction, n.folds=n.folds) #this was gbm.step
    
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    d <- cbind(obs=DataInput_test[,c(gbm.y)], preds)
    
    #currently using R2 and RMSE as model fit metrics
    R2<-cor(d[,1],d[,2])^2
    RMSE<-rmse(d[,1],d[,2])
    e<-cbind(R2,RMSE)
    colnames(e)<-c("R2","RMSE")
    e
    return(list(e, DataInput.kfolds))}
  
}

#the following is the same function, just using a fixed number of trees
#instead of using gbm step to evaluate 
eval_CV_BRT_Balanced_Fixed<- function(dataInput, gbm.x, gbm.y, lr=lr, tc=tc, family=family, nt=nt, bag.fraction=bag.fraction){
  DataInput <- dataInput
  if(family=="bernoulli"){
    DataInput_Positive_Loc<- which(DataInput[,c(gbm.y)]==1)
    DataInput_Positive<-DataInput[DataInput_Positive_Loc,]
    DataInput_Zero_Loc<- which(DataInput[,c(gbm.y)]==0)
    DataInput_Zero<-DataInput[DataInput_Zero_Loc,]
    DataInput_bound_Positive  <- floor((nrow(DataInput_Positive)/4)*3) #define % of training and test set, 75/25 here
    DataInput_bound_Zero <- floor((nrow(DataInput_Zero)/4)*3) #define % of training and test set
    DataInput_train_Positive<- DataInput_Positive[sample(nrow(DataInput_Positive),DataInput_bound_Positive),]
    DataInput_train_Zero<- DataInput_Zero[sample(nrow(DataInput_Zero),DataInput_bound_Zero),]
    DataInput_train<-rbind(DataInput_train_Positive, DataInput_train_Zero)
    DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')
    DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                  family=family, tree.complexity=tc,learning.rate = lr, n.trees = nt, bag.fraction=bag.fraction)
    
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    d <- cbind(obs=DataInput_test[,c(gbm.y)], preds)
    
    pres <- d[d[,1]==1, 2]
    abs <- d[d[,1]==0, 2]
    e <- evaluate(p=pres, a=abs)
    e
    return(list(e, DataInput.kfolds))}
  
  else if (family=="poisson"| family=="gaussian"){
    DataInput_bound <- floor((nrow(DataInput)/4)*3) #define % of training and test set
    DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
    DataInput_test<- sqldf('SELECT * FROM DataInput EXCEPT SELECT * FROM DataInput_train')
    DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                  family=family, tree.complexity=tc,
                                  learning.rate = lr, n.trees=nt, bag.fraction=bag.fraction) 
    
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    d <- cbind(obs=DataInput_test[,c(gbm.y)], preds)
    
    R2<-cor(d[,1],d[,2])^2
    RMSE<-rmse(d[,1],d[,2])
    e<-cbind(R2,RMSE)
    colnames(e)<-c("R2","RMSE")
    e
    return(list(e, DataInput.kfolds))}
  
}

#Function that generates an ensemble of models using iterative splits 
fit.brt.n_eval_Balanced <- function(data, gbm.x, gbm.y, lr,tc,family, bag.fraction=bag.fraction, n.folds=n.folds, iterations){
  Species_Models <- vector("list",iterations)
  Species_Models_Eval <- vector("list",iterations)
  
  for (i in 1:iterations){
    model <- eval_CV_BRT_Balanced(data=data, gbm.x= gbm.x, gbm.y=gbm.y,lr=lr, tc=tc, family=family, bag.fraction=bag.fraction, n.folds=n.folds)
    
    Species_Models[[i]] <- model[[2]]    
    Species_Models_Eval[[i]]<-model[[1]]
  }
  return(list(Species_Models, Species_Models_Eval))
}

fit.brt.n_eval_Balanced_Fixed <- function(data, gbm.x, gbm.y, lr,tc, family, nt, bag.fraction=bag.fraction,iterations){
  Species_Models <- vector("list",iterations)
  Species_Models_Eval <- vector("list",iterations)
  
  for (i in 1:iterations){
    model <- eval_CV_BRT_Balanced_Fixed(data=data, gbm.x= gbm.x, gbm.y=gbm.y,lr=lr, tc=tc, nt=nt, family=family, bag.fraction=bag.fraction)
    
    Species_Models[[i]] <- model[[2]]    
    Species_Models_Eval[[i]]<-model[[1]]
  }
  return(list(Species_Models, Species_Models_Eval))
}
