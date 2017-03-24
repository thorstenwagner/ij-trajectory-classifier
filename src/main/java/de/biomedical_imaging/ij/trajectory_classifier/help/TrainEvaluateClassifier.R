## Author: Thorsten Wagner, University of Applied Sciences and Arts Dortmund
## Classification and Segmentation of Nanoparticle Diffusion Trajectories in Cellular Micro Environments 
## Wagner T, Kroll A, Haramagatti CR, Lipinski HG, Wiemann M (2017) Classification and Segmentation of Nanoparticle Diffusion Trajectories in Cellular Micro Environments. PLOS ONE 12(1): e0170165. 
## doi:  10.1371/journal.pone.0170165
##
## Description:
## Trains and evaluates the random forest classifier for diffusion trajectories
##
## Instructions:
## To apply that script you need the tracks_training.Rata and tracks_test.Rata.
## Both files are available through ResearchGate and have to be in the same
## folder as this script

library("randomForest");
library("caret");
library("irr");
library(plyr)

###############
## Load Data ##
###############


selectedFeatures<-c("TYPES","FD","TRAPPED","EFFICENCY","STRAIGHTNESS","POWER","KURT",
"GAUSS","ASYM3","MSD.R");

load("tracks_training.RData");
data[,1] <- factor(data[,1])
tracks.training <- data;
tracks.training <- tracks.training[,selectedFeatures]

load("tracks_test.RData");
data[,1] <- factor(data[,1])
tracks.test <- data
tracks.test <- tracks.test[,selectedFeatures]
print("Features were loaded");

###############################
## CLASSIFICATION - Training ##
###############################
set.seed(1234);


training.features <- tracks.training[,-1,drop=FALSE];

training.types <- tracks.training[,1];
training.types <- revalue(training.types, c("FREE"="NORM. DIFFUSION","CONFINED+FLOW"="DIRECTED/ACTIVE", "ANOMALOUS+FLOW"="DIRECTED/ACTIVE","FREE+FLOW"="DIRECTED/ACTIVE","ANOMALOUS"="SUBDIFFUSION","CONFINED"="CONFINED","ACTIVE"="DIRECTED/ACTIVE","STALLED"="STALLED"))

test.features <-  tracks.test[,-1,drop=FALSE];
test.types  <- tracks.test [,1];
test.types <- revalue(test.types,  c("FREE"="NORM. DIFFUSION","CONFINED+FLOW"="DIRECTED/ACTIVE", "ANOMALOUS+FLOW"="DIRECTED/ACTIVE","FREE+FLOW"="DIRECTED/ACTIVE","ANOMALOUS"="SUBDIFFUSION","CONFINED"="CONFINED","ACTIVE"="DIRECTED/ACTIVE","STALLED"="STALLED"))

myMtry = 4;
model <- randomForest(training.types ~ ., data = training.features,mtry=myMtry, importance=TRUE)

## Save the model to disk ##
save(model,file="randomForestModel.RData")

#################################
## CLASSIFICATION - Evaluation ##
#################################
print("######################")
print("## MODEL-STATISTICS ##")
print("######################")

print(model);

folds <- createFolds(training.types,k=10);

cs_results_kappa <- lapply(folds,function(x) {
		train <- training.features[x,,drop=FALSE];
		test <- training.features[-x,,drop=FALSE];
		types.training <- training.types[x];
		types.test <- training.types[-x];
		model2 <- randomForest(types.training ~ ., data = train,mtry=myMtry);
		pred <- predict(model2,test);
		actual <- types.test;
		kappa <- kappa2(data.frame(actual,pred))$value;
		cm <- confusionMatrix(pred,actual);
		ac <- as.numeric(cm$overall["Accuracy"]);
		sensDirected <- cm$byClass[1,1]
		specDirected <- cm$byClass[1,2]
		accDirected <- cm$byClass[1,8]
		sensSub <- cm$byClass[2,1]
		specSub <- cm$byClass[2,2]
		accSub <- cm$byClass[2,8]
		sensConf <- cm$byClass[3,1]
		specConf <- cm$byClass[3,2]
		accConf <- cm$byClass[3,8]
		sensNorm <- cm$byClass[4,1]
		specNorm <- cm$byClass[4,2]
		accNorm <- cm$byClass[4,8]


		return(c(kappa,ac,sensDirected,sensSub,sensConf,sensNorm,
				specDirected,specSub,specConf,specNorm,
				accDirected,accSub,accConf,accNorm));
})

print("##############################")
print("## 10-Fold-Cross Validation ##")
print("##############################")

print(paste("Mean Kappa: ", mean(sapply(cs_results_kappa, "[[", 1)) ));
print(paste("SD Kappa: ", sd(sapply(cs_results_kappa, "[[", 1))));
print(paste("Mean Accuracy: ", mean(sapply(cs_results_kappa, "[[", 2))));
print(paste("SD Accuracy: ", sd(sapply(cs_results_kappa, "[[", 2))));
print("######");
print(paste("Mean Directed Sens.: ", mean(sapply(cs_results_kappa, "[[", 3))));
print(paste("Mean Sub. Sens.: ", mean(sapply(cs_results_kappa, "[[", 4))));
print(paste("Mean Conf. Sens.: ", mean(sapply(cs_results_kappa, "[[", 5))));
print(paste("Mean Norm. Sens.: ", mean(sapply(cs_results_kappa, "[[", 6))));
print(paste("Mean Directed Spec.: ", mean(sapply(cs_results_kappa, "[[", 7))));
print(paste("Mean Sub. Spec.: ", mean(sapply(cs_results_kappa, "[[", 8))));
print(paste("Mean Conf. Spec.: ", mean(sapply(cs_results_kappa, "[[", 9))));
print(paste("Mean Norm. Spec.: ", mean(sapply(cs_results_kappa, "[[", 10))));
print(paste("Mean Directed Acc.: ", mean(sapply(cs_results_kappa, "[[", 11))));
print(paste("Mean Sub. Acc.: ", mean(sapply(cs_results_kappa, "[[", 12))));
print(paste("Mean Conf. Acc.: ", mean(sapply(cs_results_kappa, "[[", 13))));
print(paste("Mean Norm. Acc.: ", mean(sapply(cs_results_kappa, "[[", 14))));



print("#####################")
print("## TEST-STATISTICS ##")
print("#####################")
features.predict <- predict(model,test.features);

print(confusionMatrix(features.predict,test.types));



#features.predict.prob <- predict(model,test.features,type="prob");





