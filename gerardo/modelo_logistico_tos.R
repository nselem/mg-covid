library(plotly)
library(glmnet)
library(dendextend)
library(factoextra)
library(ggplot2)
library(cluster)
library(MASS)
library(caret)
library(pROC)

#Leemos el archivo con todos los oligotipos (1405 microbiomas, 230 taxas)
oligos<-read.csv("C:/Users/gerar/Documents/oligotype_count.csv",row.names=NULL)

#Leemos el archivo con los datos de la submuestra de estudio (124 microbiomas)
study_data<-read.csv("C:/Users/gerar/Documents/study_data.csv",row.names=NULL)

#Leemos los oligotipos de la submuestra de estudio
oligos_study=matrix(nrow = length(study_data[,3]),ncol = length(oligos[,1]))
for (i in c(1:length(oligos_study[,1]))){
  for (j in c(1:length(oligos_study[1,]))){
    oligos_study[i,j]<-oligos[j,study_data[i,3]+1]
  }
}

#Normalización de los oligotipos (proporciones)
class(oligos_study)
oligos_study_pr<-oligos_study/rowSums(oligos_study)

##############################################################################################
#ESTUDIO DE LOS MICROBIOMAS Y LA FIEBRE
##############################################################################################

#Leemos los datos de presencia de fiebre
cough_study=study_data[,10]

#Fijamos una semilla
set.seed(100)

#Elegimos conjuntos de entrenamiento y prueba
train_set<-sample(1:length(cough_study), 0.8*length(cough_study))
oligos_study_pr_train<-oligos_study_pr[train_set,]
oligos_study_pr_test<-oligos_study_pr[-train_set,]
cough_study_train<-cough_study[train_set]
cough_study_test<-cough_study[-train_set]

#Datos de estudio
data_train<-data.frame(oligos_study_pr_train)
data_test<-data.frame(oligos_study_pr_test)
cough_study_train<-factor(cough_study_train)
cough_study_test<-factor(cough_study_test)

#############################################################################################
#MODELO LOGISTIC STEPWISE FORWARD
#############################################################################################

#Formula del modelo completo (saturado)
xnam<-paste("X",1:230,sep="")
fmla<-as.formula(paste("cough_study_train~",paste(xnam,collapse = "+")))

#Formulación del modelo logístico (dos clases)
null_model<-glm(formula=cough_study_train~1,family ="binomial",data=data_train)
summary(null_model)
horizonte<-formula(fmla)
model_forward<-stepAIC(null_model,trace=FALSE,direction = "forward",scope=horizonte,criterion="AIC",steps = 10000)
model_forward$anova     #Muestra proceso de seleccion
summary(model_forward)  #Descripción del modelo elegido

#Predicción y matriz de confusion
proba_test<-predict.glm(model_forward,newdata=data_test,type = "response")
confusionMatrix(data=factor(as.numeric(proba_test>0.5)),reference = cough_study_test)

#Curva ROC
curva_roc<-roc(response=cough_study_test,predictor=proba_test)
plot(curva_roc)

##############################################################################################
#MODELO LOGISTICO ELASTIC NET
##############################################################################################

#Eleccion de alfa
alfa=0.1

#Elegimos el valor de lambda por validacion cruzada
cv_net<-cv.glmnet(x=oligos_study_pr_train, y=cough_study_train,type.measure="class",alpha = alfa,family="binomial",nfolds=5)

#Grafica de errores 
plot(cv_net)

#Lambdas encontrados
cv_net

#Elegimos el valor de lambda 
lam<-cv_net$lambda.min

#Modelo logístico con la elección de alfa y lambda anterior
model_net<-glmnet(x=oligos_study_pr_train,y=cough_study_train,family="binomial",alpha=alfa,lambda=lam)
summary(model_net)

#Predicción
cough_predict<-predict.glmnet(model_net,newx=oligos_study_pr_test,type="response")

#Error de clasificacion
x=assess.glmnet(model_net,newx=oligos_study_pr_test,newy=cough_study_test,family="response")$class

#Matriz de confusion
confusion.glmnet(model_net,newx=oligos_study_pr_test,newy=cough_study_test,family="binomial")

#Curva ROC
curva_roc<-roc(response=cough_study_test,predictor=cough_predict)
plot(curva_roc)

##############################################################################################
#COMPARACION DE MODELOS: STEPWISE VS ELASTIC NET
##############################################################################################

#Numero de ensayos
num=100

#Vectores con las precisiones
accuracy_stepwise=vector(length = num)
accuracy_net=vector(length = num)

#Semilla 
set.seed(1)

for (i in c(1:num)){
  #Elegimos conjuntos de entrenamiento y prueba
  train_set<-sample(1:length(cough_study), 0.8*length(cough_study))
  oligos_study_pr_train<-oligos_study_pr[train_set,]
  oligos_study_pr_test<-oligos_study_pr[-train_set,]
  cough_study_train<-cough_study[train_set]
  cough_study_test<-cough_study[-train_set]
  
  #Datos de estudio
  data_train<-data.frame(oligos_study_pr_train)
  data_test<-data.frame(oligos_study_pr_test)
  cough_study_train<-factor(cough_study_train)
  cough_study_test<-factor(cough_study_test)
  
  #############################################################################################
  #MODELO LOGISTIC STEPWISE FORWARD
  #############################################################################################
  
  #Formula del modelo completo (saturado)
  xnam<-paste("X",1:230,sep="")
  fmla<-as.formula(paste("cough_study_train~",paste(xnam,collapse = "+")))
  
  #Formulación del modelo logístico (dos clases)
  null_model<-glm(formula=cough_study_train~1,family ="binomial",data=data_train)
  summary(null_model)
  horizonte<-formula(fmla)
  model_forward<-stepAIC(null_model,trace=FALSE,direction = "forward",scope=horizonte,criterion="AIC",steps = 10000)
  model_forward$anova     #Muestra proceso de seleccion
  summary(model_forward)  #Descripción del modelo elegido
  
  #Predicción y precision
  proba_test<-predict.glm(model_forward,newdata=data_test,type = "response")
  accuracy_stepwise[i]=unname(confusionMatrix(data=factor(as.numeric(proba_test>0.5)),reference = cough_study_test)$overall[1])
  
  ##############################################################################################
  #MODELO LOGISTICO ELASTIC NET
  ##############################################################################################
  #Buscamos el índice de la rejilla en alpha que minimiza el error de prueba en la elastic net
  index_net=1
  pres_net=0
  
  #Elegimos conjuntos de entrenamiento y prueba para elegir alfa y beta para eleastic net (usando solo el conjunto de entrenamiento)
  train_set<-sample(1:length(cough_study_train), 0.9*length(cough_study_train))
  oligos_study_pr_train_net<-oligos_study_pr_train[train_set,]
  oligos_study_pr_test_net<-oligos_study_pr_train[-train_set,]
  cough_study_train_net<-cough_study_train[train_set]
  cough_study_test_net<-cough_study_train[-train_set]
  
  for (j in c(1:20)){
    alfa=j/20
    #Elegimos el valor de lambda por validacion cruzada
    cv_net<-cv.glmnet(x=oligos_study_pr_train_net, y=cough_study_train_net,type.measure="class",alpha = alfa,family="binomial",nfolds=5)
    
    #Elegimos el valor de lambda 
    lam<-cv_net$lambda.min
    
    #Modelo logístico con la elección de alfa y lam anterior
    model_net<-glmnet(x=oligos_study_pr_train_net,y=cough_study_train_net,family="binomial",alpha=alfa,lambda=lam)
    
    #Predicción
    cough_predict<-predict.glmnet(model_net,newx=oligos_study_pr_test_net,type="response")
    
    #Error de clasificacion
    aux=unname(assess.glmnet(model_net,newx=oligos_study_pr_test_net,newy=cough_study_test_net,family="response")$class[1])
    
    #Buscamos el índice del mejor alpha 
    if (pres_net<1-aux){
      index_net=j
      pres_net=1-aux
    }
  }
  
  alfa=index_net/20
  #Elegimos el valor de lambda por validacion cruzada
  cv_net<-cv.glmnet(x=oligos_study_pr_train_net, y=cough_study_train_net,type.measure="class",alpha = alfa,family="binomial",nfolds=5)
  
  #Elegimos el valor de lambda 
  lam<-cv_net$lambda.min
  
  #Modelo logístico con la elección de alfa y lam anterior
  model_net<-glmnet(x=oligos_study_pr_train_net,y=cough_study_train_net,family="binomial",alpha=alfa,lambda=lam)
  
  #Predicción para el conjunto de prueba original
  cough_predict<-predict.glmnet(model_net,newx=oligos_study_pr_test,type="response")
  
  #Error de clasificacion
  accuracy_net[i]=1-unname(assess.glmnet(model_net,newx=oligos_study_pr_test,newy=cough_study_test,family="response")$class[1])
}

#Resumen de los errores del método stepwise
summary(accuracy_stepwise)
hist(accuracy_stepwise)

#Resumen de los errores por elastic net
summary(accuracy_net)
hist(accuracy_net)

##################################################################################################
#TAXAS SELECCIONADOS POR ELASTIC NET
##################################################################################################

#Numero de ensayos
num=100

#Semilla 
set.seed(1)

#Vector con los coeficientes seleccionados
coef_select=replicate(230,0)
num_coef=replicate(num,0)
alfa_select=vector(length = num)
lambda_select=vector(length = num)

for (i in c(1:num)){
  #Elegimos conjuntos de entrenamiento y prueba (con el fin de estimar los alfa y lambda óptimos)
  train_set<-sample(1:length(cough_study), 0.8*length(cough_study))
  oligos_study_pr_train<-oligos_study_pr[train_set,]
  oligos_study_pr_test<-oligos_study_pr[-train_set,]
  cough_study_train<-cough_study[train_set]
  cough_study_test<-cough_study[-train_set]
  
  ##############################################################################################
  #MODELO LOGISTICO ELASTIC NET
  ##############################################################################################
  index_net=1
  pres_net=0
  
  for (j in c(1:20)){
    alfa=j/20
    #Elegimos el valor de lambda por validacion cruzada
    cv_net<-cv.glmnet(x=oligos_study_pr_train, y=cough_study_train,type.measure="class",alpha = alfa,family="binomial",nfolds=5)
    
    #Elegimos el valor de lambda 
    lam<-cv_net$lambda.min
    
    #Modelo logístico con la elección de alfa y lam anterior
    model_net<-glmnet(x=oligos_study_pr_train,y=cough_study_train,family="binomial",alpha=alfa,lambda=lam)
    
    #Predicción
    cough_predict<-predict.glmnet(model_net,newx=oligos_study_pr_test,type="response")
    
    #Error de clasificacion
    aux=unname(assess.glmnet(model_net,newx=oligos_study_pr_test,newy=cough_study_test,family="response")$class[1])
    
    #Buscamos el índice del mejor alpha 
    if (pres_net<1-aux){
      index_net=j
      pres_net=1-aux
    }
  }
  alfa=index_net/20
  #Elegimos el valor de lambda por validacion cruzada
  cv_net<-cv.glmnet(x=oligos_study_pr_train, y=cough_study_train,type.measure="class",alpha = alfa,family="binomial",nfolds=5)
  
  #Elegimos el valor de lambda 
  lam<-cv_net$lambda.min
  
  #Modelo logístico con la elección de alfa y lam anterior
  model_net<-glmnet(x=oligos_study_pr_train,y=cough_study_train,family="binomial",alpha=alfa,lambda=lam)
  
  #Variables seleccionadas
  for (k in c(1:length(model_net$beta))){
    if (model_net$beta[k]!=0){
      coef_select[k]=coef_select[k]+1
      num_coef[i]=num_coef[i]+1
    } 
  }
  #Guardamos los alfas y lambdas
  alfa_select[i]=alfa
  lambda_select[i]=lam
}

#Numero de variables elegidas
hist(num_coef,main="Numero de variables seleccionadas")

#Alfas elegidas
hist(alfa_select,main="Alfas",breaks = 10)

#Lambdas elegidas
hist(lambda_select,main="Lambdas")

#Variables elegidas
barplot(coef_select,main="Selección de variables (frecuencias)")
