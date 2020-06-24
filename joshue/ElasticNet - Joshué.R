library(plotly)
library(glmnet)
library('plot.matrix')

# Pone semilla.
set.seed(1234)

# Importa datos.
datos <- read.csv("~/CIMAT/Materias/Doctorado/Semestre 2/Taller Multidisciplinario/Influenza/datos.csv", row.names=NULL)

# Respuesta: n = 124 pacientes. Las variables binarias son
#   6.  fever
#   7.  rhino
#   8.  throat
#   9.  cough
Y <- factor(datos[,6])
n <- length(Y)

# Covariables: n = 117 x p = 230 especies
X <- datos[,10:239]
class(X)

# Abundancia relativa.
R <- X / rowSums(X)
R <- data.matrix(R)

# Filas de entrenamiento
train <- sample( 1:n, 0.66*n )
R_train <- R[ train, ]
R_test <- R[ -train, ]

Y_train <- Y[ train ]
Y_test <- Y[ -train ]

# Elastic-net: Ridge 0 < ---- > 1 Lasso
alpha = 0.350

# Hiperparámetro por Validación cruzada.
cv_net <- cv.glmnet(
  x=R_train, y=Y_train, 
  type.measure = "class", 
  standarize = FALSE,
  alpha = alpha, 
  family = "binomial", 
  nfolds = 5 )

plot(cv_net)

# Factor de CV. ( aprox exp(-4.3) )
lambda <- cv_net$lambda.min

# Ajuste de modelo.
net <- glmnet( 
  x= R_train, y = Y_train, 
  family = "binomial", 
  standarize = FALSE,
  alpha = alpha, 
  lambda = lambda )

# Coeficientes no cero.
coef(net,lambda)
coef_net <- net$beta[,1]

# Predicción
Y_predict <- predict( 
  net, 
  newx = R_test, 
  s = lambda, 
  type = "response" )

# Comparación.
assess.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = "binomial", 
  s = lambda )$class

confusion.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = "binomial", 
  s = lambda )

# Porcentaje de Coeficientes no-nulos
sum( coef_net != 0 ) / length( coef_net  )

# Nombres
tax_data <- read.csv("~/CIMAT/Materias/Doctorado/Semestre 2/Taller Multidisciplinario/Influenza/tax_data.csv",row.names = 1)
tax_data <- cbind(tax_data, coef_net)

x <- as.array(tax_data[,1])
y <- as.array(tax_data[,2])
z <- as.numeric( y != 0 )

tax_data <- cbind(tax_data, z)

plot_ly(tax_data, x = ~x, y = ~y, type = 'bar')

plot_ly(tax_data, x = ~x, y = ~z)

orden <- tax_data[ order( abs(tax_data$coef_net), decreasing = TRUE ), ]
orden[1:15,1]


