library(vegan)
library(ggpubr)
library(knitr)
library(plotly)
library(phyloseq)
library(microbiome)
library(DirichletMultinomial)
library(magrittr)
library(dplyr)
library(reshape2)
library(plotly)
library('plot.matrix')
library(glmnet)
library(h2o)


# Importa datos.
datos_incompletos <- read.csv("datos.csv", row.names=NULL)

# Limpia NAs y filas incompletas
datos <- datos_incompletos[ complete.cases(datos_incompletos), ]

# Síntoma:
# 6: Fiebre.
# 7: Escurrimiento.
# 8: Dolor de Garganta.
# 9: Tos.
sint <- 7

# Covariables: n x p = 230 especies
X <- datos[,10:239]

# Abundancia relativa.
R <- X / rowSums(X)
R <- data.matrix(R)

# Con abundancia total
#R <- data.matrix(X)

# Fija semilla.
set.seed(1)

# Núm. de observaciones.
n <- length(datos[,1])

# Filas de entrenamiento
train <- sample( 1:n, 0.85*n )
R_train <- R[ train, ]
R_test <- R[ -train, ]

# Totales.
Y <- factor(datos[,sint])

# Estimación.
Y_train <- Y[ train ]

# Prueba.
Y_test <- Y[ -train ]

# Frames para h2o.
data_set <- cbind(R, datos[,sint])
colnames(data_set) <- c(colnames(R),"Y")
train_set <- data_set[ train, ]
test_set <- data_set[ -train, ]

# Hiperparámetros de los modelos de búsqueda
alpha_busqueda <- list()
lambda_busqueda <- list()

#######################################################
# Modelos en h2o para encontrar hiperparámetros       #
# candidatos, que se verifican luego con glmnet.      #
#######################################################

h2o.init()

# Asigna formato h2o
train_h2o <- as.h2o(train_set)
test_h2o <- as.h2o(test_set)

## Rejilla de hiperparámetros
hiper <- list( 
  lambda = exp( seq( log(0.001), 0, length.out = 250) ),
  alpha = seq(0,1,length.out = 20)
)

#######################################################
# Correr la parte de rejilla unas 2 o 3 de veces.     #
# Cada una de ellas, guardar los 3 primeros juegos de #
# hiperparámetros lambda y alpha que da summary.      #
# Los hiperparámetros son de distintas corridas con   #
# diferentes semillas. Una vez que se tengan varios   #
# juegos, tomar promedio de ALPHA para glmnet, correr #
# búsqueda de  lambda y verificar                     #
#######################################################

for (i in seq(1:10)) {
  # Ajuste de hiperparámetros
  rejilla <- h2o.grid(
    algorithm = "glm",
    y="Y",
    x=colnames(R),
    #seed = 1,
    training_frame = train_h2o,
    family = "binomial",
    #family = "fractionalbinomial",
    hyper_params = hiper,
    intercept = FALSE,
    search_criteria = list( 
      strategy = "RandomDiscrete",
      #strategy = "Cartesian",
      stopping_metric = "misclassification", 
      stopping_tolerance = 0.00001, 
      stopping_rounds = 5)
  )
  
  #####################################################
  # Aquí se guardan los tres mejores hiperparámetros  #
  # (mínimo error de clasificación) de las            #
  # regresiones según la rejilla propuesta para cada  #
  # una de las simulaciones.                          #
  #                                                   #
  # AQUÍ TAMBIÉN ES POSIBLE OBTENER LAS ESPECIES      #
  # QUE MEJOR PREDICEN EN CADA SIMULACIÓN.            #
  #                                                   #
  # Más adelante el programa las taxa que mejor       #
  # predicen, 'tmp', de una sóla regresión, la de     #
  # los mejores hiperparámetros; PERO tal vez sea     #
  # mejor ese conjunto 'tmp' a partir de ésta parte   #
  # del algoritmo, y a ese ajustar con los            #
  # mejores hiperparámtros.                           #
  #####################################################
  
  alpha_busqueda[[i]] <- rejilla@summary_table$alpha[1:3]
  lambda_busqueda[[i]] <- rejilla@summary_table$lambda[1:3]
}

# Cierra H2o.
h2o.shutdown(prompt = FALSE)

# Hiperparámetros (deshace lista)
alpha_busqueda <- unlist(alpha_busqueda)
lambda_busqueda <- unlist(lambda_busqueda)

# Hiperparámetros (quita brackets)
alpha_busqueda <- sub('\\[(.*)\\]','\\1',alpha_busqueda)
lambda_busqueda <- sub('\\[(.*)\\]','\\1',lambda_busqueda)

# Hiperparámetros (caracter a número)
alpha_busqueda <- as.numeric(alpha_busqueda)
lambda_busqueda <- as.numeric(lambda_busqueda)

#######################################################
# Parte de glmnet para validación cruzada.            #
#######################################################

alpha_modelo <- median(alpha_busqueda)
lambda_seq <- seq( min(lambda_busqueda), max(lambda_busqueda), length.out = 100)

# Hiperparámetro por Validación cruzada.
cv_net <- cv.glmnet(
  x=R_train, y=Y_train, 
  type.measure = "class", 
  standarize = FALSE,
  lambda = lambda_seq,
  alpha = alpha_modelo, 
  family = "binomial",
  intercept = FALSE,
  seed = 1,
  #nfolds = 10
  )

# Gráfica para CV de lambda
plot(cv_net)

# Lambda de menor error.
lambda_modelo <- cv_net$lambda.min

# Ajusta modelo.
net <- cv_net$glmnet.fit

# Comparación.
assess.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = "binomial", 
  s = lambda_modelo )$class

confusion.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = "binomial", 
  s = lambda_modelo )

# Coeficientes no-nulos
coeff_net <- coef.glmnet(net,lambda_modelo)
sum( coeff_net != 0 )
sum( coeff_net != 0 ) / length( coeff_net )

#######################################################
# Importancia de los parámetros (porcentaje que el    #
# parámetro contribuye a la norma 2 total del vector  #
# de coefficientes).                                  #
#######################################################

escala <- coeff_net^2
magnitud <- sum( escala ) 
porcent <- escala / magnitud

# Importancia de variables (individual)
iv <- data.frame( as.matrix(porcent) )
colnames(iv) <- c("Porcentaje")
iv <- iv[ order(-iv$Porcentaje), , drop = FALSE ]

# Acumulado
iv$Acumulado <- cumsum(iv$Porcentaje)


# Compocisión porcentual
plot_ly( iv,
         x= ~reorder(rownames(iv),-Porcentaje),
         y= ~Porcentaje*100,
         type = "bar") %>%
  layout(yaxis = list(title = 'Porcentaje explicado (%)'))

# Acumulado.
plot_ly( iv,
         x= ~reorder(rownames(iv),-Porcentaje),
         y= ~Acumulado*100,
         type = "bar") %>%
  layout(yaxis = list(title = 'Porcentaje acumulado explicado (%)'))

# Taxa más predictora.
tmp <- rownames( iv[iv$Porcentaje>0,] )
tmp <- sub('X(.*)','\\1',tmp)


#######################################################
# Comienza clustering.                                #
#######################################################


# OTUs: lee datos -> matriz numérica -> objeto de phyloseq
tabla_otu <- read.csv("oligotype_count.csv", row.names = "X")
otus <- data.matrix(tabla_otu)

# Quita otus fuera de tmp.
# Aquí utiliza TODOS los oligotipos de la tabla de conteo,
# de manera que clasifica todas las observaciones, sin
# importar si son pacientes o no.
otus <- otus[tmp,]
OTU <- otu_table(otus, taxa_are_rows = TRUE)

# Taxas: lee datos -> data.frame -> objeto de phyloseq
tabla_taxa <- read.csv("tax_data.csv",row.names = 1)
taxa <- as.matrix(tabla_taxa)
taxa <- data.frame(taxa[as.array(tmp),])
colnames(taxa) <- c("Taxonomy")
taxa <- as.matrix(taxa)
TAX <- tax_table(taxa)

# Crea elemento de la clase phyloseq
physeq <- phyloseq(OTU, TAX)

#Abundancia relativa.
physeq.comp <- microbiome::transform(physeq, "compositional")

# Taxa del core
taxa <- core_members(physeq.comp, detection = 0.0005/100, prevalence = 100/100)
#physeq <- prune_taxa(taxa, physeq)

# Conteo de OTU (con la info descartada)
dat <- abundances(physeq)
# Formato muestra x taxa
count <- as.matrix(t(dat))

#######################################################
# Ajusta Modelos de clustering.                       #
# 'load' ESTÁ DESCOMENTADO PARA CHECAR EL CÓDIGO.     #
# COMENTAR 'load' Y DESCOMENTAR 'fit' PARA PRUEBAS    #
# REALES.                                             #
#######################################################

set.seed(1234)
#fit <- lapply(1:8, dmn, count = count, verbose=TRUE)
load(file = "Modelo-propuesto.RData")

# Bondad de ajuste
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Número de Clusters", ylab="Ajuste de modelo")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

# Escoger MANUALMENTE mejor modelo de los que calculó.
n_cluster <- 7
best <- fit[[ n_cluster ]]

# Parámetros del mejor modelo.
mixturewt(best)
pi_modelo <- mixturewt(best)[,1]

# Clasificación de cada paciente
gpos <- apply( mixture(best), 1, which.max )
#write.csv(gpos,"grupos_r100.csv")

# Pacientes por grupo (aproximado)
round( best@mixture$Weight * length(tabla_otu) )

# Aporte de cada OTU en el cluster.
porc <- fitted(best,scale=TRUE)
#write.csv(porc,"aporte100.csv")


#######################################################
# Exploración gráfica.                                #
#######################################################

# Gráfica.
plot_ly( x = ~gpos, type = "histogram" )

# Mapa de calor
#heatmapdmn(count,fit[[1]],best,ntaxa = 15)
#dev.off()

# Nombra los OTUs de cada cluster.
datos_sn <- melt( porc )
colnames(datos_sn) <- c("otu","grupo","valor")
nombres <- cbind(tabla_taxa,row.names(tabla_taxa))
colnames(nombres) <- c("taxa","otu")

# Deja taxa mejor predictora.
nombres <- nombres[tmp,]

datos_n <- merge(datos_sn,nombres,by="otu")
datos_n <- datos_n[ order(datos_n$grupo,-datos_n$valor), ]



# Taxa más prevalecientes en cada grupo
rep <- as.matrix(datos_n[datos_n$grupo==1,"otu"], ncol=1)
for( k in seq(2,n_cluster) ){
  rep <- cbind(rep,datos_n[datos_n$grupo==k,"otu"])
}

# GRafica todo junto.
graf <- ggplot(datos_n, aes(x=grupo, y=valor) ) +
  geom_bar(stat="identity", aes(fill=otu), position = "stack")
ggplotly(graf)


#######################################################
# Califica ponderando los coeficientes de pi y los    #
# lugares en cada grupo. CHECAR MANUALMENTE QUE 'vec' #
# funcione correctamente SEGÚN EL NÚMERO DE CLUSTERS. #
#######################################################

# Calificación de taxa.
n_rep <- length(rep[,1])

# Ajuste de las calificaciones según número de Clusters.
vec <- list()
for (i in seq(1,n_cluster)) {
  vec[[i]] <- n_rep*(i-1)
}
vec <- as.array( unlist(vec) )

# Función calificadora.
calif <- function(x) as.numeric( 
  pi_modelo %*% (
    which(rep==x) - vec
    ) 
  )

calificaciones <- apply(as.array(nombres[tmp,]$otu), 1, calif)
nombres <- cbind(nombres, calificaciones )
nombres <- nombres[ nombres$calificaciones > 0, ]
nombres <- nombres[ order( nombres$calificaciones), ]
#nombres[1:15, c("taxa")]
#write.csv(nombres,"clusters100.csv")


# Composición de los 20 mejores clasificados.
mejores <- row.names( nombres[1:20,] )
ind <- function(i) row.names(datos_n[ datos_n$otu == i, ])
mejores_ind <- apply(as.array(mejores), 1, ind)
mejor <- datos_n[ c(mejores_ind), ]

# GRafica mejores.
graf <- ggplot(mejor, aes(x=grupo, y=valor) ) +
  geom_bar(stat="identity", aes(fill=otu), position = "stack")
ggplotly(graf)

# Tabla de phylosecuencias clusterizada
sample_data(physeq) <- data.frame(gpos)

#######################################################
# Alfa diversidad                                     #
#######################################################

alfa <- microbiome::alpha(physeq, index = "all")
#kable(head(alfa))

# Metadatos para gráfica
m_physeq <- meta(physeq)
m_physeq$shannon <- alfa$diversity_shannon
m_physeq$chao1 <- alfa$chao1
kable(head(m_physeq))


# Gráfica de violín (shannon) con plotly.
fig <- m_physeq %>%
  plot_ly(
    x = ~gpos,
    y = ~shannon,
    split = ~gpos,
    type = 'violin',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  ) 
fig <- fig %>%
  layout(
    xaxis = list(
      title = "Grupos"
    ),
    yaxis = list(
      title = "Shannon",
      zeroline = F
    )
  )
fig

# Pruebas (Shannon). SE USA AJUSTE POR BENJAMINI.
kruskal.test( shannon ~ gpos, data = m_physeq )
pairwise.wilcox.test( m_physeq$shannon, m_physeq$gpos, p.adjust.method = "BH" ) 


# Gráfica de violín (chao1) con plotly.
fig <- m_physeq %>%
  plot_ly(
    x = ~gpos,
    y = ~chao1,
    split = ~gpos,
    type = 'violin',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  ) 
fig <- fig %>%
  layout(
    xaxis = list(
      title = "Grupos"
    ),
    yaxis = list(
      title = "Chao1",
      zeroline = F
    )
  )
fig

# Pruebas (Chao1). SE USA AJUSTE POR BENJAMINI
kruskal.test( chao1 ~ gpos, data = m_physeq )
pairwise.wilcox.test( m_physeq$chao1, m_physeq$gpos, p.adjust.method = "BH" ) 

#######################################################
# Beta diversidad                                     #
#######################################################

clusters <- as.list(1:n_cluster)
bray <- as.list(1:n_cluster)
jaccard <- as.list(1:n_cluster)


for (k in seq(1:n_cluster)) {
  eval(
    parse(
      text=sprintf(
        'clusters[[%s]] <- subset_samples(physeq, sample_data(physeq)$gpos == %s)', 
        k, 
        k
      )
    )
  )
  
  # Referencia a la mediana.
  ref <- apply(abundances(clusters[[k]]), 1, median)
  
  # Bray.
  bray[[k]] <- array( as.numeric( 
    unlist( divergence(clusters[[k]], ref, method = "bray") )
  ) )
  
  # Jaccard.
  jaccard[[k]] <- array( as.numeric( 
    unlist( divergence(clusters[[k]], ref, method = "jaccard") )
  ) )
  
}


# Grafica Bray.
fig <- plot_ly(
  y = bray[[1]], 
  type = "box", 
  quartilemethod="linear",
  name="Grupo 1"
)

for (k in c(seq(2:n_cluster) + 1) ) {
  eval(
    parse(
      text=sprintf(
        'fig <- fig %%>%% add_trace( y = bray[[%s]], quartilemethod="linear", name="Grupo %s")', 
        k, 
        k
      )
    )
  )
}

fig <- fig %>% layout(title = "Disimilaridad Bray-Curtis")
fig

# Pruebas Bray-Curtis.
bc_dist <- phyloseq::distance(physeq,method = "bray")
met <- data.frame(sample_data(physeq))
adonis( bc_dist ~ factor(gpos), data = met)


# Grafica Jaccard.
fig <- plot_ly(
  y = jaccard[[1]], 
  type = "box", 
  quartilemethod="linear",
  name="Grupo 1"
)

for (k in c(seq(2:n_cluster) + 1) ) {
  eval(
    parse(
      text=sprintf(
        'fig <- fig %%>%% add_trace( y = jaccard[[%s]], quartilemethod="linear", name="Grupo %s")', 
        k, 
        k
      )
    )
  )
}

fig <- fig %>% layout(title = "Disimilaridad Jaccard")
fig

# Jaccard.
jac_dist <- phyloseq::distance(physeq,method = "jaccard")
adonis( jac_dist ~ factor(gpos), data = met)


#######################################################
# Regresión sobre los grupos.                         #
#######################################################

# Data frame para regresión.
data_reg <- datos[,c(2,sint)]

# ID de la muestra para pacientes del estudio.
#set.seed(5678)
g_id <- function( id ) gpos[id]
sample_gps <- apply( as.array(data_reg$SampleID) , 1, g_id )
data_reg$gps <- factor(sample_gps)

# Muestras de estimación y validación
train_idx <- apply(
  as.array(datos$SampleID[ train ]), 
  1, 
  function( id ) which( data_reg$SampleID == id )
  )
train_gps <- data_reg[train_idx,]
test_gps <- data_reg[-train_idx,]

# Modelo
reg <- glm( train_gps[,2] ~ gps, 
            data = train_gps, 
            family=binomial(link="logit")
            )

summary(reg)

# Predicción
prob <- predict.glm(reg,test_gps,type = "response")
pred <- rep(0,length(test_gps[,1]))
pred[ prob > 0.5 ] <- 1

# Validación
valida <- data.frame( prob, pred, test_gps[,2] )
valida$Coincidencia <- equals( pred[], test_gps[,2] )
colnames(valida) <- c("Probabilidad","Predicción","Real","Coincidencia")

# Resumen y Porcentaje de aciertos
kable(valida)
countMatches(TRUE,valida$Coincidencia)/length(valida[,1])














