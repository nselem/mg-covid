---
title: "Reporte 2 - Joshué Ricalde"
author: "Joshué Helí Ricalde Guerrero"
date: "1/7/2020"
output: 
  html_document:
  toc: true
runtime: shiny
---



# Introducción.

\par
En este reporte se presentan los avances realizados durante las 2 primeras semanas del taller multidisciplinario. El documento se estructura de la siguiente manera: Se comienza explicando el trabajo realizado con el modelo personal propuesto, que relaciona directamente especies de taxa con sintomatología viral, usando como datos los del artículo de Lee et. al. y se obtienen del repositorio  https://deepblue.lib.umich.edu/data/concern/generic_works/sb3979224?locale¼en. En la sección posterior, se retoma el mismo artículo y se procede a replicar los resultados, enfocándose únicamente en el problema de detección de comunidades de Taxa. La última sección corresponde a las Conclusiones y trabajo futuro. 


# Regresión Penalizada (Modelo original)

\par
A continuación se presenta una propuesta para relacionar la abundancia relativa del microbioma con los síntomas de enfermedades virales. Se asume que el problema es de alta dimensionalidad, ya que espera que el número de observaciones $n$ sea menor que el número de parámetros $p$ (el número de especies) en el microbioma. Por tal motivo, el modelo propuesto consiste en aplicar técnicas de \emph{Regresión Logística Penalizadas}, específicamente con \emph{Elastic-Net} o *Red Elástica*. 

\par
La razón de aplicar éste tipo de regresión es porque no sólo actua como un predictor de posibles síntomas, sino que por la penalización misma, este tipo de modelo puede aplicarse como un método de selección de variables, pues la restricción en norma $l_1$ y $l_2$ evitan que todos los coeficientes (cada uno asociado a un OTU) crezcan, de manera que únicamente sobrevien aquellos que realmente aportan información al modelo. Así, de poder implementar correctamente este tipo de regresión, se estaría reduciendo la dimensión del problema, considerando únicamente los OTU cuyo coeficiente sea no nulo, al mismo tiempo que se tendría un clasificador "elemental", gracias a la parte logística.


## Datos

Primero, se importan las librerías y datos correspondientes. Se aprovecha para señalar que la razón de haber elegido los datos de ésta manera es porque parece ser que los datos reales que se obtendrán serán similares a estos.

```r
# Las primeras 4 para gráficos, y la última para regresión
library(knitr)
library(reshape2)
library(plotly)
library('plot.matrix')
library(glmnet)

# Importa datos.
datos_incompletos <- read.csv("datos.csv", row.names=NULL)

# Limpia NAs
datos <- datos_incompletos[ complete.cases(datos_incompletos), ]
```


\par
Cada fila de la variable `datos` contiene una de las $n$ observaciones, mientras que las colmunas se ordenan de la siguente manera:

 1. ID del participante (alfanumérica).

 2. ID de la muestra (alfanumérica).

 3. Grupo al que pertenece su micriobioma (categórica, se pretende corroborar que el modelo presenta resultados similares).

 4. Tipo de influenza según PCR (categórica).

 5. Tipo de influenza según Hemaglutinación (categórica).

 6. Presencia de Fiebre (Binaria).

 7. Presencia de Rinorrea o escurrimiento nasal (Binaria).

 8. Presencia de dolor de garganta (Binaria).

 9. Presencia de Tos (Binaria).

El resto de las columnas, 10 a 239, corresponde a cada una de las especies observadas y su abundancia total (numérica).

\par
Antes de pasar a las regresiones en sí, se preporcesan los datos, de manera que se estén trabajando con *Abundancias relativas*. En este caso, se denotará por $X$ a la matriz con abundancias relativas ($n$ observaciones $\times p$ especies) y por $R$ a la matriz $X$ normalizada por las filas.

```r
# Covariables: n x p = 230 especies
X <- datos[,10:239]

# Abundancia relativa.
R <- X / rowSums(X)
R <- data.matrix(R)

# Con abundancia total
#R <- data.matrix(X)
```

\par
**IMPORTANTE**: una observación que se hizo durante las reuniones es la de que las variables independientes en el modelo propuesto son *porcentajes*; entonces, la regresión logística usual no aplíca, sino que *deben hacerse modificaciones*. Al momento de escritura, no se ha podido solventar este problema, pero se reportan los avances tenidos hasta ahora.

\par
Después, se procede a seleccionar la muestra de estimación ($R_{train}$) y la de prueba ($R_{test}$). Para la primera, se toma el $90\%$ de la muestra total, y para la segunda el porcentaje restante.

```r
# Fija semilla.
set.seed(1)

# Nùm. de observaciones.
n <- length(datos[,1])

# Filas de entrenamiento
train <- sample( 1:n, 0.90*n )
R_train <- R[ train, ]
R_test <- R[ -train, ]
```


## Síntomas individuales (Elastic net)

### Fiebre

\par
Se comienza aplicando una regresión a cada uno de los síntomas por separado. Se detalla el procedimiento para el síntoma de Fiebre; el resto es análogo. En este caso, se denota por $Y$, $Y_{train}$ e $Y_{test}$ al vector total de respuestas, respuestas de entrenamiento y respuestas de prueba, respectivamente. También se escoge la familia de regresión como `binomial` al tratarse de una variable respuesta categórica de dos niveles (binaria).



```r
# Familia de regresión.
familia <- "binomial"

# Totales.
Y <- factor(datos[,6])

# Estimación.
Y_train <- Y[ train ]

# Prueba.
Y_test <- Y[ -train ]
```


\par
Se fija el hiperparámetro $\alpha = 0.15$ para red elástica, ya que durante las pruebas se observó que el aporte de la norma $l_1$ (penalización Lasso) era demasiado severo, de manera que en muchas simulaciones, *todos* los coeficientes de regresión eran cero. Después se procede a realizar el procedimiento de *Validación Cruzada*, para hallar el parámetro de penalización adecuada.

```r
# Elastic-net: Ridge 0 < ---- > 1 Lasso
alpha = 0.150

# Fija semilla.
set.seed(2)

# Número de CV
no_cv <- 250

# Hiperparámetro por Validación cruzada.
cv_net <- cv.glmnet(
  x=R_train, y=Y_train, 
  type.measure = "class", 
  standarize = TRUE,
  alpha = alpha, 
  family = familia, 
  nlamda = no_cv,
  #lambda = rejilla,
  #type.multinomial = "grouped",
  #nfolds = 5,
  )

plot(cv_net)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

\par
La gráfica anterior muestra los cambios en el error de clasificación (eje $y$) que se tienen al incrementar el parámetro de penalización $\lambda$ (eje $x$). Cabe resaltar que la gráfica también muestra (en la parte superior) a cuantos coeficientes *no nulos* corresponde cada valor de $\lambda$.


\par
Con respecto a los valores punteados, éstos corresponden a los valores de penalización en que se tiene el *menor error de clasificación* ($\lambda_{\min}$) y *una desviación después del menor error* ($\lambda_{1DE}$)). Comúnmente, éstos son los candidatos a tomar como parámetros en el modelos final. 

\par
Para éste ejemplo, estos valores se encuentran en la siguiente tabla, y de la gráfica previa, se deduce que a éstos valores de $\lambda$ corresponderán modelos con a lo más $6$ coeficientes *no cero* (más aún, $\lambda_{1DE}$ tendrá un modelo con vector de coeficientes cero).

|          |Error mínimo      |1 Desv. del Error mínimo |
|:---------|:-----------------|:------------------------|
|$\lambda$ |0.695251948018439 |0.919082158397839        |


\par
Para ajustar los modelos y recuperar los vectores de coeficientes se usan los siguientes comandos.

```r
# Se ajustan los modelos para TODOS LOS PARÁMETROS DE VALIDACIÓN CRUZADA.
net <- cv_net$glmnet.fit

# Coeficientes no cero.
coef_net <- net$beta
```

\par
Ahora, se pasa a validar la eficacia del modelo. Para eso, se toman las submuestras de prueba y se calcula el *error de clasificación*. Unicamente se verifican $\lambda_{\min}$ y $\lambda_{1DE}$.

```r
#Error de clasificación
assess.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = familia, 
  s = c( cv_net$lambda.min, cv_net$lambda.1se )
  )$class
```

```
##         1         2 
## 0.4545455 0.4545455 
## attr(,"measure")
## [1] "Misclassification Error"
```

\par
Como puede verse, el modelo es apenas mejor que un volado. Para saber exactamente cuales fueron los errores, se tiene el siguiente comando.

```r
confusion.glmnet( 
  net, 
  newx = R_test, newy = Y_test, 
  family = familia, 
  s = c( cv_net$lambda.min, cv_net$lambda.1se )
  )
```

```
## $`1`
##          True
## Predicted 0 1 Total
##     0     6 5    11
##     Total 6 5    11
## 
##  Percent Correct:  0.5455 
## 
## $`2`
##          True
## Predicted 0 1 Total
##     0     6 5    11
##     Total 6 5    11
## 
##  Percent Correct:  0.5455
```

\par
De todo esto se puede llegara pensar que éste tipo de modelo puede no aportar mucho debido a lo severo que puede penalizar. Sin embargo, tambien se puede deber a que no hay suficiente información para clasificar (no hay una buena proporción de VERDADEROS/FALSOS en éste síntoma, la muestra de prueba es muy pequeña), a que se tuvo una mala semilla ó, como se dijo arriba, al *efecto de porcentajes como variables independientes*. 

\par
Por el momento, se revisará el comportamiento de éste modelo para el resto de los síntomas; sin embargo, antes de pasar a ello, se obtiene la cantidad de especies que aportan información para la predicción del síntoma.



|                     |$\lambda_{\min}$ |$\lambda_{1DE}$ |
|:--------------------|:----------------|:---------------|
|Coeficientes no cero |4                |0               |

Procediendo al resto de los síntomas, los resultados son los que siguen.

### Escurrimiento nasal.

\par
En éste caso, desde la gráfica se puede ver que el modelo será ralo; sin embargo, es curioso que en éste caso el porcentaje de aciertos (equivalentemente, error de clasificación) en la muestra de prueba sea alto, mayor al $70\%$. De cierta manera, el modelo dice que no elegir algo es la manera óptima de clasificar.

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)



|                     |$\lambda_{\min}$  |$\lambda_{1DE}$   |
|:--------------------|:-----------------|:-----------------|
|Valor                |0.702571743434693 |0.702571743434693 |
|Coeficientes no cero |0                 |0                 |

\par
Error de clasificación.

```
##         1         2 
## 0.2727273 0.2727273 
## attr(,"measure")
## [1] "Misclassification Error"
```

\par
Porcentaje de aciertos.

```
## $`1`
##          True
## Predicted 0 1 Total
##     0     8 3    11
##     Total 8 3    11
## 
##  Percent Correct:  0.7273 
## 
## $`2`
##          True
## Predicted 0 1 Total
##     0     8 3    11
##     Total 8 3    11
## 
##  Percent Correct:  0.7273
```

### Presencia de dolor de garganta.

\par
El comportamiento de éste síntoma es símilar al de Fiebre.

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)



|                     |$\lambda_{\min}$  |$\lambda_{1DE}$   |
|:--------------------|:-----------------|:-----------------|
|Valor                |0.700605403469955 |0.805525963548342 |
|Coeficientes no cero |5                 |0                 |

\par
Error de clasificación.

```
##         1         2 
## 0.1818182 0.1818182 
## attr(,"measure")
## [1] "Misclassification Error"
```

\par
Porcentaje de aciertos.

```
## $`1`
##          True
## Predicted 0 1 Total
##     0     9 2    11
##     Total 9 2    11
## 
##  Percent Correct:  0.8182 
## 
## $`2`
##          True
## Predicted 0 1 Total
##     0     9 2    11
##     Total 9 2    11
## 
##  Percent Correct:  0.8182
```

### Tos.

\par
A primera vista, éste síntoma parece tener un comportamiento más regular, pues el modelo trabaja como se esperaba: seleccionando sólo las especies que aportan a la aparición del síntoma, pero sin penalizar severamente al resto. 
Aunque habría que considerar que tanto se puede creer en el error de clasificación, debido a la submuestra de prueba tan pequeña (11 observaciones).

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)


|                     |$\lambda_{\min}$  |$\lambda_{1DE}$   |
|:--------------------|:-----------------|:-----------------|
|Valor                |0.277875287002677 |0.334701766887149 |
|Coeficientes no cero |64                |53                |

\par
Error de clasificación.

```
##         1         2 
## 0.6363636 0.6363636 
## attr(,"measure")
## [1] "Misclassification Error"
```

\par
Porcentaje de aciertos.

```
## $`1`
##          True
## Predicted 0 1 Total
##     0     1 3     4
##     1     4 3     7
##     Total 5 6    11
## 
##  Percent Correct:  0.3636 
## 
## $`2`
##          True
## Predicted 0 1 Total
##     0     1 3     4
##     1     4 3     7
##     Total 5 6    11
## 
##  Percent Correct:  0.3636
```



# Clustering (Réplica de resultados)

\par
Para la resolución de este problema, se hace uso del modelo probabilístico de *Mezclas Dirichlet-Multinomial*. Éste toma en cuenta varias características propias de los datos biológicos con los que se estrá trabajando: Primero, que la información con la que se cuenta (Abundancia) es en esencia discreta, pese a que pueda ser aproximada a variables continuas (Abundancia relativa); la alta diversidad, en comparación con el número de muestras, resulta en datos ralos; tercero, las observaciones varían en profundidad (número de lecturas.)

\par
Para superar estas dificultades, el modelo cambia el enfoque usual de pensar a la muestra como representativa de la comunidad, sino como una realización de ésta. De ahí el uso de herramientas propias de estadística Bayesiana (en la práctica, la aplicación del modelo hace uso de algoritmos de *Esperanza-Maximización* y *MCMC*).

\par
Como su nombre sugiere, el modelo está relacionado con distribuciones Dirichlet y Multinomial. Desde un punto de vista técnico, éstas son la distribución *previa* y *posterior*, respectivamente. En el contexto del problema, la primera se refiere a las metacomunidades a partir de las cuales las se obtienen las observaciones; luego, al tomarse no una sola distribución sino una mezcla de ellas, el modelo captura el hecho de que las observaciones no se generan por una sóla metacomunidad, sino por una mezcla de varias metacomunidades. Con respecto a la parte mutlinomial, ésta es la que se encarga de clasificar cada una de las observaciones en la comunidad correspondiente, de acuerdo a la metacomunidad más probable, obtenida a partir de la mezcla Dirichlet y de la información de la muestra.


## Algoritmo

\par
A continuación, se presenta el algorítmo usado para la replicación de los resultados, junto con el código fuente. A diferencia del artículo original, el lenguaje implementado es R en vez de Mothur.


### Preparación de los datos

\par
Se comienza llamando las librerías necesarias para el modelo: `vegan,  phyloseq, microbiome, DirichletMultinomial`. El resto sirve para la visualización. También se fija la semilla del generador de números aleatorios para replicación.


```r
library(knitr)
library(reshape2)
library(plotly)
#library(dplyr)
#library(magrittr)
#library(ggpubr)
library(vegan)
library(phyloseq)
library(microbiome)
library(DirichletMultinomial)

# Semilla
set.seed(1234)
```

\par
Ahora, se lee la información de los archivos csv y se les asigna da el formato `phyloseq`, necesario para las funciones de la paquetería. Los ingredientes necesarios para esto son la tabla de abundancias totales (`tabla_otu`) y la tabla de nombres de taxa (`tabla_taxa`).


```r
### Ingredientes (¿mínimos?) para objeto phyloseq. ###
# 1.  Tabla de OTUs.
# 2.  Nombres de Taxa.

# OTUs: lee datos -> matriz numérica -> objeto de phyloseq
tabla_otu <- read.csv("oligotype_count.csv", row.names = "X")
otus <- data.matrix(tabla_otu)
OTU <- otu_table(otus, taxa_are_rows = TRUE)

# Taxas: lee datos -> data.frame -> objeto de phyloseq
tabla_taxa <- read.csv("tax_data.csv",row.names = 1)
taxas <- as.matrix(tabla_taxa)
TAX <- tax_table(taxas)

# Crea elemento de la clase phyloseq
physeq <- phyloseq(OTU, TAX)
```

\par
No es necesario tener directamente las abundancias relativas, pues la paquetería permite obtenerlas mediante la transformación `compositional`.
(Lo que realiza esta transformación es normalizar, pasando cada abundancia a su respectivo porcentaje en la muestra.)


```r
#Abundancia relativa.
physeq.comp <- microbiome::transform(physeq, "compositional")
```


\par
El siguiente código se utiliza para identificar la taxa del core. En este caso, se detectan las que están en más de un 0.05\% en toda la muestra. Si se busca filtrar taxa según el porcentaje de prevalencia, se cambia el factor de `prevalence`.


```r
# Sólo taxa del core, prevalente en 0.05% abundancia relativa, 
# en el 50% de las muestras.
# (¿informaciòn descartada? ¿Agregar taxa rara?)
taxa <- core_members(physeq.comp, detection = 0.0005/100, prevalence = 100/100)
#physeq <- prune_taxa(taxa, physeq)
```

\par
Finalmente, se obtiene la matriz de abundancias del core y se le da el formato para pasarlo al modelo. En éste caso conviene tenerla como matriz, pero para el resto convienen elementos de la clase `phyloseq`.

```r
# Conteo de OTU (con la info descartada)
dat <- abundances(physeq)
# Formato muestra x taxa
count <- as.matrix(t(dat))
```


### Ajuste de modelo: Mezclas Dirichlet-Multinomial (DMM)

\par
En la lista `fit` se ajustan el modelo DMM, considerando de 1 a 5 clusters. **La manera correcta de proceder es aplicar el modelo considerando de 1 hasta $n$ clusters, que es lo que se hizo originalmente, pero por tiempos de compilación se deja hasta 5, que ya se vió que es el óptimo**.

```r
# Ajusta modelo.
fit <- lapply(1:5, dmn, count = count, verbose=TRUE)
```



\par
La manera "correcta" de escoger el número de clusters es tomando algún criterio de información. En la gráfica siguiente se toman en cuenta tres: de Akaike (línea rayada), de Bayes (punteada), y de Laplace (sólida). Debido al problema mencionado antes sobre el tiempo de computo no es posible aprecierlo, pero **se escoge como cinco el número de clusters óptimo porque los tres criterios (el paper se fija en Laplace) mostraban que a partir de éste, el aporte que hacía incluir más comunidades no era significante**.

```r
lplc <- sapply(fit, laplace) # Laplace
aic  <- sapply(fit, AIC) # Akaike
bic  <- sapply(fit, BIC) # Bayes
plot(lplc, type="b", xlab="Número de Clusters", ylab="Ajuste de modelo")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30-1.png)

\par
El programa elige el mejor modelo segín el argumento que minimiza el criterio de Laplace; en este caso se sabe de antemano que es cinco, pero con los datos reales **conviene revisarlo personalmente**.


```r
# Escoge mejor modelo de los que calculó.
best <- fit[[which.min(unlist(lplc))]]
```

\par
Una vez elegido el modelo, los parámetros de éste se interpretan como la probabilidad de que una cierta observación haya venido de una comunidad dada (`pi`), y como la variablidad (en el sentido estadístico) que existe dentro de cada grupo (`theta`, un valor menor se interpreta como menor cohesión dentro del grupo). En el ejemplo, se tiene que de las cinco comunidades, las tres primeras son más abundantes (`pi > 0.20`), pero la 1 y la 3 son más homogéneas que el resto (`theta > 120`). Por otro lado, el cluster 5 es el más raro y más heterogéneo (`pi < 14, theta < 25`).

|   |     $\pi$|  $\theta$|
|:--|---------:|---------:|
|1  | 0.2528649| 134.16117|
|2  | 0.2488784|  83.25385|
|3  | 0.2108701| 120.76073|
|4  | 0.1571346|  86.85013|
|5  | 0.1302520|  23.62787|

\par
Con respecto a la *clusterización* en sí, esta se obtiene asignando a cada observación la comunidad que le es más probable de la *mejor* mezcla. Desde aquí también es posible tener una aproximación del tamaño de los clusters sin hacer la clusterización (puede ahorrar tiempo si hay muchos datos).

```r
# Clasificación de cada paciente
gpos <- apply( mixture(best), 1, which.max )
#write.csv(gpos,"grupos_r100.csv")

# Tabla de phylosecuencias clusterizada
sample_data(physeq) <- data.frame(gpos)

# Pacientes por grupo (aproximado)
kable(
  round( best@mixture$Weight * length(tabla_otu) ),
  row.names = TRUE,
  col.names = c("Núm. de observaciones (aprox.)"))
```



|   | Núm. de observaciones (aprox.)|
|:--|------------------------------:|
|1  |                            355|
|2  |                            350|
|3  |                            296|
|4  |                            221|
|5  |                            183|

\par
Explorando la composición de los clusters de manera visual se deducen varias cosas.

  1. Del histograma, la distribción de observaciones en los grupos concuerda casi perfectamente con el aproximado anterior (salvo una clasificación erronea). Con respecto a los reportados en el artículo, los errores no son mayores al orden de 10.

  2. El mapa de calor compara las observaciones (columnas delgadas) con las abundancias promedio de cada grupo (columnas gruesas, resultados de las mezclas Dirichlet en el modelo). Cuando se compara este mapa con los parámetros del modelo y su interpretación se observa que, en efecto, los tres primeros grupos albergan la mayor cantidad de observaciones (hecho tambien claro del histograma), pero además se aprecia que estas observaciones *individualmente* varían en intensidad respecto al promedio. En cambio, el cluster 5 no sólo es el más pequeño (en número de observaciones), sino que individualmente la intensidad de cada observación es muy similar al resto del grupo, lo que coincide con la interpretación del valor $\theta$ en los parámetros del modelo.

  3. Por último, la gráfica de abundancias relativas de cada grupo muestra claramente el aporte promedio de cada especie individual de OTU a cada grupo.


```
## Error in loadNamespace(name): there is no package called 'webshot'
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34-1.png)

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

\par
Ahora, se pasa a revisar sólo la taxa más representativa. Para ello, cada OTU se califica según su *posición promedio* como contribuyente en cada cluster.
Por ejemlo, en la tabla de abajo se muestran los ID de los 10 OTUs que más abundan en cada grupo, ordenados de más abundante a menos. 

```r
# Taxa más prevalecientes en cada grupo
rep <- as.matrix(datos_n[datos_n$grupo==1,"otu"], ncol=1)
for( k in seq(2,5) ){
  rep <- cbind(rep,datos_n[datos_n$grupo==k,"otu"])
}
kable(
  data.frame( rep[1:10,] ), row.names = TRUE
  )
```



|   |   X1|   X2|   X3|   X4|   X5|
|:--|----:|----:|----:|----:|----:|
|1  | 5593| 5593| 5593| 5089| 5593|
|2  | 4359| 5003| 5003| 5904| 5089|
|3  | 5089| 5089| 4575| 4101| 5003|
|4  | 4575| 5034| 4359| 5669| 5090|
|5  | 5003| 4359| 4627| 3085| 5034|
|6  | 4627| 4575| 5034| 1726| 1726|
|7  | 5904| 4101| 5813| 4359| 4359|
|8  | 3085| 1726| 5089| 5034| 2765|
|9  | 5034| 2055| 3322| 5090| 4575|
|10 | 4101| 2765| 4432| 2765| 5904|
Entonces, la calificación que le corresponde al OTU `5089` es
$$
  \frac{3 + 3 + 8 + 1 +2}{5}=3.4.
$$
En general, la calificación de cada OTU está dada por
$$
  \frac{1}{\# \mathrm{Clusters}} 
    \sum_{ i=1 }^{\# \mathrm{Clusters} } P_i( \mathrm{OTU} ), 
$$
donde $P_i( \mathrm{OTU} )$ denota la posición del OTU en el $i$-ésimo cluster. Aplicando ésta calificación a todos los OTU de la muestra, se obtiene que los 20 más representativos son los siguientes. (Una calificación menor equivale a una mayor abundancia a lo largo de toda la muestra.)


|     |taxa                                                                                                                                                                                                                                                                                                                        | calificaciones|
|:----|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------:|
|5593 |Veillonella dispar / Veillonella atypica / Veillonella parvula / Veillonella rogosae                                                                                                                                                                                                                                        |            3.4|
|5089 |Streptococcus sp. / Streptococcus dentisani / Streptococcus mitis / Streptococcus oralis / Streptococcus infantis / Streptococcus tigurinus / Streptococcus lactarius / Streptococcus peroris / Streptococcus pneumoniae                                                                                                    |            3.4|
|4359 |Prevotella melaninogenica / Prevotella scopos / Prevotella sp. / Prevotella histicola / Prevotella veroralis                                                                                                                                                                                                                |            5.0|
|5003 |Streptococcus vestibularis / Streptococcus salivarius / Streptococcus gordonii / Streptococcus sp.                                                                                                                                                                                                                          |            5.2|
|5034 |Streptococcus australis / Streptococcus parasanguinis II / Streptococcus parasanguinis I / Streptococcus sp. / Streptococcus oligofermentans / Streptococcus cristatus / Streptococcus sinensis / Streptococcus sanguinis / Streptococcus gordonii / Streptococcus lactarius / Streptococcus peroris / Streptococcus oralis |            6.4|
|4575 |Prevotella histicola / Prevotella sp. / Prevotella veroralis / Prevotella scopos / Prevotella fusca / Prevotella melaninogenica                                                                                                                                                                                             |           10.0|
|4101 |Haemophilus parainfluenzae / Haemophilus parahaemolyticus / Haemophilus paraphrohaemolyticus / Haemophilus sputorum / Haemophilus sp. / Haemophilus haemolyticus / Haemophilus influenzae                                                                                                                                   |           10.2|
|1726 |Gemella haemolysans / Gemella sanguinis / Gemella morbillorum / Gemella bergeri                                                                                                                                                                                                                                             |           10.4|
|5904 |Neisseria subflava / Neisseria flavescens / Neisseria flava / Neisseria sicca / Neisseria pharyngis / Neisseria mucosa / Neisseria polysaccharea / Neisseria weaveri / Neisseria meningitidis / Neisseria lactamica                                                                                                         |           11.6|
|5090 |Streptococcus peroris / Streptococcus lactarius / Streptococcus sp. / Streptococcus tigurinus / Streptococcus infantis / Streptococcus dentisani / Streptococcus oralis / Streptococcus mitis                                                                                                                               |           12.2|
|2765 |Granulicatella adiacens / Enterococcus italicus / Enterococcus faecalis                                                                                                                                                                                                                                                     |           14.0|
|5398 |Actinomyces sp. / Actinomyces odontolyticus / Actinomyces meyeri / Actinomyces cardiffensis / Actinomyces lingnae [NVP] / Actinomyces georgiae / Actinomyces gerencseriae / Actinomyces massiliensis                                                                                                                        |           14.2|
|2055 |Rothia mucilaginosa                                                                                                                                                                                                                                                                                                         |           15.6|
|2423 |Actinomyces graevenitzii                                                                                                                                                                                                                                                                                                    |           17.2|
|5669 |Veillonella parvula / Veillonella rogosae / Veillonella atypica / Veillonella denticariosi / Veillonella dispar                                                                                                                                                                                                             |           18.6|
|4627 |Prevotella sp. / Prevotella veroralis / Prevotella fusca / Prevotella histicola / Prevotella scopos / Prevotella melaninogenica                                                                                                                                                                                             |           19.6|
|4432 |Prevotella salivae                                                                                                                                                                                                                                                                                                          |           20.0|
|3085 |Fusobacterium periodonticum / Fusobacterium nucleatum subsp. animalis / Fusobacterium sp. / Fusobacterium nucleatum subsp. vincentii / Fusobacterium nucleatum subsp. polymorphum / Fusobacterium naviforme / Fusobacterium nucleatum subsp. nucleatum                                                                      |           20.2|
|5397 |Actinomyces sp. / Actinomyces odontolyticus / Actinomyces meyeri / Actinomyces cardiffensis / Actinomyces lingnae [NVP] / Actinomyces georgiae                                                                                                                                                                              |           24.4|
|3322 |Megasphaera micronuciformis                                                                                                                                                                                                                                                                                                 |           24.8|

\par
Más aún, se puede ver que tan sólo éstos aportan más del $40\%$ de la abundancia de cada cluster. Queda pendiente revisar la composición de cada grupo individualmente y obtener los subconjutos mínimos de OTUs que explican el $50\%$ de cada comunidad, cotejar con el resto de los clusters, y comparar con la lista del artículo.

```
## Error in loadNamespace(name): there is no package called 'webshot'
```


### Diversidad alfa

\par
Los resultados anteriores son alentadores, pues a primera vista parece que replican adecuadamente el proceso de *Clustering*, sin embargo aún queda realizar el análisis estadístico correspondiente. Se comienza revisando la diversidad $\alpha$ dentro de los grupos. Usando el paquete `microbiome`, basta utilizar el comando `alpha`.


```r
### Alfa diversidad. ###
alfa <- microbiome::alpha(physeq, index = "all")
#kable(head(alfa))
```

\par
El comando anterior calcula todos los índices de riqueza, diversidad, rareza, entre otros. Para efectos de este trabajo, sólo se analizan los índices de `shannon` y `chao1`. La siguiente tabla es un ejemplo de tales índices para las primeras observaciones.


|      | Cluster|  Shannon|    Chao1|
|:-----|-------:|--------:|--------:|
|T0001 |       1| 3.475495| 167.2857|
|T0002 |       4| 3.224170| 116.0625|
|T0003 |       2| 2.291909| 112.0000|
|T0004 |       2| 3.304302| 154.5714|
|T0005 |       5| 2.249061|  77.5000|
|T0006 |       1| 3.615769| 165.3333|

\par
Analizando gráficamente el comportamiento del índice de Shannon (gráficas de violín de abajo), se tiene en general los 4 primeros grupos tienen una diversidad de OTU similar. Caso aparte el del grupo 5, que posee una diversidad media más baja. Recordando las deducciones hechas antes, se puede concluir que los primeros grupos son más abundantes, quedan determinados por una mayor cantidad de especies, y por ende, sus elementos tienden a variar más entre sí, dentro del mismo grupo. Por otro lado, el grupo 5 es más raro (menos recurrente), queda determinado por una menor cantidad de especies, e individualmente sus integrantes suelen ser muy parecidos.


```
## Error in loadNamespace(name): there is no package called 'webshot'
```

\par
Para dar formalidad al estudio, se realizan pruebas estadísticas para la hipótesis nula $H_0$: las medianas de los grupos no son distintas; i.e., la $\alpha$-diversidad no es distinta entre grupos. Las pruebas que se utilizan son **Kruskal-Wallis** (pruebas de análisis de varianza *conjunta* entre más de 2 grupos) y **Wilcoxon Rank Sum** (pruabas similares pero a *pares*). Se usan estas pruebas y no ANOVA porque éstas son *no-paramétricas*. También se aprovecha para señalar que en la prueba (Wilcoxon), se ajusta por pruebas multiples usando la corrección de *Benjamini–Hochberg* para asegurar el *FDR* (*False Discovery Rate*), que controla la tasa de descubrimientos falsos; es decir, evita falsos positivos en cuanto a comunidades con diversidad diferente.


```r
# Shannon. SE USA AJUSTE POR FDR.
kruskal.test( shannon ~ gpos, data = m_physeq )
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  shannon by gpos
## Kruskal-Wallis chi-squared = 704.45, df = 4, p-value < 2.2e-16
```

```r
pairwise.wilcox.test( m_physeq$shannon, m_physeq$gpos, p.adjust.method = "BH" ) 
```

```
## 
## 	Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
## 
## data:  m_physeq$shannon and m_physeq$gpos 
## 
##   1       2       3       4      
## 2 < 2e-16 -       -       -      
## 3 < 2e-16 8.7e-09 -       -      
## 4 9.6e-16 < 2e-16 1.7e-06 -      
## 5 < 2e-16 < 2e-16 < 2e-16 < 2e-16
## 
## P value adjustment method: BH
```

\par
Los resultados de las pruebas anteriores son alentadores, pues tanto la prueba en conjunto (Kruskall-Wallis) como las pruebas a pares (Wilcoxon) indican que hay evidencia estadística suficiente para rechazar que todas las comunidades, o algún par de ellas, tengan la misma distribución, respectivamente. Más aún, los resultados son congruentes con lo reportado por Lee et. al.

\par
Realizando un procedimiento análogo (sólo se muestran resultados) para el índice de Chao1, se obtiene que *en conjunto* las 5 comunidades sí tienen una $\alpha$-diversidad distinta; sin embargo, en las pruebas *a pares*, no se encontró evidencia estadística suficiente para separar al grupo 3 del 4. En caso de encontrar éste hecho con los datos reales, se recomienda analizar individualmente las comunidades que no pasen esta prueba.


```
## Error in loadNamespace(name): there is no package called 'webshot'
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  chao1 by gpos
## Kruskal-Wallis chi-squared = 622.43, df = 4, p-value < 2.2e-16
```

```
## 
## 	Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
## 
## data:  m_physeq$chao1 and m_physeq$gpos 
## 
##   1       2       3       4      
## 2 < 2e-16 -       -       -      
## 3 < 2e-16 6.7e-14 -       -      
## 4 < 2e-16 1.4e-08 0.32    -      
## 5 < 2e-16 < 2e-16 < 2e-16 < 2e-16
## 
## P value adjustment method: BH
```

### Diversidad beta

\par
Por último, la $\beta$-diversidad (diversidad entre grupos) también se explora gráficamente primero, y luego se pasa a las pruebas estadísticas. Las distancias que se consideran son *Bray-Curtis* y *Jaccard*, siguiendo los resultados del paper. En la práctica, se utiliza la función `divergence` de la paquetería `microbiome`, que ya tiene incluidos dichos métodos.

```r
# Listas donde se guardan las distancias.
clusters <- as.list(1:5)
bray <- as.list(1:5)
jaccard <- as.list(1:5)

# Cálculo de distancias por Cluster.
for (k in seq(1:5)) {
  
  # Recupera individualmente cada cluster.
  eval(
    parse(
      text=sprintf(
        'clusters[[%s]] <- subset_samples(physeq, sample_data(physeq)$gpos == %s)', 
        k, 
        k
        )
      )
    )
  
  # Referencia a la mediana del cluster.
  ref <- apply(abundances(clusters[[k]]), 1, median)
  
  # Bray-Curtis.
  bray[[k]] <- array( as.numeric( 
    unlist( divergence(clusters[[k]], ref, method = "bray") )
    ) )
  
  # Jaccard.
  jaccard[[k]] <- array( as.numeric( 
    unlist( divergence(clusters[[k]], ref, method = "jaccard") )
  ) )
}
```

\par
Ya calculados los índices de disimilaridad, se procede a analizar cómo se distribuyen. Comenzando con el índice de Bray-Curtis, de la gráfica de cajas parece ver que la variación entre los grupos no es la misma, lo cuál indicaría que la clusterización fue adecuada.

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

\par
Para probar ahora la hipótesis $H_0$ de que no hay diferencias entre 2 o más vectores de medias (las abundancias promedio de los vectores de especies según cada grupo), i.e., para ver que hay una diferencia significativa entre grupos, se pasa a realizar una prueba PERMANOVA. Ésta, como su nombre sugiere, tiene que ver con las pruebas de análisis de varianza entre grupos: la parte de MANOVA se refiere a *Multivariate ANOVA*, que es el análisis de varianza para variables respuesta multivariadas o vectores; la parte de PER se refiere a pruebas basadas en permutación de los datos. Éstas son pruebas *no paramétricas*, que se basan en permutar (aleatoriamente) submuestras de los datos para obtener una aproximación de la función de distribución empírica y trabajar con ella, en vez de asumir supuestos distribucionales.

\par
Aplicando la prueba de PERMANOVA a las distancias Bray-Curtis, se obtiene un coeficiente de correclación $R^2 = 0.203$, es decir, $20.3\%$ de la variación total se explica por los clusters (un $0.4\%$ menor al reportado en el artículo de influenza). El $p$-valor correspondiente es menor a $0.001$, así que considerando el nivel de significancia usual de $0.05$, se concluye que hay evidencia significativa para rechazar la hipótesis $H_0$; dicho de otra manera, la diferencia entre grupos es significativa.

```r
bc_dist <- phyloseq::distance(physeq, method = "bray")
met <- data.frame(sample_data(physeq))
adonis( bc_dist ~ factor(gpos), data = met, method = "bray")
```

```
## 
## Call:
## adonis(formula = bc_dist ~ factor(gpos), data = met, method = "bray") 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## factor(gpos)    4     58.54 14.6350  89.151 0.20301  0.001 ***
## Residuals    1400    229.82  0.1642         0.79699           
## Total        1404    288.36                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

\par
Realizando un procedimiento análogo con el índice de Jaccard, se obtienen resultados similares (con $R^2 = 0.1412$).

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

```
## 
## Call:
## adonis(formula = jac_dist ~ factor(gpos), data = met, method = "jaccard") 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## factor(gpos)    4     58.35 14.5879  57.553 0.14122  0.001 ***
## Residuals    1400    354.86  0.2535         0.85878           
## Total        1404    413.21                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# Conclusiones y trabajo futuro.

\par
Del modelo original propuesto, lo primero que hay que hacer es hacer las modificaciones necesarias para que sea posible introducir las abundancias relativas (proporciones, variables aleatorias acotadas) como variables predictoras en el contexto de regresión logística. Revisando la literatura se encontró que esta modificación existe (*Fractional logit model*), pero aún no se ha encontrado paquetería que ligue dicho modelo a regresión penalizada. Se están considerando como opciones programar el modelo, o modificar el propuesto de manera que se haga uso de los resultados de Clustering para la parte de replicación.

\par
Con respecto a la parte de replicación, parece ser que la parte del problema de detección de comunidades quedó cubierto. Queda cotejar con los datos del artículo que se presentó el pasado lunes (29 de junio) para ver si en efecto se tratan de las mismas comunidades que presentaron los autores en el su otro paper. También queda revisar la teoría de *modelos de falla acelerada*, que son con los que Lee et al. correlacionaron la sintomatología viral de la Influenza con la microbiota respiratoria.

# Anexos

\par
Los anexos se pueden hallar en la siguiente liga: https://github.com/nselem/mg-covid/blob/master/joshue/Joshu%C3%A9%20Ricalde%20-%20Reporte%202%20(Anexos).pdf


# Referencias

1. Lee, K. H., Foxman, B., Kuan, G., López, R., Shedden, K., Ng, S., ... & Gordon, A. (2019). The respiratory microbiota: associations with influenza symptomatology and viral shedding. Annals of epidemiology, 37, 51-56.
2. Meier, L., Van De Geer, S., & Bühlmann, P. (2008). The group lasso for logistic regression. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(1), 53-71.
3. Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
4. Efron, B., & Hastie, T. (2016). Computer age statistical inference (Vol. 5). Cambridge University Press.
5. Giraud, C. (2014). Introduction to high-dimensional statistics (Vol. 138). CRC Press.
6. Algamal, Z. Y., & Lee, M. H. (2015). Penalized logistic regression with the adaptive LASSO for gene selection in high-dimensional cancer classification. Expert Systems with Applications, 42(23), 9326-9332.
7. Holmes, I., Harris, K., & Quince, C. (2012). Dirichlet multinomial mixtures: generative models for microbial metagenomics. PloS one, 7(2), e30126.
8. Xia, Y., Sun, J., & Chen, D. G. (2018). Statistical analysis of microbiome data with R. Singapore:: Springer.
