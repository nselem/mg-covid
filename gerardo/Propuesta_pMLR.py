# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 22:36:26 2020

@author: Gerardo
"""

"""PROPUESTA: REGRESIÓN LOGÍSTICA MULTINOMIAL PENALIZADA"""
import numpy as np
import math
import scipy
import copy
import matplotlib.pyplot as plt
from scipy.stats import norm, probplot, uniform, gamma, kstest
import seaborn
import csv
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing

#DATOS DE MICROBIOMAS
datos=np.loadtxt(open(r"C:\Users\gerar\Documents\oligotype_count.csv", "rb"),dtype=str, delimiter=',',skiprows=1)
aux=[[float(datos[i][j]) for j in range(1,len(datos[i])) ]for i in range(len(datos))]
mbiomas=np.transpose(aux)

#NORMALIZACIÓN DE LOS DATOS
mbiomasN=[[] for i in range(len(mbiomas))]
for i in range(len(mbiomasN)):
    mbiomasN[i]=mbiomas[i]/np.std(mbiomas[i])
    
#DATOS DE GRAVEDAD DE LA ENFERMEDAD: DURACIÓN DE LA FIEBRE
#Índices de submuestra de estudio
arch=np.loadtxt(open(r"C:\Users\gerar\Documents\study_data.csv", "rb"),dtype=str, delimiter=',',skiprows=1)
index=[int(arch[i][2]) for i in range(len(arch))]
#Elección de submuestra
mbiomasNS=[mbiomasN[index[i]-1] for i in range(len(index))]
#Lectura de datos de duración de fiebre
arch1=np.loadtxt(open(r"C:\Users\gerar\Documents\fever_dur.csv", "rb"),dtype=str, delimiter=',',skiprows=1)
index1=[[int(arch1[i][4]), int(arch1[i][5]) ]for i in range(len(arch1))]
#Gravedad de enfermedad para la submuestra
gravedad=[0 for i in range(len(mbiomasNS))]
for i in range(len(index)):
    for j in range(len(index1)):
        if index[i]==index1[j][0]:
            gravedad[i]=index1[j][1]
            break
    
#VISUALIZACIÓN RALEZA DE LOS DATOS
#Mapa de calor
plt.figure()
plt.imshow(np.transpose(mbiomasN),cmap='hot',interpolation='nearest')
plt.xlabel('Muestras')
plt.ylabel('Taxas')
plt.title("Matriz de datos")
plt.show()
#Otro mapa de calor
plt.figure()
fig=seaborn.heatmap(np.transpose(mbiomasN),vmin=np.min(mbiomasN),vmax=np.max(mbiomasN))
plt.xlabel('Muestras')
plt.ylabel('Taxas')
plt.title("Matriz de datos")
plt.show()

#MATRIZ DE COVARIANZAS Y CORRELACIÓN
#Cálculo matriz de covarianzas
cov=np.cov(np.transpose(mbiomasN))
plt.figure()
fig=seaborn.heatmap(cov,vmin=np.min(cov),vmax=np.max(cov))
plt.xlabel('Taxas')
plt.ylabel('Taxas')
plt.title("Matriz de covarianzas")
plt.show()
#Un acercamiento covarianzas
zoom=50
cov1=[[cov[i][j] for j in range(zoom)] for i in range(zoom)]
plt.figure()
fig=seaborn.heatmap(cov1,vmin=np.min(cov1),vmax=np.max(cov1))
plt.xlabel('Taxas')
plt.ylabel('Taxas')
plt.title("Matriz de covarianzas")
plt.show()
#Cálculo de matriz de correlaciones
corr=np.corrcoef(np.transpose(mbiomasN))
plt.figure()
fig=seaborn.heatmap(corr,vmin=np.min(corr),vmax=np.max(corr))
plt.xlabel('Taxas')
plt.ylabel('Taxas')
plt.title("Matriz de correlaciones")
plt.show()
#Un acercamiento correlaciones
zoom=50
corr1=[[corr[i][j] for j in range(zoom)] for i in range(zoom)]
plt.figure()
fig=seaborn.heatmap(corr1,vmin=np.min(corr),vmax=np.max(corr))
plt.xlabel('Taxas')
plt.ylabel('Taxas')
plt.title("Matriz de correlaciones")
plt.show()

#REGRESIÓN LOGÍSTICA MULTINOMIAL PENALIZADA
lr=LogisticRegression(max_iter=1000,penalty='l2',C=1.0).fit(mbiomasNS,gravedad)
#Gráfica
plt.figure()
for i in range(len(lr.coef_)):
    plt.bar([k for k in range(len(lr.coef_[i]))], lr.coef_[i],align="center")
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial")
plt.show()
#Primer coeficiente (sanos)
plt.figure()
plt.bar([k for k in range(len(lr.coef_[0]))], lr.coef_[0],align="center")
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial (sanos)")
plt.show()
#Tercer coeficiente (graves)
plt.figure()
plt.bar([k for k in range(len(lr.coef_[2]))], lr.coef_[2],align="center")
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial (graves)")
plt.show()
#Ejemplo con umbrales
#Para el primer coeficiente (sanos)
plt.figure()
plt.bar([k for k in range(len(lr.coef_[0]))], lr.coef_[0],align="center")
plt.hlines(y=0.05,xmin=1,xmax=230)
plt.hlines(y=-0.05,xmin=1,xmax=230)
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial (sanos)")
plt.show()
#Para el tercer coeficiente (graves)
plt.figure()
plt.bar([k for k in range(len(lr.coef_[2]))], lr.coef_[2],align="center")
plt.hlines(y=0.05,xmin=1,xmax=230)
plt.hlines(y=-0.05,xmin=1,xmax=230)
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial (graves)")
plt.show()
#Para el segundo coeficiente (leves)
plt.figure()
plt.bar([k for k in range(len(lr.coef_[1]))], lr.coef_[1],align="center")
plt.hlines(y=0.05,xmin=1,xmax=230)
plt.hlines(y=-0.05,xmin=1,xmax=230)
plt.xlabel('Taxas')
plt.ylabel('Coeficientes estimados')
plt.title("Regresión logística multinomial (leves)")
plt.show()
#Variables con coeficientes más grandes
#Sanos
print("Coeficientes sanos:")
for i in range(len(lr.coef_[0])):
    if lr.coef_[0][i]>0.05:
        print(i+2)
#Graves
print("Coeficientes graves:")
for i in range(len(lr.coef_[2])):
    if lr.coef_[2][i]>0.05:
        print(i+2)
#Leves
print("Coeficientes leves:")
for i in range(len(lr.coef_[1])):
    if lr.coef_[1][i]>0.05:
        print(i+2)


