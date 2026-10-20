# Referencia espacial bayesiana sin exposición


## Metodología


Se ajusta el ln(1+conteo de imagen) con anatomía, bases espaciales anisotrópicas y efectos propios de corte con distribuciones previas propias. La actualización conjugada Normal–Gamma inversa produce intervalos t de Student distintos para media latente y nueva célula. Se fijan genes, longitudes y distribuciones previas; se omiten bloques completos de cortes y vecinos antes de aprender escalas, centros o coeficientes. No se utilizan etiquetas celulares derivadas de expresión.


## Resultados


HY: se evaluaron las 37,811 células reales de 20 cortes físicos y 10 genes medidos fijos. Mediana del RMSE por gen = 0.140; mediana de cobertura predictiva 95% = 97.3%, frente a 0.0% entre células positivas. La mediana de probabilidad predictiva gaussiana bajo cero es 41.1%. El RMSE espacial mejora la referencia anatómica en 1/10 genes. Se conservan estos resultados aunque la calibración sea deficiente.

ACA: se evaluaron las 14,582 células reales de 18 cortes físicos y 10 genes medidos fijos. Mediana del RMSE por gen = 0.373; mediana de cobertura predictiva 95% = 94.0%, frente a 66.5% entre células positivas. La mediana de probabilidad predictiva gaussiana bajo cero es 38.2%. El RMSE espacial mejora la referencia anatómica en 4/10 genes. Se conservan estos resultados aunque la calibración sea deficiente.


## Interpretación y límites


La verosimilitud gaussiana es una referencia insuficiente para predecir células positivas en estos datos, no una reconstrucción validada ni un modelo generador de conteos. Las medias ajustadas y sus intervalos latentes también pueden ser negativos. La masa predictiva imposible bajo cero, la escasa mejora del RMSE y la recuperación deficiente de la cola positiva limitan conjuntamente su utilidad, aunque la cobertura global parezca próxima al 95%. La inclusión del intervalo entre positivos es un diagnóstico condicionado al resultado; no exige cobertura nominal 95%en ese subconjunto. Los intervalos creíbles son condicionales a bases y distribuciones previas; no describen incertidumbre entre animales ni precisión de registro. Omitir cortes evalúa transporte espacial dentro del donante; no valida una respuesta PM. Las matrices ajustadas a toda la referencia se separan de predicciones y métricas omitidas. No se generan ni desplazan células ni se les atribuye una respuesta de exposición medida.


## Fuentes


[Atlas original](https://doi.org/10.1038/s41586-023-06808-9); [análisis gaussiano conjugado](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf); [implementación t de Student](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).
