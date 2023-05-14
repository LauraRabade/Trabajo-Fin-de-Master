# Trabajo Fin de Master

En este repositorio aparece todo el código generado durante el Trabajo Fin de Máster y los archivos se encuentran divididos en dos partes; una carpeta que contiene los archivos generados para R y otra los archivos generados para Bash. El uso de los archivos para los diferentes análisis durante el desarrollo del trabajo ha sido el siguiente:

# Analisis de metilación del ADN
## 1. Estimación de la proporción de células epiteliales: EpithelialCells.R
## 2. Análisis con el software de ChAMP: ChAMP_Analysis.R
## 3. Función para el uso de DMRCate en el pipeline de ChAMP: champDMR_funct.R
## 4. Filtrado de las DMPs para cada comparación y obtención de gráficas con las características de las sondas: DMPs_filtering.R

# Variación común asociada a las DMRs
## 1. Función para obtener las listas de gene sets correspondientes a las DMRs: DMR_DataCleaning.R
## 2. Preparación de los summary statistics de cada GWAS: script_summary_prep.sh
## 3. Análisis de los gene sets con MAGMA: script_magma.sh 

# Uso del reloj epignético EPIC clock
## 1. Estimación de la edad gestacional y de la aceleración de los pacientes, así como el cálculo de estadísticos t-Student: GestationalDNAmAge.R
