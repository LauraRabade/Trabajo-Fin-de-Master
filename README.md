# Trabajo Fin de Master

En este repositorio aparece todo el código generado durante el Trabajo Fin de Máster y los archivos se encuentran divididos en dos partes; una carpeta que contiene los archivos generados para R y otra los archivos generados para Bash. El uso de los archivos para los diferentes análisis durante el desarrollo del trabajo ha sido el siguiente:

1. Analisis de metilación del ADN
- Estimación de la proporción de células epiteliales: EpithelialCells.R
- Análisis con el software de ChAMP: ChAMP_Analysis.R
- Función para el uso de DMRCate en el pipeline de ChAMP: champDMR_funct.R
- Filtrado de las DMPs para cada comparación y obtención de gráficas con las características de las sondas: DMPs_filtering.R

2. Variación común asociada a las DMRs
- Función para obtener las listas de gene sets correspondientes a las DMRs: DMR_DataCleaning.R
- Preparación de los summary statistics de cada GWAS: script_summary_prep.sh
- Análisis de los gene sets con MAGMA: script_magma.sh 

3. Uso del reloj epignético EPIC clock
- Estimación de la edad gestacional y de la aceleración de los pacientes, así como el cálculo de estadísticos t-Student: GestationalDNAmAge.R
