# Proyecto Final — Genómica Humana  
## Análisis de la distribución genómica de mutaciones asociadas a enfermedades mendelianas

---

# Objetivo general del proyecto

Determinar qué porcentaje de las mutaciones asociadas a enfermedades mendelianas se localizan en regiones codificantes, intrónicas e intergénicas del genoma humano, y posteriormente evaluar el posible impacto funcional de variantes representativas mediante herramientas bioinformáticas.

---

# Flujo de trabajo general del proyecto

El flujo de trabajo completo del proyecto se divide en dos fases principales:

## Fase I. Análisis global de distribución genómica

Esta fase responde la pregunta principal del proyecto:

> ¿Qué porcentaje de las mutaciones asociadas a enfermedades mendelianas se encuentran en regiones codificantes, intrones y regiones intergénicas?

### Pipeline general

```text
ClinVar (dataset completo)
        ↓
Filtrado de variantes relevantes
        ↓
Extracción de coordenadas genómicas
        ↓
Mapeo contra anotación del genoma humano (GRCh38)
        ↓
Clasificación de variantes:
    - Exónicas (codificantes)
    - Intrónicas
    - Intergénicas
        ↓
Cálculo de porcentajes
        ↓
Visualización de resultados

# Fase II 
# Análisis funcional de variantes representativas

Selección de variantes representativas
        ↓
Análisis funcional según categoría:
    - variantes codificantes
    - variantes intrónicas
    - variantes regulatorias/intergénicas
        ↓
Predicción de efecto molecular
        ↓
Propuesta de mecanismo biológico asociado a enfermedad

# 1. Obtención del dataset de variantes clínicas

Se descargó el dataset completo de variantes clínicas desde ClinVar (NCBI), utilizando la versión basada en el ensamblaje de referencia humano GRCh38.

Archivo descargado:

clinvar.vcf.gz

Posteriormente, el archivo fue descomprimido para su procesamiento:

clinvar.vcf

Formato:

VCF (Variant Call Format)

Este archivo contiene información sobre:

cromosoma
posición genómica
alelo de referencia
alelo alternativo
significancia clínica
anotaciones clínicas
asociaciones con enfermedades

# 2. Filtrado de variantes de interés

Se implementó un script en Python (clean_clinvar.py) para filtrar el dataset completo según criterios biológicamente relevantes.

Criterios de filtrado

Se conservaron únicamente variantes que cumplieran:

asociación con enfermedades registradas en OMIM
clasificación clínica Pathogenic
variantes tipo SNV (Single Nucleotide Variant)
Justificación
OMIM

OMIM se utilizó como proxy para identificar enfermedades mendelianas.

Pathogenic

Se seleccionaron variantes con evidencia clínica fuerte de patogenicidad.

SNVs

Se restringió el análisis a variantes puntuales para evitar ambigüedad en la asignación de coordenadas genómicas de variantes estructurales.

# 3. Resultados del filtrado
Dataset inicial
Total variantes en ClinVar: 4,434,137
Variantes asociadas a OMIM
1,438,127
Variantes patogénicas
159,853
SNVs finales
85,237
# 4. Generación de archivo simplificado

A partir del dataset filtrado se generó:

clinvar_simp.txt

Columnas conservadas:

Name
Gene
Chromosome
Start
End

Este archivo fue utilizado como input para el análisis genómico en R.

# 5. Análisis de localización genómica

Se utilizó un script en R (mapping.r) para mapear las variantes filtradas contra la anotación del genoma humano.

Archivos utilizados
Dataset de variantes
clinvar_simp.txt
Anotación del genoma humano
GRCh38_latest_genomic.gff

Formato:

GFF (General Feature Format)

Contiene:

genes
exones
CDS
transcritos
mRNA
otras anotaciones genómicas
# 6. Metodología de mapeo
Conversión a coordenadas genómicas

Las variantes se transformaron a objetos GRanges.

Esto permitió representar cada variante como una coordenada genómica formal.

Intersección con anotación del genoma

Se utilizaron funciones de Bioconductor:

GenomicRanges
GenomicAlignments

Particularmente:

findOverlaps()

para determinar en qué elementos genómicos cae cada variante.

## Clasificación de variantes

Las variantes se clasificaron como:

Exónicas

Si intersectaban con:

exon
CDS
Intrónicas

Inferidas cuando la variante intersectaba con:

gene
transcript
mRNA

pero no con regiones exónicas.

Intergénicas

Variantes que no se clasificaron en las categorías anteriores.

## Resultados obtenidos
Distribución genómica de variantes mendelianas
Categoría	Porcentaje
Exónicas	82.75%
Intrónicas	16.50%
Intergénicas	0.75%
Interpretación preliminar
Regiones exónicas

La mayoría de las variantes mendelianas reportadas en ClinVar se encuentran en regiones codificantes.

Esto es consistente con:

mayor detectabilidad en estudios clínicos
predominio histórico de secuenciación de exoma
mayor facilidad para interpretar consecuencias proteicas
Regiones intrónicas

Una fracción importante corresponde a variantes intrónicas.

Estas podrían incluir:

variantes de splicing
alteraciones en sitios donador/aceptor
cambios regulatorios intragénicos
Regiones intergénicas

Representan una proporción muy baja.

Esto probablemente refleja sesgo de anotación clínica, más que ausencia de relevancia funcional.

Visualización generada

Se generó un gráfico circular mostrando la distribución porcentual de variantes.

Archivo pendiente de exportación:

mutaciones_mendelianas_piechart.png