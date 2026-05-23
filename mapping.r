##Se llama a las librerías correspondientes##
library(GenomicRanges)
library(GenomicAlignments)
library(ggplot2)
library(stringr)
##Se lee el archivo clinvar_simp.txt con encabezados
clinvar <- read.table('./clinvar_simp.txt', header = T, sep = "\t")

#Lista que contiene los nombres de todos los cromosomas en formato de NCBI RefSeq
ChromosomeNCBI <- c(
  "1" = "NC_000001.11",
  "2" = "NC_000002.12",
  "3" = "NC_000003.12",
  "4" = "NC_000004.12",
  "5" = "NC_000005.10",
  "6" = "NC_000006.12",
  "7" = "NC_000007.14",
  "8" = "NC_000008.11",
  "9" = "NC_000009.12",
  "10" = "NC_000010.11",
  "11" = "NC_000011.10",
  "12" = "NC_000012.12",
  "13" = "NC_000013.11",
  "14" = "NC_000014.9",
  "15" = "NC_000015.10",
  "16" = "NC_000016.10",
  "17" = "NC_000017.11",
  "18" = "NC_000018.10",
  "19" = "NC_000019.10",
  "20" = "NC_000020.11",
  "21" = "NC_000021.9",
  "22" = "NC_000022.11",
  "X" = "NC_000023.11",
  "Y" = "NC_000024.10",
  "MT" = "NC_012920.1"
)

clinvar$Chromosome <- ChromosomeNCBI[as.character(clinvar$Chromosome)]
#Creación un objeto de rangos genómicos a partir de las variaciones encontradas en OMIM y clinvar
mendeliandsGR <- GRanges(
  seqnames = clinvar$Chromosome,
  ranges = IRanges(
    start = as.numeric(clinvar$Start),
    end = as.numeric(clinvar$End)
  )
)

#Lectura de la anotación del genoma de referencia (puede editarse la dirección del archivo)
reference <- read.table('./GRCh38_latest_genomic.gff', sep = "\t")
#Solo se utilizan las primeras 7 líneas correspondientes al cromosoma, [], biotipo, inicio, fin, . y cadena
reference <- unique(reference[, 1:7])
#Los biotipos correspondientes a "region" y "match" son redundantes, pues representan a los cromosomas o regiones genómicas y estos biotipos no nos interesan
reference <- reference[!(reference$V3 %in% c("region", "match")), ]
#Para simplificar la transformación de la anotación del genoma de referencia se nombran cada una de sus columnas
colnames(reference) <- c(
  "Seqnames",
  "Source",
  "type",
  "Start",
  "End",
  ".",
  "Strand"
)
#creación del objeto Granges con el genoma de referencia
genome_ref <- as(reference, "GRanges")

#Conteo de los sobrelapes para saber cuantas variantes mapearon y cuantas no
sum(overlapsAny(mendeliandsGR, genome_ref))

sum(!overlapsAny(mendeliandsGR, genome_ref))

#Con findOverlaps Obtenemos los índices de los sobrelapes en una matriz
overlaps <- as.matrix(findOverlaps(mendeliandsGR, genome_ref))


#Con la matríz de sobrelapes que tiene los índices podemos buscarlos en el genoma de referencia y colocarlos en la lista de biotipos a los que mapea cada variante (ya que pueden ser varios se crea una lista)
temp <- unique(overlaps[, 1])
biotype <- vector(mode = "list", length = length(temp))
for (i in temp) {
  biotype[[i]] <- as.vector(unique(
    genome_ref[overlaps[overlaps[, 1] == i, 2], ]$type
  ))
}
#Lista temporal para guardar los biotipos ya que la lista será modificada mas adelante
biotypos2 <- biotype


#Debido a que no existe el biotipo de "intrón" consideramos a todo lo que contenga el biotipo "gene", "mRNA" y "transcript" como un intrón si no contiene alguna otra notación (Como CDC o exón)
categories <- c("gene", "mRNA", "transcript")
for (i in temp) {
  if (length(biotypos2[[i]]) == 1) {
    if (biotypos2[[i]] == "CDS") {
      biotypos2[i] <- "exon"
    }
  }
  if (all(biotypos2[[i]] %in% categories)) {
    biotypos2[[i]] <- c(biotypos2[[i]], "intron")
  }
}
###si existe la notacion exon solo se mantiene esa notación, si existe CDS se reemplaza por exon y solo se mantiene esa notación, si existe la notación intron se mantiene solo esa y si existe alguna notación que tenga "gene","mRNA" y "transcript" pero no exon ni CDS y además tiene otra anotación se cambia por intron y solo se mantiene esa anotación.
for (i in temp) {
  if ("exon" %in% biotypos2[[i]]) {
    biotypos2[[i]] <- "exon"
  }
  if ("CDS" %in% biotypos2[[i]]) {
    biotypos2[[i]] <- "exon"
  }
  if ("intron" %in% biotypos2[[i]]) {
    biotypos2[[i]] <- "intron"
  }
  if (any(grepl(categories, biotypos2[[i]]))) {
    biotypos2[[i]] <- "intron"
  }
}


#No todas las notaciones contienen biotipo, por lo cual pueden quedar algunos Null que son eliminados
biotypos2 <- biotypos2[!(sapply(biotypos2, is.null))]

#se calcula el numero de variantes que tienen exon, intron y ninguna de las dos (intergenic)
percent_cod <- c(
  length(grep('intron', biotypos2)),
  length(grep('exon', biotypos2)),
  length(biotypos2) -
    length(grep('intron', biotypos2)) -
    length(grep('exon', biotypos2))
)
names(percent_cod) <- c('Intron', 'Exon', 'Intergenic')
#con prop.table se calculan las proporciones de la tabla y se multiplica por 100 para tener la frecuencia relativa.
percent_cod <- as.data.frame(prop.table((as.matrix(percent_cod))) * 100)
#se le agrega el nombre freq a la columna para su uso en ggplot
colnames(percent_cod) <- 'freq'


percent_plot_intergenic <- ggplot(
  percent_cod,
  aes(x = "", y = freq, fill = (rownames(percent_cod)))
) +
  geom_bar(width = 1, stat = "identity", position = "stack") +
  coord_polar("y") +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Porcentaje de mutaciones en regiones 
       genicas e intergenicas"
  ) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  theme(
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 16)
  )

percent_plot_intergenic
