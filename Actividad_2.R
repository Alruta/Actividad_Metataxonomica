setwd("/media/alvaro/DISCO_DURO_ART/Master_Bioinformatica/Poblaciones/Actividad_2")
library(qiime2R)
library(vegan)
library(microbiome)  
library(phyloseq)

# Carga de datos y creación del objeto phyloseq
SVs<-read_qza("dada2_table.qza") # Tabla de dada2 con conteos (feature table)
ASV <- otu_table(SVs$data, taxa_are_rows = TRUE) 

taxonomy <- read_qza("taxonomy.qza")
taxonomy <- parse_taxonomy(taxonomy$data) 
write.table(taxonomy, file='taxonomy.tsv', quote=FALSE, sep='\t', col.names = TRUE)
taxmat <- as.matrix(read.table("taxonomy.tsv", sep="\t"))
TAX = tax_table(taxmat)

metadata <- read_q2metadata("metadata.tsv")
rownames(metadata) <- metadata$SampleID
metadata <- sample_data(metadata)

tree<-read_qza("rooted-tree.qza")

physeq <- phyloseq(ASV, TAX, metadata, tree$data)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))

#  1. Haz un subset de la orden Clostridiales, ¿Cuántos taxones había y cuántos taxones hay ahora?
physeq.clos = subset_taxa(physeq, Order == "Clostridiales")

length(taxa_names(physeq)) # Taxones originales
length(taxa_names(physeq.clos)) # Taxones tras hacer el subset


# 2. Representa para cada muestra su abundancia para cada una de las familias.
plot_bar(physeq.clos, fill = "Family", title = "Abundancia de familias del orden Clostridiales") # Familias en cada muestra del orden clostridiales
plot_bar(physeq, fill = "Family", title = "Abundancias de familias") # Familias en cada muestra de todos los taxones 


# 3. Calcula el valor medio de diversidad alfa (Índice de Shannon) en las muestras según la variable donor_status, 
# que diferencia a las muestras según tengan un trasplante de microbioma de personas sanas o con la enfermedad de 
# Parkinson, haciendo uso del paquete microbiome. ¿ Hay diferencias significativas entre los grupos? 
# ¿Qué grupo de muestras tiene una mayor diversidad alfa según este índice? Incluye el código utilizado para 
# dar tu respuesta.

d <- meta(physeq)
d$diversity <- microbiome::diversity(physeq, "shannon")$shannon
d$diversity # Calculo del indice de shannon para cada muestra e introducción en una variable con metadatos

donor_status <- levels(d$donor_status)

# Combinaciones posibles de comparaciones
donor.pairs <- combn(seq_along(donor_status), 2, simplify = FALSE, FUN = function(i)donor_status[i])


# Represtacion de la diversidad respecto a distintas variables
ggplot(d, aes(x= donor_status, y=diversity, fill=donor_status))+geom_violin(trim=FALSE)

# Test para ver si hay diferencias entre los grupos 
spl <- split(d$diversity, d$donor_status)
mean(spl$Healthy) # Media del indice de shannon en ratones con implantes sanos
mean(spl$PD) # Media del indice de shannon en ratones con impantes de enfermos con Parkinson
ks.test(spl$Healthy, spl$PD)$p.value


# 4. Utiliza el paquete microbiome para obtener los índices de diversidad alfa, ¿Cuál es la muestra con menor 
# índice de Shannon y qué valor tiene?

alpha <- microbiome::alpha
tab <- alpha(physeq, index = "all")
colnames(tab) # Todos los distintos indices calculados

rownames(tab[which.max(tab$diversity_shannon),]) # Muestra con mayor indice de shannon
tab[which.max(tab$diversity_shannon),"diversity_shannon"] # y el valor que tiene
d[rownames(tab[which.max(tab$diversity_shannon),]),] # Metadados de la muestra

rownames(tab[which.min(tab$diversity_shannon),]) # Muestra con menor indice de shannon
tab[which.min(tab$diversity_shannon),"diversity_shannon"] # y el valor que tiene
d[rownames(tab[which.min(tab$diversity_shannon),]),]


# 5. Haz un subset de la clase Betaproteobacteria y representa en un árbol filogenético cada 
# uno de los géneros contenidos en esta clase, distinguiéndolos con un marcador que indique la 
# familia en el árbol y coloreando las muestras según la variable “genotype_and_donor_status”. 
# Adjunta el código que utilizas para crearlo, así como el gráfico generado y qué conclusiones 
# extraes de la observación de este árbol.  

physeq.beta = subset_taxa(physeq, Class == "Betaproteobacteria")

plot_tree(physeq.beta, color = "genotype_and_donor_status", shape = "Family", 
          label.tips = "Genus", size = "abundance", plot.margin = 0.5, ladderize = TRUE)






