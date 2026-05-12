################## RNA-Seq Veri Analizi E??itimi ###########################

##################### R ??le Paket ??ndirme ve ??a????rma #####################

# Paketleri ??ndirme

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Sc.sgd.db")

#  org.Sc.sgd.db paketi R'da olmad?????? i??in Bioconductor'dan al??nd??.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# clusterProfiler paketi R'da olmad?????? i??in Bioconductor'dan al??nd??.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

# pathview paketi R'da olmad?????? i??in Bioconductor'dan al??nd??.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

# EnhancedVolcano paketi R'da olmad?????? i??in Bioconductor'dan al??nd??.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# DESeq2 paketi R'da olmad?????? i??in Bioconductor'dan al??nd??.   

install.packages("tidyverse")
install.packages("pheatmap")

# Paketleri Aktif Etme

library(DESeq2)
library(org.Sc.sgd.db)
library(pheatmap)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(EnhancedVolcano)

# DESeq2 diferansiyel ekspresyon analizi i??in,
# org.Sc.sgd.db Saccharomyces cerevisiae gen veritaban?? i??in, 
# pheatmap ??s?? haritalar?? i??in, 
# tidyverse veri manip??lasyonu i??in, 
# clusterProfiler ve pathview ise zenginle??tirme analizleri i??in kullan??l??r.

######################## Veriyi Y??kleme ##################################

# Bu kod, DESeq2 analizi i??in gerekli gen say??m verilerini "gene_count_cutted.txt" dosyas??ndan okuyor.
# Bu i??lem i??in read.delim() fonksiyonu kullan??l??yor.
# ??lk s??tun sat??r isimleri olmas?? i??in row.names=1 parametresi giriliyor.

count_data <- read.delim("gene_count_cutted.txt", sep= "\t", row.names=1)

##################### Deney Ko??ullar??n??n Tan??mlanmas?? ######################

# DESeq objesi olu??turmak i??in gerekli ??zellik dosyas?? (colData) haz??rlan??yor.

sample_info <- data.frame(condition=c( "knockout", "knockout", "wild_type", "wild_type"))

rownames(sample_info) <- colnames(count_data)

# Bu kod, ko??ul verilerinin sat??r isimlerini, say??m verilerinin s??tun isimleriyle e??le??tiriyor.

## Veri Uyumlulu??unun Kontrol?? ##

# Bu sat??r, ko??ul verileri ile say??m verilerinin uyumlu olup olmad??????n?? kontrol ediyor.

all(rownames(sample_info)==colnames(count_data))

##################### DESeqDataSet Olu??turma ###############################

# DESeq2 analizi i??in gerekli olan DESeqDataSet nesnesi DESeqDataSetFromMatrix() fonksiyonu ile olu??turuluyor. 
# Gerekli olan parametreler; countData, colData ve design giriliyor.

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~condition)

dds

## D??????k Say??ml?? Genlerin Filtrelenmesi ## 

# Default olarak 10???dan az okuma say??s??na sahip genler kald??r??l??r.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Referans Seviyesinin Ayarlanmas?? ##

# Bu kod, ???wild_type??? ko??ulunu referans seviyesi olarak ayarl??yor.

dds$condition <- relevel(dds$condition, ref = "wild_type") 

############# Diferansiyel Ekspresyon Analizini ??al????t??ral??m ###############

# DESeq() fonksiyonu ile analiz ger??ekle??tiriliyor.

dds <- DESeq(dds) 

# DESeq2 analizinin sonu??lar?? results() fonksiyonu ile al??n??yor.

res <- results(dds)

res 

###################### PCA Grafiklerinin Olu??turulmas?? #####################

# PCA, boyut azaltma tekni??idir. Gen ekspresyon verisetindeki varyans?? a????klamak i??in kullan??l??r. 
# PCA olu??turmak i??in vst fonksiyonu ile varyans dengeleyici d??n??????m yap??l??r. (variance stabilizing transformation)

## Varianca Stabilizing Transformation ## 

# Bu blok, varyans stabilize edilmi?? verileri kullanarak bir PCA (Temel Bile??en Analizi) grafi??i olu??turuyor.

vstdata <- vst(dds, blind = F) 

# PCA Olu??tural??m

plotPCA(vstdata, intgroup=c("condition")) 

## MA Plot Olu??tural??m ## 

# Bu kod, bir MA (Mean-Average) grafi??i olu??turuyor. 
# alpha parametresi anlaml??l??k de??erini belirtiyor.

plotMA(dds, alpha=0.05) 

############## Sonu??lar??n Veri ??er??evesine D??n????t??r??lmesi ###################

# as.data.frame() fonksiyonu ile sonu??lar bir data frame d??n????t??r??l??r.
res_df <- as.data.frame(res) 

# Bu kod, d??zeltilmi?? p-de??eri 0.05???ten k??????k olan genleri anlaml?? olarak i??aretliyor. 
# Bunun i??in ifelse() fonksiyonundan faydalan??l??yor.

res_df$significant <- ifelse( res_df$padj < 0.05, "yes", "no") 

## plotMA Grafi??inin ggplot2 ile Tasarlan?????? ##

ggplot(res_df, aes(log(baseMean), log2FoldChange, color=significant)) +
  geom_point()

## Enhanced Volcano Grafi??i ##

# Volcano grafi??i i??in EnhancedVolcano() fonksiyonu kullan??l??r. 
# x ve y parametreleri x ekseni ve y eksenine hangi s??tunlar??n yerle??ece??ini belirtir. 
# lab parametresi etiket olarak kullan??lacak s??tunu belirtir.

EnhancedVolcano(res_df, x="log2FoldChange", y="padj", lab = rownames(res_df)) 

############### Normalize Edilmi?? Say??mlar??n Elde Edilmesi ###############

# 1. Diferansiyel ekspresyon sonu??lar??n?? (res_df) ve 
# VST ile stabilize edilmi?? verileri (vst_mat) haz??rlayal??m.
# assay() fonksiyonu, DESeq2 nesnesi gibi bir objenin i??indeki say??sal ifade matrisini (VST ile normalize edilmi?? gen ekspresyon de??erlerini) d??z bir matris olarak ????kar??r.

vst_mat <- assay(vstdata)

# 2. En anlaml?? ilk 20 geni filtreleyelim (padj ve log2FC kriterine g??re)

top_20_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(20) %>% 
  rownames()

# 3. Sadece bu 20 genin VST de??erlerini ??s?? haritas?? i??in ??ekelim
vst_top20 <- vst_mat[top_20_genes, ]

# 4. S??tun isimlerini SRR kodlar??ndan kurtar??p grup isimlerine (KO/WT) d??n????t??relim
colnames(vst_top20) <- sample_info$condition
colnames(vst_top20) <-sample_info$condition
######################## Heatmap Olu??turma ###############################

# Profesyonel Heatmap: VST verisi kullan??r, sat??r bazl?? scale yapar.

# scale="row"; Her geni kendi i??inde standartla??t??rarak (z-score) genler aras??ndaki de??i??im desenini vurgular.
# clustering_distance_rows = "euclidean"; Genleri (sat??rlar??) birbirine benzerliklerine g??re Euclidean mesafesi ile k??meleyerek gruplar.
# clustering_distance_cols = "euclidean"; ??rnekleri (s??tunlar??) yine Euclidean mesafesine g??re benzerliklerine g??re k??meler.

pheatmap(vst_top20, 
         scale = "row",                     
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Heatmap Plot")

#################### Gene Ontology ve KEGG Pathway #######################

# Gen seti zenginle??tirme analizi, gen ekspresyonu sonucunda ??retilen veri setlerinin kar????la??t??r??lmas?? ve yorumlanmas??n?? sa??lamaktad??r.
# Bu sayede gen s??n??flar??n??n fenotopik ili??kiler ile ba??lant??lar?? g??zlemlenebilir.

## Gene Ontology Enrichment ## 

# Gen ontoloji analizi, ifade edilen gen ve gen ??r??nleri i??levlerini a????klamay?? sa??lar. 
# Biyolojik anlamda ???? a????dan ele al??n??r: 
# 1) Moleculer Function: Molek??ler d??zeyde meydana gelen aktiviteleri tan??mlar. 
# 2) Cellular Component: Gen ??r??nlerinin konumlar??n?? tan??mlar. 
# 3) Biological Process: Biyolojik s??re??ler hakk??nda bilgi sa??lar.

res_df$significant <- ifelse(res_df$padj < 0.05, "yes", "no")

res_go <- res_df %>%
  filter(significant == "yes")

# Anlaml?? olan de??erleri sadece filtreledik. Analize anlaml?? olan genler ile devam edece??iz. 

gene_list <- rownames(res_go) 
head(gene_list)

## Biyolojik Process Ontoloji Analizi ##

# Gen ontoloji analizi i??in clusterProfiler paketinin sa??lam???? oldu??u enrichGO() fonksiyonu kullan??l??yor.
# gene parametresi gene listesini ister.
# OrgDb parametresi ??al??????lan veritaban??n?? ister. 
# keyType kullan??lan ID???lerin hangi veritabana ait oldu??unu ister. 
# ont ger??ekle??tirilecek ontoloji analizinin tipini ister. 
# Biyolojik process i??in BP girilir.

# OrgDb = org.Sc.sgd.db: Saccharomyces cerevisiae gen anotasyon veritaban??n?? kullan??r.
# keyType = "ORF": Kullan??lan gen ID tipinin ORF oldu??unu belirtir.
# ORF, bir organizmadaki protein kodlayan gen b??lgesinin standart isimlendirilmi?? kimli??idir (maya i??in YDR123W gibi).

ego1 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF", 
                 ont = "BP")

# Sonu??lar barplot() fonksiyonu ile g??rselle??tirilir.

barplot(ego1, title = "Biological Process")

## Molek??ler Function Ontoloji Analizi ## 

ego2 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "MF")

# Sonu??lar barplot() fonksiyonu ile g??rselle??tirilir.

barplot(ego2, title = "Molecular Function")

## Cellular Component Ontoloji Analizi ##

ego3 <- enrichGO(gene = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "CC")

# Sonu??lar barplot() fonksiyonu ile g??rselle??tirilir.

ego3_plot <- barplot(ego3, title = "Cellular Component") 

ego3_plot

# Bu kod, elde edilen grafiklerin kaydedilmesini sa??l??yor.

png("ego3_plot.png", width = 1200, height = 600, res = 150)
print(ego3_plot)
dev.off()

# ggrepel, ggplot2 grafiklerinde metin etiketlerini otomatik olarak ??ak????may?? ??nleyecek ??ekilde yerle??tirir.
# max.overlaps = Inf ile etiketlerin ??ak????mas?? g??z ard?? edilir, yani t??m etiketler mutlaka g??sterilir. 
# Normalde ??ok fazla ??ak????ma olursa baz?? etiketler gizlenir.
# goplot() fonksiyonu, GO terimlerini ve gen ili??kilerini bir a?? grafi??i ??eklinde g??rselle??tirir. 
# goplot () fonksiyonu GO analiz sonu??lar??n?? a?? grafi??i olarak ??iziyor.

options(ggrepel.max.overlaps = Inf)
goplot(ego1) 

## KEGG Yolu Zenginle??tirme Analizi ##

# KEGG, gen fonksiyonlar??n??n sistematik analizine y??nelik, biyolojik s??re??ler, hastal??klar, ila?? ara??t??rmas?? gibi ??nemli konular ??zerindeki ili??kilerin incelenmesini sa??lamaktad??r. 
# KEGG analizi de clusterProfiler paketi taraf??ndan desteklenmektedir. 
# ??al??????lan organizma KEGG taraf??ndan desteklenip desteklenmedi??ini ????renmek i??in search_kegg_organism() fonksiyonu kullan??l??r.

search <- search_kegg_organism("Saccharomyces cerevisiae") 

search

# Enrichment analizi i??in enrichKEGG() fonksiyonu kullan??l??r.
# organism = "sce".

kegg_enrich <- enrichKEGG(gene = gene_list,
                          organism = "sce",
                          keyType = "kegg")

head(kegg_enrich)

# G??rselle??tirme i??in browseKEGG() ve pathview() fonksiyonlar?? kullan??labilir. 

# ???sce00190???: ??ncelenecek KEGG yolak ID???si. 

browseKEGG(kegg_enrich, "sce00190")

# 1. Genlerin log2FoldChange de??erlerini alal??m
gen_verisi <- res_go$log2FoldChange

# 2. Bu de??erlere gen isimlerini (YDR... gibi) etiket olarak atayal??m
names(gen_verisi) <- rownames(res_go)

# 3. Pathview'i bu veriyle ??al????t??ral??m

sce00190<-pathview(gene.data = gen_verisi, 
         species = "sce", 
         pathway.id = "sce00190", 
         gene.idtype = "orf")

##############################################################################
##############################################################################