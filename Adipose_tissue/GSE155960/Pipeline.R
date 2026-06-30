´´´r
########################################################
## GSE155960 - Adipose Tissue (CD45+, lean vs. obese) scRNA-seq
## Vollstaendige Pipeline: Einlesen -> QC -> SCTransform -> Clustering ->
## Makrophagen-Annotation -> CaM-Score -> lean vs. obese -> Venn/Heatmap
##
## Methodik angelehnt an die GSE179640/GSE163973-Pipeline (CaM-Projekt).
## Wichtigste Abweichungen, da GSE155960 anders aufgebaut ist:
##  - Input: CSV.gz (Gene x Zellen, ENSEMBL-IDs) statt Read10X()/Read10X_h5()
##  - 6 Samples (CD45+ sortiert, 3x lean, 3x obese)
##  - ENSEMBL-IDs als rownames -> Symbol-Mapping noetig fuer alle Marker
##  - KEIN Harmony (UMAP/Sample-Plots zeigten bereits gute Durchmischung)
##  - Kein PTPRC-Cleanup (bereits CD45+ vorsortiert)
##  - Myeloid-Identifikation: UCell-Score allein war hier NICHT robust
##    (keine klare Bimodalitaet) -> finale Entscheidung ueber Marker-
##    Tabelle pro Cluster (CD68/CD14/CSF1R/LYZ/ITGAM/MRC1 etc.)
##  - percent.mt gilt hier als BIOLOGISCHES Signal -> keine Regression
########################################################


########################################################
## 0) Pakete laden
########################################################

library(Seurat)
library(Matrix)
library(data.table)
library(tidyverse)
library(future)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(UCell)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

if (!requireNamespace("effsize", quietly = TRUE)) install.packages("effsize")
library(effsize)

if (!requireNamespace("ggVennDiagram", quietly = TRUE)) install.packages("ggVennDiagram")
library(ggVennDiagram)

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

set.seed(42)

raw_dir <- "C:/Users/janin/Documents/In_Silico_CaM_Analysis/GSE155960_RAW"  # ggf. anpassen


########################################################
## 1) Sample-Metadaten
########################################################

sample_meta <- data.frame(
  gsm        = c("GSM4717158", "GSM4717159", "GSM4717160",
                 "GSM4717161", "GSM4717162", "GSM4717163"),
  file       = c("GSM4717158_CD45P-L1.csv.gz",
                 "GSM4717159_CD45P-L2.csv.gz",
                 "GSM4717160_CD45P-L3.csv.gz",
                 "GSM4717161_CD45P-O1.csv.gz",
                 "GSM4717162_CD45P-O2.csv.gz",
                 "GSM4717163_CD45P-O3.csv.gz"),
  sample     = c("L1", "L2", "L3", "O1", "O2", "O3"),
  condition  = c("lean", "lean", "lean", "obese", "obese", "obese"),
  stringsAsFactors = FALSE
)


########################################################
## 2) CSV.gz einlesen -> sparse Matrix -> Seurat-Objekte (6 Samples)
########################################################

read_csv_counts <- function(filepath) {
  message("Lese: ", filepath)
  dt <- fread(filepath, header = TRUE, sep = ",", showProgress = TRUE)

  gene_ids <- dt[[1]]
  cell_barcodes <- colnames(dt)[-1]
  dt[[1]] <- NULL

  mat <- Matrix(as.matrix(dt), sparse = TRUE)
  rownames(mat) <- gene_ids
  colnames(mat) <- cell_barcodes

  rm(dt)
  gc(verbose = FALSE)
  return(mat)
}

seurat_list <- list()

for (i in seq_len(nrow(sample_meta))) {
  gsm   <- sample_meta$gsm[i]
  fpath <- file.path(raw_dir, sample_meta$file[i])

  mat <- read_csv_counts(fpath)
  colnames(mat) <- paste0(gsm, "_", colnames(mat))

  obj <- CreateSeuratObject(counts = mat, project = gsm,
                             min.cells = 0, min.features = 0)

  obj$gsm       <- gsm
  obj$sample    <- sample_meta$sample[i]
  obj$condition <- sample_meta$condition[i]

  seurat_list[[gsm]] <- obj
  message(gsm, ": ", ncol(obj), " Zellen x ", nrow(obj), " Gene eingelesen.")
}


########################################################
## 3) Mergen zu einem Objekt
########################################################

all <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
all$condition <- factor(all$condition, levels = c("lean", "obese"))

message("Gesamt gemergt: ", ncol(all), " Zellen x ", nrow(all), " Gene")
table(all$sample)
table(all$condition)


########################################################
## 4) QC-Metriken (percent.mt ueber ENSEMBL-IDs)
########################################################

mt_ensembl_ids <- c(
  "ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712",
  "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840",
  "ENSG00000212907", "ENSG00000198886", "ENSG00000198786", "ENSG00000198695",
  "ENSG00000198727"
)
mt_present <- intersect(mt_ensembl_ids, rownames(all))
message(length(mt_present), " von 13 MT-Genen im Datensatz gefunden.")

all[["percent.mt"]] <- PercentageFeatureSet(all, features = mt_present)

VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "condition")

FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


########################################################
## 5) Filtering (einheitliche Cutoffs ueber alle Samples)
########################################################

all <- subset(
  all,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 4000 &
    percent.mt < 20
)

message("Nach Filterung: ", ncol(all), " Zellen x ", nrow(all), " Gene")
table(all$sample)
table(all$condition)


########################################################
## 6) Normalisierung (SCTransform)
## Bewusst OHNE vars.to.regress = "percent.mt": der MT%-Unterschied
## lean/obese gilt hier als biologisches Signal, nicht als Artefakt.
########################################################

options(future.globals.maxSize = 8 * 1024^3)
all <- SCTransform(all, assay = "RNA", verbose = TRUE)


########################################################
## 7) PCA
########################################################

all <- RunPCA(all, npcs = 30, verbose = TRUE)
ElbowPlot(all, ndims = 30)
DimPlot(all, reduction = "pca", group.by = "condition")
DimPlot(all, reduction = "pca", group.by = "sample")


########################################################
## 8) Clustering & UMAP (KEIN Harmony)
########################################################

dims_use <- 1:20

all <- FindNeighbors(all, dims = dims_use, reduction = "pca")
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = dims_use, reduction = "pca")

DimPlot(all, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(all, reduction = "umap", group.by = "sample")
DimPlot(all, reduction = "umap", group.by = "condition")

table(all$seurat_clusters)


########################################################
## 9) ENSEMBL -> Symbol Mapping
########################################################

ensembl_ids <- rownames(all[["RNA"]])

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_map_df <- data.frame(
  ensembl_id = ensembl_ids,
  symbol = symbol_map,
  stringsAsFactors = FALSE
)
saveRDS(gene_map_df, file = "ensembl_to_symbol_map_GSE155960.rds")

message(sum(!is.na(symbol_map)), " von ", length(ensembl_ids), " Genen gemappt.")

symbol_to_ensembl <- function(symbols) {
  hits <- gene_map_df %>% dplyr::filter(symbol %in% symbols)
  hits$ensembl_id
}


########################################################
## 10) Myeloid-Score (UCell) - EXPLORATIV
## ACHTUNG: dieser Score allein war NICHT ausreichend bimodal, um einen
## verlaesslichen globalen Cutoff zu setzen (siehe Diagnose unten) - er
## dient hier nur als Orientierung, die FINALE Entscheidung erfolgt ueber
## die Marker-Tabelle in Schritt 11.
########################################################

myeloid_markers_symbol_v2 <- c("CD68", "CD14", "LYZ", "CSF1R", "ITGAM")  # ohne FCGR3A (NK-Kontamination)
myeloid_ensembl_v2 <- symbol_to_ensembl(myeloid_markers_symbol_v2)
names(myeloid_ensembl_v2) <- gene_map_df$symbol[match(myeloid_ensembl_v2, gene_map_df$ensembl_id)]

Myeloid_sig_v2 <- list(Myeloid_v2 = myeloid_ensembl_v2)
DefaultAssay(all) <- "SCT"
all <- AddModuleScore_UCell(all, features = Myeloid_sig_v2)

hist(all$Myeloid_v2_UCell, breaks = 200,
     main = "Myeloid_v2_UCell - 200 Bins", xlab = "Myeloid_v2_UCell")

VlnPlot(all, features = "Myeloid_v2_UCell", group.by = "seurat_clusters", pt.size = 0) +
  ggtitle("Myeloid_v2_UCell pro Cluster")

score_per_cluster <- all@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    median_score = round(median(Myeloid_v2_UCell), 3),
    mean_score   = round(mean(Myeloid_v2_UCell), 3)
  ) %>%
  arrange(desc(median_score))

print(score_per_cluster, n = 25)


########################################################
## 11) Marker-Validierung auf ALLEN Clustern (Gesamt-Datensatz)
## Das ist die ROBUSTE Entscheidungsgrundlage fuer die Cluster-Auswahl,
## nicht der UCell-Score allein (der teils irrefuehrend war, z.B. weil
## ein urspruenglich "myeloid-hoher" Cluster sich als CD1C+ DC entpuppte).
########################################################

DefaultAssay(all) <- "RNA"

check_symbols <- c(
  "CD68", "CD14", "CSF1R", "ITGAM", "LYZ", "MRC1",   # Myeloid/Makrophage
  "CD1C", "CLEC9A", "LAMP3",                          # DC
  "CD3D", "CD3E",                                     # T-Zellen
  "MS4A1", "CD79A",                                   # B-Zellen
  "NCAM1", "GNLY", "NKG7"                             # NK
)

check_ensembl <- symbol_to_ensembl(check_symbols)
names(check_ensembl) <- gene_map_df$symbol[match(check_ensembl, gene_map_df$ensembl_id)]

marker_check <- FetchData(all, vars = c(check_ensembl, "seurat_clusters"))
colnames(marker_check)[seq_along(check_ensembl)] <- names(check_ensembl)

marker_summary <- marker_check %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    across(all_of(names(check_ensembl)),
           list(mean = ~round(mean(.x), 2), pct_pos = ~round(mean(.x > 0) * 100, 1)))
  ) %>%
  arrange(as.numeric(as.character(seurat_clusters)))

print(marker_summary, width = Inf, n = 25)
write.csv(marker_summary, "marker_check_GSE155960_allclusters.csv", row.names = FALSE)

# ACHTUNG: dplyr::select() explizit, da AnnotationDbi::select() (von org.Hs.eg.db)
# den dplyr-select() im Namespace ueberschattet und sonst einen Fehler wirft.
pct_pos_cols <- grep("_pct_pos$", colnames(marker_summary), value = TRUE)
compact_view <- marker_summary %>%
  dplyr::select(seurat_clusters, n_cells, all_of(pct_pos_cols))
colnames(compact_view) <- gsub("_pct_pos$", "", colnames(compact_view))
print(compact_view, width = Inf, n = 25)

p_dot <- DotPlot(all, features = check_ensembl, group.by = "seurat_clusters") +
  scale_x_discrete(labels = names(check_ensembl)) +
  RotatedAxis() +
  ggtitle("Marker-Expression pro Cluster (Gesamt-Datensatz)")
print(p_dot)

## Auswertung der marker_summary-Tabelle (GSE155960, 21 Cluster):
##   Cluster 3  (n=4262): CD68 88.4%, CSF1R 64.1%, LYZ 92.4%, MRC1 61.1%  -> Makrophage
##   Cluster 9  (n=1846): CD68 64.5%, CSF1R 45.1%, LYZ 99.1%, MRC1 59.9%  -> Makrophage
##   Cluster 15 (n=825):  CD68 91.8%, CSF1R 62.4%, LYZ 88.0%, MRC1 1.5%   -> Makrophage (inflamm., MRC1-low)
##   Cluster 19 (n=154):  CLEC9A 95.9% -> klassische DC1, KEIN Makrophage, separat halten
##   Cluster 20 (n=71):   CD68 68.2%, CD14 25.4%, CSF1R 0%, LYZ 13.6% -> schwaches/gemischtes
##                        Signal, trotzdem auf Nutzerwunsch aufgenommen
myeloid_cluster_ids <- c("3", "9", "15", "20")


########################################################
## 12) Finale Makrophagen-Extraktion + Lymphoid-Kontrollgruppe
########################################################

Idents(all) <- "seurat_clusters"

macrophages_final <- subset(all, idents = myeloid_cluster_ids)
lymphoid_control  <- subset(all, idents = setdiff(unique(all$seurat_clusters), myeloid_cluster_ids))

macrophages_final$coarse_cluster <- macrophages_final$seurat_clusters

message(ncol(macrophages_final), " Zellen als finale Makrophagen extrahiert.")
message(ncol(lymphoid_control), " Zellen als lymphoid_control (Rest) extrahiert.")

DimPlot(macrophages_final, reduction = "umap", group.by = "sample") +
  ggtitle("Finale Makrophagen - nach Sample")
DimPlot(macrophages_final, reduction = "umap", group.by = "condition") +
  ggtitle("Finale Makrophagen - nach Condition")

table(macrophages_final$sample)
table(macrophages_final$condition)
table(macrophages_final$coarse_cluster)

saveRDS(macrophages_final, file = "GSE155960_macrophages_final.rds")
saveRDS(lymphoid_control, file = "GSE155960_lymphoid_control.rds")


########################################################
## 13) CaM-Signatur definieren (40 Gene)
## 5 Paper-Anker (Murthy et al. 2022) + 35 funktional gruppierte,
## statistisch abgeleitete Gene (DEG CaM vs. GM-CSF UND CaM vs. M-CSF,
## aus der GSE180113-Re-Analyse)
########################################################

anchors <- c("SPP1", "TREM2", "FABP5", "CD63", "FABP4")
calcium_genes <- c("ATP2A2", "ATP2B1", "ATP2B4", "SLC8A1", "RCAN1")
migration_genes <- c("SEPTIN4", "SCIN", "DST", "KIF16B")
lysosomal_genes <- c("RAB38", "RAB7B", "SEC11C", "CCPG1", "GMPPB")
inflammation_genes <- c("GSDME", "BID", "IL18BP", "HPGDS", "NQO1")
immune_genes <- c("ABCG2", "CD300LB", "CLEC7A", "CD109", "PTPN22", "PHLDA1",
                   "HLA-DQB1", "DOK2", "DDAH2", "BST1", "ID3", "LHFPL2",
                   "LBH", "FN1", "HAMP", "NUCB2")

cam_genes_symbol <- c(anchors, calcium_genes, migration_genes, lysosomal_genes,
                       inflammation_genes, immune_genes)

length(cam_genes_symbol)          # erwartet: 40
length(unique(cam_genes_symbol))  # falls < 40 -> Duplikate

cam_lookup <- gene_map_df %>% dplyr::filter(symbol %in% cam_genes_symbol)
cam_ensembl_present <- intersect(cam_lookup$ensembl_id,
                                  rownames(macrophages_final[["RNA"]]))
names(cam_ensembl_present) <- gene_map_df$symbol[match(cam_ensembl_present,
                                                         gene_map_df$ensembl_id)]

message(length(cam_ensembl_present), " von ", length(cam_genes_symbol),
        " CaM-Genen im Datensatz vorhanden.")

CaM_sig <- list(CaM_Mac = cam_ensembl_present)


########################################################
## 14) CaM-Score berechnen + visualisieren
########################################################

DefaultAssay(macrophages_final) <- "SCT"
macrophages_final <- AddModuleScore_UCell(macrophages_final, features = CaM_sig)

hist(macrophages_final$CaM_Mac_UCell, breaks = 50,
     main = "CaM Module Score - finale Makrophagen GSE155960",
     xlab = "CaM_Mac_UCell Score")

FeaturePlot(macrophages_final, features = "CaM_Mac_UCell", reduction = "umap") +
  ggtitle("CaM Score - finale Makrophagen GSE155960")

VlnPlot(macrophages_final, features = "CaM_Mac_UCell",
        group.by = "coarse_cluster", pt.size = 0.1) +
  ggtitle("CaM Score pro Cluster (3/9/15/20) - Makrophagen GSE155960")


########################################################
## 15) Spezifitaetskontrolle: Makrophagen vs. lymphoid_control
########################################################

DefaultAssay(lymphoid_control) <- "SCT"
lymphoid_control <- AddModuleScore_UCell(lymphoid_control, features = CaM_sig)

specificity_df <- bind_rows(
  data.frame(score = macrophages_final$CaM_Mac_UCell, group = "Macrophage"),
  data.frame(score = lymphoid_control$CaM_Mac_UCell,  group = "Lymphoid_control")
)

ggplot(specificity_df, aes(x = group, y = score, fill = group)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  ggtitle("CaM-Score Spezifitaet: Makrophagen vs. Rest (lymphoid_control)") +
  theme_minimal()

wilcox_specificity <- wilcox.test(score ~ group, data = specificity_df)
print(wilcox_specificity)

specificity_summary <- specificity_df %>%
  group_by(group) %>%
  summarise(n = n(), median_score = round(median(score), 3), mean_score = round(mean(score), 3))
print(specificity_summary)

saveRDS(macrophages_final, file = "GSE155960_macrophages_final_with_CaM.rds")


########################################################
## 16) CaM-Score: lean vs. obese (Zell- und Sample-Ebene)
## ACHTUNG: dplyr::select() explizit (siehe Namespace-Hinweis oben)
########################################################

stopifnot("CaM_Mac_UCell" %in% colnames(macrophages_final@meta.data))

VlnPlot(macrophages_final, features = "CaM_Mac_UCell", group.by = "condition", pt.size = 0.1) +
  ggtitle("CaM-Score: lean vs. obese (Zellebene, alle Makrophagen)")

VlnPlot(macrophages_final, features = "CaM_Mac_UCell", group.by = "sample", pt.size = 0.1) +
  ggtitle("CaM-Score pro Sample (L1-L3 = lean, O1-O3 = obese)")

cam_per_cell <- macrophages_final@meta.data %>%
  dplyr::select(condition, sample, CaM_Mac_UCell)

# Zellebene (ACHTUNG Pseudoreplikation - immer mit Cliff's Delta berichten)
wilcox_zellebene <- wilcox.test(CaM_Mac_UCell ~ condition, data = cam_per_cell)
cliff_zellebene <- cliff.delta(CaM_Mac_UCell ~ condition, data = cam_per_cell)
print(wilcox_zellebene)
print(cliff_zellebene)

# Sample-Ebene (n=6, die eigentlich aussagekraeftige Stichprobengroesse)
cam_per_sample <- cam_per_cell %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(), mean_CaM = mean(CaM_Mac_UCell), median_CaM = median(CaM_Mac_UCell),
            .groups = "drop")
print(cam_per_sample)

ggplot(cam_per_sample, aes(x = condition, y = mean_CaM, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  geom_text(aes(label = sample), vjust = -1, size = 3) +
  ggtitle("Mittlerer CaM-Score pro Sample - lean vs. obese") +
  theme_minimal()

wilcox_sampleebene <- wilcox.test(mean_CaM ~ condition, data = cam_per_sample)
cliff_sampleebene <- cliff.delta(mean_CaM ~ condition, data = cam_per_sample)
print(wilcox_sampleebene)
print(cliff_sampleebene)

# Technischer Konfundierungs-Check
cor.test(macrophages_final$CaM_Mac_UCell, macrophages_final$nCount_RNA, method = "spearman")
cor.test(macrophages_final$CaM_Mac_UCell, macrophages_final$nFeature_RNA, method = "spearman")

cat("\n========================================\n")
cat("ZUSAMMENFASSUNG: CaM-Score lean vs. obese\n")
cat("========================================\n\n")
cat("Zellebene (n =", nrow(cam_per_cell), "Zellen, ACHTUNG Pseudoreplikation):\n")
cat("  Wilcoxon p-Wert:", format.pval(wilcox_zellebene$p.value), "\n")
cat("  Cliff's Delta:", round(cliff_zellebene$estimate, 3), "(", as.character(cliff_zellebene$magnitude), ")\n\n")
cat("Sample-Ebene (n = 6 Samples, 3 lean / 3 obese):\n")
cat("  Wilcoxon p-Wert:", format.pval(wilcox_sampleebene$p.value), "\n")
cat("  Cliff's Delta:", round(cliff_sampleebene$estimate, 3), "(", as.character(cliff_sampleebene$magnitude), ")\n\n")

write.csv(cam_per_sample, "GSE155960_CaM_score_per_sample.csv", row.names = FALSE)


########################################################
## 17) VENN-DIAGRAMM A: CaM-Gene detektierbar lean vs. obese
########################################################

DefaultAssay(macrophages_final) <- "RNA"
detect_threshold <- 0.01  # mind. 1% der Zellen muessen das Gen exprimieren

expr_by_condition <- FetchData(macrophages_final, vars = c(cam_ensembl_present, "condition"))
colnames(expr_by_condition)[seq_along(cam_ensembl_present)] <- names(cam_ensembl_present)

pct_detected <- expr_by_condition %>%
  group_by(condition) %>%
  summarise(across(all_of(names(cam_ensembl_present)), ~mean(.x > 0)))

genes_lean  <- names(cam_ensembl_present)[
  as.numeric(pct_detected[pct_detected$condition == "lean", names(cam_ensembl_present)]) > detect_threshold]
genes_obese <- names(cam_ensembl_present)[
  as.numeric(pct_detected[pct_detected$condition == "obese", names(cam_ensembl_present)]) > detect_threshold]

message("Lean: ", length(genes_lean), " / ", length(cam_ensembl_present), " Gene detektierbar")
message("Obese: ", length(genes_obese), " / ", length(cam_ensembl_present), " Gene detektierbar")

venn_data_condition <- list(lean = genes_lean, obese = genes_obese)

p_venn_condition <- ggVennDiagram(venn_data_condition, label = "count",
                                    category.names = c("lean", "obese")) +
  scale_fill_gradient(low = "#F8F8F8", high = "#4292C6") +
  ggtitle(paste0("CaM-Signatur (", length(cam_ensembl_present),
                 " Gene): Detektierbarkeit lean vs. obese\n(Cutoff: >",
                 detect_threshold*100, "% der Zellen exprimieren das Gen)"))
print(p_venn_condition)

only_obese <- setdiff(genes_obese, genes_lean)
only_lean  <- setdiff(genes_lean, genes_obese)
message("NUR in obese detektierbar: ", paste(only_obese, collapse = ", "))
message("NUR in lean detektierbar: ", paste(only_lean, collapse = ", "))


########################################################
## 18) VENN-DIAGRAMM B: CaM-Gene ueber die 4 Makrophagen-Cluster
########################################################

expr_by_cluster <- FetchData(macrophages_final, vars = c(cam_ensembl_present, "coarse_cluster"))
colnames(expr_by_cluster)[seq_along(cam_ensembl_present)] <- names(cam_ensembl_present)

pct_detected_cluster <- expr_by_cluster %>%
  group_by(coarse_cluster) %>%
  summarise(across(all_of(names(cam_ensembl_present)), ~mean(.x > 0)))

get_genes_for_cluster <- function(cl) {
  row <- pct_detected_cluster[pct_detected_cluster$coarse_cluster == cl, names(cam_ensembl_present)]
  names(cam_ensembl_present)[as.numeric(row) > detect_threshold]
}

venn_data_cluster <- list(
  Cluster3  = get_genes_for_cluster("3"),
  Cluster9  = get_genes_for_cluster("9"),
  Cluster15 = get_genes_for_cluster("15"),
  Cluster20 = get_genes_for_cluster("20")
)
lengths(venn_data_cluster)

p_venn_cluster <- ggVennDiagram(venn_data_cluster, label = "count",
                                  category.names = c("Cluster 3", "Cluster 9",
                                                      "Cluster 15", "Cluster 20")) +
  scale_fill_gradient(low = "#F8F8F8", high = "#74C476") +
  ggtitle(paste0("CaM-Signatur (", length(cam_ensembl_present),
                 " Gene): Detektierbarkeit ueber Makrophagen-Cluster"))
print(p_venn_cluster)


########################################################
## 19) HEATMAP: 40 CaM-Gene x lean/obese (Mittelwert-Expression)
## WICHTIG: assay = "SCT" (nicht "RNA"!) - der RNA-Assay hat hier nur
## rohe Counts, keinen normalisierten "data"-Layer (wir nutzen ja
## durchgehend SCTransform statt NormalizeData auf RNA).
########################################################

DefaultAssay(macrophages_final) <- "SCT"

avg_expr_condition <- AverageExpression(
  macrophages_final,
  features = cam_ensembl_present,
  group.by = "condition",
  assay = "SCT"
)$SCT

stopifnot(
  "avg_expr_condition ist leer - Assay/Layer pruefen!" = nrow(avg_expr_condition) > 0,
  "avg_expr_condition enthaelt nur NA-Werte!" = any(!is.na(avg_expr_condition))
)
message("avg_expr_condition: ", nrow(avg_expr_condition), " Gene x ",
        ncol(avg_expr_condition), " Bedingungen erfolgreich berechnet.")

rownames(avg_expr_condition) <- names(cam_ensembl_present)[
  match(rownames(avg_expr_condition), cam_ensembl_present)]
avg_expr_condition <- avg_expr_condition[!is.na(rownames(avg_expr_condition)), , drop = FALSE]

gene_order_available <- intersect(cam_genes_symbol, rownames(avg_expr_condition))
avg_expr_condition <- avg_expr_condition[gene_order_available, , drop = FALSE]

heatmap_matrix <- t(scale(t(avg_expr_condition)))

gene_group_lookup <- c(
  setNames(rep("Anchor", length(anchors)), anchors),
  setNames(rep("Calcium", length(calcium_genes)), calcium_genes),
  setNames(rep("Migration", length(migration_genes)), migration_genes),
  setNames(rep("Lysosomal", length(lysosomal_genes)), lysosomal_genes),
  setNames(rep("Inflammation", length(inflammation_genes)), inflammation_genes),
  setNames(rep("Immune", length(immune_genes)), immune_genes)
)

row_annotation <- data.frame(Funktion = gene_group_lookup[rownames(heatmap_matrix)])
rownames(row_annotation) <- rownames(heatmap_matrix)

pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation,
  main = "CaM-Signatur (40 Gene): lean vs. obese Makrophagen\n(z-Score pro Gen)",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  fontsize_row = 8,
  border_color = NA
)

print(round(avg_expr_condition, 3))
write.csv(avg_expr_condition, "GSE155960_CaM_genes_mean_expr_lean_obese.csv")

´´´
