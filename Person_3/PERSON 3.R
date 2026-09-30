# ==========================================
# Phase 1: Subsetting Malignant Cells & Re-clustering
# Person 3 - GSE131907 Analysis
# ==========================================

library(Seurat)
library(tidyverse)

# 1. ضبط مسار العمل
setwd("D:/Person 3")

# 2. قراءة كائن Seurat وملف الشخص الثاني
cat("--- جاري تحميل البيانات ---\n")
seu_obj <- readRDS("GSE131907_malignant_portable.rds")
handoff <- read.csv("25_Person2_final_malignancy_handoff.csv")

# 3. دمج تصنيف الشخص الثاني داخل Metadata الخاص بـ Seurat
cat("--- جاري دمج تصنيفات التسرطن ---\n")
metadata_merged <- seu_obj@meta.data %>%
  left_join(handoff %>% select(cell_id, Patient_ID, final_malignancy_class, CNV_score, CNV_label), by = "cell_id")

# إعادة تعيين الـ rownames لتطابق أسماء الخلايا في Seurat
rownames(metadata_merged) <- colnames(seu_obj)
seu_obj@meta.data <- metadata_merged

# 4. تصفية الخلايا الخبيثة المؤكدة فقط (Malignant supported)
cat("--- جاري تصفية الخلايا الخبيثة المؤكدة ---\n")
seu_mal <- subset(seu_obj, subset = final_malignancy_class == "Malignant supported")

cat("عدد الخلايا قبل التصفية:", ncol(seu_obj), "\n")
cat("عدد الخلايا الخبيثة المؤكدة بعد التصفية:", ncol(seu_mal), "\n")

# 5. إعادة المعالجة وخفض الأبعاد (Pipeline Standard)
cat("--- جاري معالجة البيانات وبناء UMAP ---\n")
seu_mal <- NormalizeData(seu_mal)
seu_mal <- FindVariableFeatures(seu_mal, selection.method = "vst", nfeatures = 2000)
seu_mal <- ScaleData(seu_mal)
seu_mal <- RunPCA(seu_mal, npcs = 30, verbose = FALSE)
seu_mal <- RunUMAP(seu_mal, dims = 1:20, verbose = FALSE)
seu_mal <- FindNeighbors(seu_mal, dims = 1:20, verbose = FALSE)
seu_mal <- FindClusters(seu_mal, resolution = 0.4, verbose = FALSE)

# TRY ANOTHER RESOULATION 
seu_mal <- FindClusters(seu_mal, resolution = 0.1, verbose = FALSE)
table(Idents(seu_mal))

# 6. حفظ الكائن النظيف للاستخدام في المراحل التالية
cat("--- جاري حفظ كائن Seurat النظيف ---\n")
saveRDS(seu_mal, "GSE131907_Person3_Malignant_Processed.rds")

cat("=== تم إنجاز المرحلة الأولى بنجاح! ===\n")


# 1. حساب الـ Cell Cycle Scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu_mal <- CellCycleScoring(seu_mal, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# 2. فحص توزيع مرحلة الدورة الخلوية بين الخلايا
cat("=== توزيع مراحل الدورة الخلوية (Cell Cycle Phase) ===\n")
print(table(seu_mal$Phase))

# 3. فحص نسبة الميتوكندريا (الإجهاد)
cat("\n=== متوسط نسبة الميتوكندريا (Stress level) ===\n")
print(summary(seu_mal$percent_mt))
# فحص عدد الخلايا من كل مريض داخل كل Cluster
cat("=== توزيع المرضى داخل الـ Clusters ===\n")
print(table(seu_mal$seurat_clusters, seu_mal$Patient_ID))

# ==========================================
# Patient Batch Correction (Harmony Integration)
# ==========================================

library(Seurat)
library(harmony)

# 1. تشغيل Harmony لإلغاء تأثير المريض (Patient_ID)
seu_mal <- RunHarmony(seu_mal, group.by.vars = "Patient_ID", plot_convergence = FALSE)

# 2. إعادة حساب UMAP و Clustering بناءً على أبعاد Harmony المدمجة
seu_mal <- RunUMAP(seu_mal, reduction = "harmony", dims = 1:20)
seu_mal <- FindNeighbors(seu_mal, reduction = "harmony", dims = 1:20)
seu_mal <- FindClusters(seu_mal, resolution = 0.1, verbose = FALSE)

# 3. فحص توزيع المرضى الجديد داخل الـ Clusters
cat("=== توزيع المرضى بعد الـ Harmony Integration ===\n")
print(table(seu_mal$seurat_clusters, seu_mal$Patient_ID))

# ==========================================
# Biological Program Scoring & Marker Detection
# ==========================================

library(Seurat)
library(tidyverse)

# 1. تعريف قوائم الجينات المرجعية (Signatures)
signatures <- list(
  EMT = c("VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "FN1", "MMP2"),
  Stemness = c("SOX2", "POU5F1", "NANOG", "PROM1", "ALDH1A1", "MYC", "KLF4"),
  Plasticity = c("EPCAM", "VIM", "CDH1", "CDH2", "SOX2", "ALDH1A1"),
  Proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "MCM6"),
  Lung_Diff = c("SFTPC", "SFTPA1", "SFTPB", "KRT5", "TP63", "NKX2-1")
)

# 2. حساب الـ Module Scores
cat("--- جاري حساب الـ Module Scores البيولوجية ---\n")
for(sig_name in names(signatures)) {
  genes <- list(intersect(signatures[[sig_name]], rownames(seu_mal)))
  if(length(genes[[1]]) > 0) {
    seu_mal <- AddModuleScore(seu_mal, features = genes, name = paste0(sig_name, "_Score"))
    col_idx <- grep(paste0("^", sig_name, "_Score"), colnames(seu_mal@meta.data))
    colnames(seu_mal@meta.data)[col_idx[length(col_idx)]] <- paste0(sig_name, "_Score")
  }
}

# 3. تحديد الجينات المميزة (Markers) لكل Cluster
cat("--- جاري تحديد الـ Marker Genes لتوصيف الحالات النسيجية ---\n")
cluster_markers <- FindAllMarkers(seu_mal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "GSE131907_Person3_Cluster_Markers.csv", row.names = FALSE)

cat("=== تم حساب الـ Scores واستخراج الـ Markers بنجاح! ===\n")

# ==========================================
# Biological States Assignment & Summary
# ==========================================

library(tidyverse)

# حساب متوسط كل برنامج بيولوجي لكل Cluster
scores_summary <- seu_mal@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Cell_Count = n(),
    Mean_EMT = mean(EMT_Score, na.rm = TRUE),
    Mean_Stemness = mean(Stemness_Score, na.rm = TRUE),
    Mean_Plasticity = mean(Plasticity_Score, na.rm = TRUE),
    Mean_Proliferation = mean(Proliferation_Score, na.rm = TRUE),
    Mean_Lung_Diff = mean(Lung_Diff_Score, na.rm = TRUE)
  )

print(as.data.frame(scores_summary))

# ==========================================
# Step 1: Assigning Biological States & Validation
# ==========================================

library(Seurat)
library(tidyverse)

# 1. تعريف مسميات الحالات البيولوجية بناءً على جدول الـ Scores
state_labels <- c(
  "0" = "Malignant_Hybrid_EMT_tS1",
  "1" = "Malignant_Differentiated_tS2",
  "2" = "Malignant_Cycling_tS3",
  "3" = "Exploratory_State_3",
  "4" = "Exploratory_High_EMT_4",
  "5" = "Exploratory_High_LungDiff_5",
  "6" = "Exploratory_Cycling_6",
  "7" = "Exploratory_State_7",
  "8" = "Exploratory_State_8",
  "9" = "Exploratory_High_Stemness_9",
  "10" = "Exploratory_State_10"
)

# 2. إضافة العمود للكائن
seu_mal$Malignant_State <- unname(state_labels[as.character(seu_mal$seurat_clusters)])

# 3. توثيق شرط نقطة التسليم (Core Shared vs Exploratory)
seu_mal$State_Validation <- ifelse(
  seu_mal$seurat_clusters %in% c("0", "1", "2"), 
  "Core_Shared_State", 
  "Patient_Specific_Exploratory"
)


cat("=== تم إسناد وتوثيق الحالات البيولوجية داخل الكائن بنجاح! ===\n")   


# ==========================================
# Step 2: Generating & Saving Visualizations
# ==========================================

library(Seurat)
library(ggplot2)

# 1. رسم UMAP حسب المرضى وحسب الحالات البيولوجية
p1 <- DimPlot(seu_mal, reduction = "umap", group.by = "Patient_ID") + 
  ggtitle("Malignant Cells by Patient (Post-Harmony)") +
  theme(plot.title = element_text(size = 12, face = "bold"))

p2 <- DimPlot(seu_mal, reduction = "umap", group.by = "Malignant_State", label = TRUE, repel = TRUE) + 
  ggtitle("Malignant Biological States (tS1-tS3 & Exploratory)") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# 2. رسم توزيع الـ Scores للبرامج البيولوجية
p3 <- FeaturePlot(
  seu_mal, 
  features = c("EMT_Score", "Stemness_Score", "Proliferation_Score", "Lung_Diff_Score"), 
  ncol = 2, 
  cols = c("lightgrey", "red")
)

# عرض الصور في شاشة RStudio
print(p1)
print(p2)
print(p3)

# 3. حفظ الصور بأعلى جودة في مجلد العمل للتقرير
ggsave("01_UMAP_by_Patient_PostHarmony.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("02_UMAP_by_Malignant_States.png", plot = p2, width = 9, height = 6, dpi = 300)
ggsave("03_FeaturePlot_Biological_Scores.png", plot = p3, width = 10, height = 8, dpi = 300)

cat("=== تم رسم وحفظ جميع الصور بنجاح في مجلد العمل! ===\n") 

# ==========================================
# Step 3: Pathway Enrichment Analysis for Core States
# ==========================================

library(tidyverse)

# 1. قراءة الجينات المميزة التي تم استخراجها سابقاً
markers <- read.csv("GSE131907_Person3_Cluster_Markers.csv")

# 2. تصفية أعلى الجينات تعبيراً لكل حالة معتمدة (p_val_adj < 0.05 & avg_log2FC > 0.5)
top_markers <- markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

cat("=== أهم الجينات المميزة للعناقيد الرئيسية ===\n")
print(top_markers %>% select(cluster, gene, avg_log2FC, p_val_adj), n = 30)

# حفظ القائمة المصفاة للملاحظات
write.csv(top_markers, "GSE131907_Person3_Top10_Markers_per_Cluster.csv", row.names = FALSE)


# ==========================================
# True Pathway Enrichment Analysis (GO / KEGG)
# ==========================================

# 1. تثبيت وتشغيل مكتبة enrichR لو مش موجودة
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)
library(tidyverse)

# 2. تحديد قواعد البيانات المرجعية (GO Biological Process & KEGG)
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")

# 3. قراءة الجينات المميزة للـ Core States (Clusters 0, 1, 2)
markers <- read.csv("GSE131907_Person3_Top10_Markers_per_Cluster.csv")

# أخذ الجينات المميزة للحالة tS1 (Cluster 0)
genes_tS1 <- markers %>% filter(cluster == 0) %>% pull(gene)
# أخذ الجينات المميزة للحالة tS3 (Cluster 2)
genes_tS3 <- markers %>% filter(cluster == 2) %>% pull(gene)

# 4. تشغيل الـ Enrichment لـ tS1 (EMT / Plasticity)
enriched_tS1 <- enrichr(genes_tS1, dbs)
# 5. تشغيل الـ Enrichment لـ tS3 (Cycling / Proliferation)
enriched_tS3 <- enrichr(genes_tS3, dbs)

# 6. عرض أعلى 5 مسارات حيوية لكل حالة
cat("=== Top Pathways for tS1 (Cluster 0) ===\n")
print(head(enriched_tS1[["GO_Biological_Process_2023"]][, c("Term", "Adjusted.P.value")], 5))

cat("\n=== Top Pathways for tS3 (Cluster 2) ===\n")
print(head(enriched_tS3[["GO_Biological_Process_2023"]][, c("Term", "Adjusted.P.value")], 5))

# 7. حفظ نتائج الـ Pathways في ملف للتقرير
write.csv(enriched_tS1[["GO_Biological_Process_2023"]], "GSE131907_Person3_Pathways_tS1.csv", row.names = FALSE)

# ==========================================
# Enrichment Analysis for tS2 (Cluster 1)
# ==========================================

# 1. أخذ الجينات المميزة للحالة tS2 (Cluster 1)
genes_tS2 <- markers %>% filter(cluster == 1) %>% pull(gene)

# 2. تشغيل الـ Enrichment لـ tS2
enriched_tS2 <- enrichr(genes_tS2, dbs)

# 3. عرض أعلى 5 مسارات حيوية لـ tS2
cat("=== Top Pathways for tS2 (Cluster 1 - Differentiated) ===\n")
print(head(enriched_tS2[["GO_Biological_Process_2023"]][, c("Term", "Adjusted.P.value")], 5))

# 4. حفظ نتائج الـ Pathways الخاصة بـ tS2
write.csv(enriched_tS2[["GO_Biological_Process_2023"]], "GSE131907_Person3_Pathways_tS2.csv", row.names = FALSE)
write.csv(enriched_tS3[["GO_Biological_Process_2023"]], "GSE131907_Person3_Pathways_tS3.csv", row.names = FALSE)

# ==========================================
# Generate Visual Plots for Pathway Enrichment
# ==========================================

library(ggplot2)
library(dplyr)

# 1. تجهيز بيانات المسارات للعناقيد الثلاثة
plot_pathways <- function(enrich_res, title, filename) {
  df <- enrich_res[["GO_Biological_Process_2023"]] %>%
    head(5) %>%
    mutate(Term = gsub(" \\(GO:.*\\)", "", Term)) # تنظيف اسم المسار
  
  p <- ggplot(df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = title, x = "Biological Process (GO)", y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 11, face = "bold"))
  
  ggsave(filename, plot = p, width = 8, height = 4, dpi = 300)
  return(p)
}

# 2. رسم وحفظ صور المسارات الحيوية
p_path_tS1 <- plot_pathways(enriched_tS1, "Top Biological Pathways for tS1 (EMT / Plastic)", "05_Pathways_tS1_EMT.png")
p_path_tS2 <- plot_pathways(enriched_tS2, "Top Biological Pathways for tS2 (Differentiated)", "06_Pathways_tS2_Diff.png")
p_path_tS3 <- plot_pathways(enriched_tS3, "Top Biological Pathways for tS3 (Cycling)", "07_Pathways_tS3_Cycling.png")

# عرض الصور في RStudio
print(p_path_tS1)
print(p_path_tS2)
print(p_path_tS3)

cat("=== تم رسم وحفظ صور الـ Pathway Enrichment بنجاح! ===\n")

# ==========================================
# Step 4: Trajectory Inference & Pseudotime
# ==========================================

library(Seurat)
library(ggplot2)

# 1. حساب المسافة النسخية (Pseudotime) انطلاقاً من حالة التمايز tS2
emb <- Embeddings(seu_mal, "harmony")
tS2_cells <- which(seu_mal$Malignant_State == "Malignant_Differentiated_tS2")

# تحديد مركز حالة tS2 كرباط أو نقطة بداية (Root)
root_point <- colMeans(emb[tS2_cells, 1:5])
dist_from_root <- apply(emb[, 1:5], 1, function(x) sqrt(sum((x - root_point)^2)))

# إضافة Pseudotime إلى Metadata
seu_mal$Pseudotime <- dist_from_root

# 2. رسم خريطة الـ Pseudotime على UMAP
p_traj <- FeaturePlot(seu_mal, features = "Pseudotime", cols = c("blue", "yellow", "red")) +
  ggtitle("Malignant Cell Trajectory (Pseudotime Progression from tS2)") +
  theme(plot.title = element_text(size = 11, face = "bold"))

# عرض الصورة في RStudio
print(p_traj)

# ==========================================
# Proper Trajectory Inference using Slingshot in R
# ==========================================

# 1. تثبيت وتشغيل مكتبة Slingshot لو مش عندك
if (!requireNamespace("slingshot", quietly = TRUE)) {
  BiocManager::install("slingshot")
}
library(slingshot)
library(Seurat)
library(ggplot2)

# 2. تحويل كائن Seurat إلى SingleCellExperiment وتطبيق Slingshot
sce <- as.SingleCellExperiment(seu_mal)

# تشغيل Slingshot وتحديد بداية المسار (Root) من حالة tS2 (Cluster 1)
sce_sling <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'HARMONY', start.clus = '1')

# 3. استخراج قيم الـ Pseudotime الحقيقية وإضافتها لـ Seurat
seu_mal$Pseudotime_Slingshot <- slingPseudotime(sce_sling)[, 1]

# 4. رسم الـ Pseudotime الانسيابي الحقيقي وحفظه
p_traj_correct <- FeaturePlot(seu_mal, features = "Pseudotime_Slingshot", cols = c("blue", "yellow", "red")) +
  ggtitle("Malignant Trajectory Inference (Slingshot Pseudotime)") +
  theme(plot.title = element_text(size = 11, face = "bold"))

print(p_traj_correct)
ggsave("04_Malignant_Trajectory_Slingshot.png", plot = p_traj_correct, width = 8, height = 6, dpi = 300)

cat("=== تم حساب ورسم الـ Trajectory المنهجي المظبوط بنجاح! ===\n")

# 3. حفظ الرسمة كصورة رسمية
ggsave("04_Malignant_Trajectory_Pseudotime.png", plot = p_traj, width = 8, height = 6, dpi = 300)

cat("=== تم حساب ورسم مسار التطور النسيجي (Trajectory / Pseudotime) بنجاح! ===\n")  


# ==========================================
# Final Export of Person 3 Deliverables
# ==========================================

library(Seurat)
library(tidyverse)

# 1. تصدير جدول تركيب العناقيد حسب المريض والموقع (Cluster Composition)
composition_table <- table(seu_mal$Malignant_State, seu_mal$Patient_ID)
write.csv(as.data.frame.matrix(composition_table), "GSE131907_Person3_Cluster_Composition_by_Patient.csv")

# 2. تجهيز جدول التسليم للشخص الرابع (Handoff Scores & Labels Table)
handoff_p4 <- seu_mal@meta.data %>%
  select(
    cell_id, Patient_ID, Sample, Sample_Origin, seurat_clusters,
    Malignant_State, State_Validation, Phase, percent_mt,
    EMT_Score, Stemness_Score, Plasticity_Score, Proliferation_Score, Lung_Diff_Score, Pseudotime_Slingshot
  )

write.csv(handoff_p4, "30_Person3_final_scores_handoff.csv", row.names = FALSE)

# 3. حفظ كائن Seurat النهائي المكتمل للدور الثالث
saveRDS(seu_mal, "GSE131907_Person3_Malignant_Final.rds")

cat("\n=======================================================\n")
cat("🎉 تم تصدير كافة مخرجات وتكليفات الشخص الثالث بنجاح تام! 🎉\n")
cat("1. Handoff File: 30_Person3_final_scores_handoff.csv\n")
cat("2. Cluster Composition: GSE131907_Person3_Cluster_Composition_by_Patient.csv\n")
cat("3. Final RDS Object: GSE131907_Person3_Malignant_Final.rds\n")
cat("=======================================================\n")



# ==========================================
# Person 3 - Validation & Core Shared Verification
# ==========================================

library(Seurat)
library(tidyverse)

setwd("D:/Person 3")

# 1. تحميل الكائن النهائي لـ Person 3 وحجم بيانات Person 2
seu_mal <- readRDS("GSE131907_Person3_Malignant_Processed.rds")
handoff_p2 <- read.csv("25_Person2_final_malignancy_handoff.csv")

# أ) التحقق من مطابقة الـ 889 خلية مع Person 2
p2_supported_count <- sum(handoff_p2$final_malignancy_class == "Malignant supported", na.rm = TRUE)
cat("=== 1. التحقق من عدد الخلايا الخبيثة ===\n")
cat("عدد الخلايا الداعمة للتسرطن لدى Person 2:", p2_supported_count, "\n")
cat("عدد الخلايا في كائن Person 3 الحالي:", ncol(seu_mal), "\n")
if(p2_supported_count == ncol(seu_mal)) {
  cat("--> النتيجة: مطابقة كاملة 100%! لا يوجد فقدان في الـ cell_id.\n\n")
} else {
  cat("--> تنبيه: يوجد اختلاف يتطلب مراجعة عملية الـ Join!\n\n")
}

# ب) فحص توزيع المرضى على الـ Clusters (Core Shared Validation)
cat("=== 2. جدول توزيع المرضى على العناقيد (Cluster vs Patient) ===\n")
patient_cluster_tbl <- table(Idents(seu_mal), seu_mal$Patient_ID)
print(patient_cluster_tbl)

cat("\n=== 3. نسبة المساهمة النسبية للمرضى لكل Cluster (%) ===\n")
print(round(prop.table(patient_cluster_tbl, margin = 1) * 100, 1))

# ج) فحص توزيع موقع العينات (Primary vs Metastatic)
if("Sample_Origin" %in% colnames(seu_mal@meta.data)) {
  cat("\n=== 4. توزيع موقع الورم على العناقيد (Cluster vs Sample Origin) ===\n")
  print(table(Idents(seu_mal), seu_mal$Sample_Origin))
}


# ==============================================================================
# === 5. تحليل الحساسية ومقارنة الاتجاهات (Full Dataset vs. 889 Malignant) ===
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# 1. ربط الكائن الكلي (31,136 خلية)
seu_full <- seu_obj

# 2. التأكد الصريح من تجهيز الـ Layers وتفعيل الـ Normalization للداتا الكلية
seu_full <- NormalizeData(seu_full, verbose = FALSE)

# 3. تحديد وقوائم الجينات وتصفيتها لتناسب rownames الخاصة بـ seu_full و seu_mal
emt_genes_full        <- intersect(signatures$EMT, rownames(seu_full))
plasticity_genes_full <- intersect(signatures$Plasticity, rownames(seu_full))

emt_genes_mal        <- intersect(signatures$EMT, rownames(seu_mal))
plasticity_genes_mal <- intersect(signatures$Plasticity, rownames(seu_mal))

cta_genes <- c("MAGEA3", "MAGEA4", "NY-ESO-1", "CTAG1B", "CTAG2", "GAGE1", "SSX2")
cta_genes_full <- intersect(cta_genes, rownames(seu_full))
cta_genes_mal  <- intersect(cta_genes, rownames(seu_mal))

cat("--- جاري حساب الـ Signatures للداتا الكاملة (31k) ---\n")
if(length(emt_genes_full) > 0) {
  seu_full <- AddModuleScore(seu_full, features = list(emt_genes_full), name = "EMT_Score_Full")
}
if(length(plasticity_genes_full) > 0) {
  seu_full <- AddModuleScore(seu_full, features = list(plasticity_genes_full), name = "Plasticity_Score_Full")
}
if(length(cta_genes_full) > 0) {
  seu_full <- AddModuleScore(seu_full, features = list(cta_genes_full), name = "CTA_Score_Full")
}

cat("--- جاري حساب الـ Signatures للخلايا الخبيثة (889) ---\n")
if(length(emt_genes_mal) > 0) {
  seu_mal <- AddModuleScore(seu_mal, features = list(emt_genes_mal), name = "EMT_Score_Mal")
}
if(length(plasticity_genes_mal) > 0) {
  seu_mal <- AddModuleScore(seu_mal, features = list(plasticity_genes_mal), name = "Plasticity_Score_Mal")
}
if(length(cta_genes_mal) > 0) {
  seu_mal <- AddModuleScore(seu_mal, features = list(cta_genes_mal), name = "CTA_Score_Mal")
}

cat("\n=== تم حساب الـ Signatures بنجاح للمجموعتين ===\n")

# --- الخطوة الثانية: استخراج ومقارنة المتوسطات ---

cat("\n=== 5.1 مقارنة متوسطات الدرجات وحجم التأثير ===\n")

# استخراج اسم العمود المضاف تلقائياً بواسطة Seurat
col_emt_full <- grep("^EMT_Score_Full", colnames(seu_full@meta.data), value = TRUE)[1]
col_emt_mal  <- grep("^EMT_Score_Mal", colnames(seu_mal@meta.data), value = TRUE)[1]

emt_full_vals <- seu_full@meta.data[[col_emt_full]]
emt_mal_vals  <- seu_mal@meta.data[[col_emt_mal]]

emt_full_mean <- mean(emt_full_vals, na.rm = TRUE)
emt_mal_mean  <- mean(emt_mal_vals, na.rm = TRUE)

cat("EMT Score - Full (31k):", round(emt_full_mean, 4), "\n")
cat("EMT Score - Malignant (889):", round(emt_mal_mean, 4), "\n")

# حساب Cohen's d
cohen_d_emt <- (emt_mal_mean - emt_full_mean) / sqrt((sd(emt_full_vals, na.rm=T)^2 + sd(emt_mal_vals, na.rm=T)^2) / 2)
cat("Cohen's d (Effect Size) for EMT:", round(cohen_d_emt, 4), "\n")

# --- الخطوة الثالثة: تجميع البيانات ورسم المخطط البياني (Violin Plots) ---

col_plas_full <- grep("^Plasticity_Score_Full", colnames(seu_full@meta.data), value = TRUE)[1]
col_plas_mal  <- grep("^Plasticity_Score_Mal", colnames(seu_mal@meta.data), value = TRUE)[1]

df_comp <- bind_rows(
  data.frame(Score = emt_full_vals, Type = "EMT", Dataset = "Full (31k)"),
  data.frame(Score = emt_mal_vals, Type = "EMT", Dataset = "Malignant (889)"),
  data.frame(Score = seu_full@meta.data[[col_plas_full]], Type = "Plasticity", Dataset = "Full (31k)"),
  data.frame(Score = seu_mal@meta.data[[col_plas_mal]], Type = "Plasticity", Dataset = "Malignant (889)")
)

col_cta_full <- grep("^CTA_Score_Full", colnames(seu_full@meta.data), value = TRUE)
col_cta_mal  <- grep("^CTA_Score_Mal", colnames(seu_mal@meta.data), value = TRUE)

if(length(col_cta_full) > 0 && length(col_cta_mal) > 0) {
  df_cta <- bind_rows(
    data.frame(Score = seu_full@meta.data[[col_cta_full[1]]], Type = "CTA", Dataset = "Full (31k)"),
    data.frame(Score = seu_mal@meta.data[[col_cta_mal[1]]], Type = "CTA", Dataset = "Malignant (889)")
  )
  df_comp <- bind_rows(df_comp, df_cta)
}

# رسم الشكل المجمع
p_sensitivity <- ggplot(df_comp, aes(x = Dataset, y = Score, fill = Dataset)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("Full (31k)" = "#2b5c8f", "Malignant (889)" = "#d95f02")) +
  labs(title = "Sensitivity Analysis: Full Dataset vs. Malignant Supported",
       y = "Module Score / Signature Expression",
       x = "") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold", size = 12))

print(p_sensitivity)

# حفظ الرسم البياني
ggsave("Sensitivity_Analysis_Full_vs_Malignant.png", plot = p_sensitivity, width = 9, height = 5, dpi = 300)
cat("\n=== تم تنفيذ تحليل الحساسية وحفظ الرسم البياني بنجاح! ===\n")
}

cat("--- جاري حساب الـ Signatures للخلايا الخبيثة (889) ---\n")
seu_mal <- AddModuleScore(seu_mal, features = list(emt_genes), name = "EMT_Score_Mal")
seu_mal <- AddModuleScore(seu_mal, features = list(plasticity_genes), name = "Plasticity_Score_Mal")

if(length(cta_genes_present_mal) > 0) {
  seu_mal <- AddModuleScore(seu_mal, features = list(cta_genes_present_mal), name = "CTA_Score_Mal")
}

cat("\n=== تم حساب الـ Signatures بنجاح للمجموعتين ===\n")

# --- الخطوة الثانية: اختبار اتجاهات الـ CTA و EMT و Plasticity ---

cat("\n=== 5.1 مقارنة متوسطات الدرجات بين المجموعتين ===\n")

# حساب المتوسطات (Seurat يضيف رقم 1 تلقائياً لاسم العمود)
emt_full_mean <- mean(seu_full$EMT_Score_Full1, na.rm = TRUE)
emt_mal_mean  <- mean(seu_mal$EMT_Score_Mal1, na.rm = TRUE)

cat("EMT Score - Full (31k):", round(emt_full_mean, 4), "\n")
cat("EMT Score - Malignant (889):", round(emt_mal_mean, 4), "\n")

# حساب حجم التأثير (Cohen's d)
cohen_d_emt <- (emt_mal_mean - emt_full_mean) / sqrt((sd(seu_full$EMT_Score_Full1, na.rm=T)^2 + sd(seu_mal$EMT_Score_Mal1, na.rm=T)^2) / 2)
cat("Cohen's d (Effect Size) for EMT:", round(cohen_d_emt, 4), "\n")

# --- الخطوة الثالثة: الرسم البياني المقارن (Sensitivity Violin Plots) ---

df_comp <- bind_rows(
  data.frame(Score = seu_full$EMT_Score_Full1, Type = "EMT", Dataset = "Full (31k)"),
  data.frame(Score = seu_mal$EMT_Score_Mal1, Type = "EMT", Dataset = "Malignant (889)"),
  data.frame(Score = seu_full$Plasticity_Score_Full1, Type = "Plasticity", Dataset = "Full (31k)"),
  data.frame(Score = seu_mal$Plasticity_Score_Mal1, Type = "Plasticity", Dataset = "Malignant (889)")
)

if("CTA_Score_Full1" %in% colnames(seu_full@meta.data) && "CTA_Score_Mal1" %in% colnames(seu_mal@meta.data)) {
  df_cta <- bind_rows(
    data.frame(Score = seu_full$CTA_Score_Full1, Type = "CTA", Dataset = "Full (31k)"),
    data.frame(Score = seu_mal$CTA_Score_Mal1, Type = "CTA", Dataset = "Malignant (889)")
  )
  df_comp <- bind_rows(df_comp, df_cta)
}

# رسم المخطط البياني
p_sensitivity <- ggplot(df_comp, aes(x = Dataset, y = Score, fill = Dataset)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("Full (31k)" = "#2b5c8f", "Malignant (889)" = "#d95f02")) +
  labs(title = "Sensitivity Analysis: Full Dataset vs. Malignant Supported",
       y = "Module Score / Signature Expression",
       x = "") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold", size = 12))

print(p_sensitivity)

# حفظ الرسم البياني
ggsave("Sensitivity_Analysis_Full_vs_Malignant.png", plot = p_sensitivity, width = 9, height = 5, dpi = 300)
cat("\n=== تم تنفيذ تحليل الحساسية وحفظ الرسم البياني بنجاح! ==\n")






# تصدير ملف الـ Scores Handoff بالأعمدة المطلوبة فقط
# تصدير ملف الـ Scores Handoff بالأعمدة المطلوبة بالضبط
# ==============================================================================
# 1. حساب باقي الـ Scores والدورة الخلوية للـ 31,136 خلية بالكامل
# ==============================================================================
library(Seurat)

cat("--- جاري تكميل حساب باقي الـ Scores لجميع الخلايا (31k) ---\n")

# أ) حساب الدورة الخلوية (Phase)
seu_full <- CellCycleScoring(
  seu_full, 
  s.features = cc.genes$s.genes, 
  g2m.features = cc.genes$g2m.genes, 
  set.ident = FALSE
)

# ب) حساب باقي الـ Module Scores
seu_full <- AddModuleScore(seu_full, features = list(signatures$Stemness), name = "Stemness_Score_Full")
seu_full <- AddModuleScore(seu_full, features = list(signatures$Proliferation), name = "Proliferation_Score_Full")
seu_full <- AddModuleScore(seu_full, features = list(signatures$Lung_Diff), name = "Lung_Diff_Score_Full")

cat("✔ تم الحساب بنجاح لكل الخلايا!\n")

# ==============================================================================
# 2. تصدير ملف الـ Handoff المكتمل 100% (صفر NA)
# ==============================================================================
meta <- seu_full@meta.data

# استخراج الأعمدة المحسوبة حديثاً
col_stem  <- grep("^Stemness_Score_Full", colnames(meta), value = TRUE)[1]
col_prolif <- grep("^Proliferation_Score_Full", colnames(meta), value = TRUE)[1]
col_lung   <- grep("^Lung_Diff_Score_Full", colnames(meta), value = TRUE)[1]

df_handoff_clean <- data.frame(
  cell_id               = rownames(meta),
  Patient_ID            = meta$Patient_ID,
  Sample                = meta$Sample,
  Sample_Origin         = meta$Sample_Origin,
  original_Cell_subtype = meta$Cell_subtype,
  EMT_Score             = meta$EMT_Score_Full1,
  Stemness_Score        = meta[[col_stem]],
  Plasticity_Score      = meta$Plasticity_Score_Full1,
  Proliferation_Score   = meta[[col_prolif]],
  Lung_Diff_Score       = meta[[col_lung]],
  Phase                 = meta$Phase,
  percent_mt            = meta$percent_mt,
  stringsAsFactors      = FALSE
)

# حفظ ملف الـ CSV المكتمل
write.csv(df_handoff_clean, "GSE131907_Person3_AllCells_Scores_Handoff.csv", row.names = FALSE)

# حفظ كائن Seurat الكلي المكتمل
saveRDS(seu_full, "GSE131907_Person3_AllCells_Final.rds")



# حفظ الملف
write.csv(df_handoff, "GSE131907_Person3_AllCells_Scores_Handoff.csv", row.names = FALSE)
cat("✔ تم حفظ الملف بنجاح بـ 31,136 صف وبالأعمدة المطلوبة بالضبط!\n")


# 2. حفظ كائن Seurat الكامل
saveRDS(seu_full, "GSE131907_Person3_AllCells_Final.rds")

# 3. حفظ جدول الجينات
df_signatures <- data.frame(
  signature_name = c(rep("EMT", length(emt_genes)), rep("Plasticity", length(plasticity_genes)), rep("CTA", length(cta_genes_present_full))),
  gene = c(emt_genes, plasticity_genes, cta_genes_present_full),
  reference = "Published Literature"
)
write.csv(df_signatures, "Person3_signature_gene_sets.csv", row.names = FALSE)
