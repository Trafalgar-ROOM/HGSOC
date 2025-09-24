library(org.Hs.eg.db)
library(clusterProfiler) # clusterProfiler v4.11.0

# group2
peak2_annot = read.table("group2_annotation.txt", header = T, sep = "\t")
peak2_ensembl = peak2_annot$geneId
peak2_entrez <- mapIds(org.Hs.eg.db,
                          keys = gsub("\\..*", "", peak2_ensembl),
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")
# GO enrichment (pvalues adjusted by 'BH' with cutoff < 0.05)
peak2_GO = enrichGO(
  gene = peak2_entrez,
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL", # includes BP (biological process), MF (molecular function), and CC (cellular component)
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# save the output
write.table(as.data.frame(peak2_GO), "peak2_GO_table.txt", sep = "\t")

# group1
peak1_annot = read.table("group1_annotation.txt", header = T, sep = "\t")
peak1_ensembl = peak1_annot$geneId
peak1_entrez <- mapIds(org.Hs.eg.db,
                       keys = gsub("\\..*", "", peak1_ensembl),
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
# GO enrichment
peak1_GO = enrichGO(
  gene = peak1_entrez,
  keyType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL", # includes BP (biological process), MF (molecular function), and CC (cellular component)
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# save the output
write.table(as.data.frame(peak1_GO), "peak1_GO_table.txt", sep = "\t")

group2_categories = c("norepinephrine transport", "axonogenesis", "regulation of norepinephrine secretion", "calcium ion transmembrane transport", "peptide hormone secretion", "manganese ion transmembrane transport", "potassium ion transmembrane transporter activity", "voltage-gated channel activity", "calcium ion transmembrane transporter activity", "high voltage-gated calcium channel activity")
group1_categories = c("leukocyte chemotaxis", "leukocyte migration", "leukocyte cell-cell adhesion", "myeloid leukocyte migration", "lymphocyte proliferation", "leukocyte proliferation", "lymphocyte differentiation", "response to chemokine", "regulation of leukocyte differentiation", "granulocyte chemotaxis")

# draw dotplot
categories = c(group1_categories, group2_categories)
pdf("Group1_immune_group2_ion_GO.pdf", width = 10, height = 8)
merge_result(list(Cluster_1 = peak1_GO, Cluster_2 = peak2_GO)) %>% dotplot(., showCategory = categories)
dev.off()
