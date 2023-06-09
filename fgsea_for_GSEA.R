library(presto)
Treg.genes <- wilcoxauc(CD4Tnew3NAME, 'CellType')
head(Treg.genes)
dplyr::count(Treg.genes, group)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
msigdbr_show_species()
m_df<- msigdbr(species = "Homo sapiens", category = "C7")
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets$ALCALA_APOPTOSIS
Treg.genes %>%
  dplyr::filter(group == "Resting_Treg") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster2.genes<- Treg.genes %>%
  dplyr::filter(group == "Resting_Treg") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)
library("tibble")
ranks<- deframe(cluster2.genes)
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()
ggplot(fgseaResTidy %>% filter(padj < 0.005) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 7.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal()
plotEnrichment(fgsea_sets[["GSE15659_RESTING_VS_ACTIVATED_TREG_DN"]],
               ranks) + labs(title="ALCALA_APOPTOSIS")
MyPathways <- c("KEGG_APOPTOSIS","PID_FAS_PATHWAY","WP_TNFRELATED_WEAK_INDUCER_OF_APOPTOSIS_TWEAK_SIGNALING_PATHWAY","BIOCARTA_TNFR1_PATHWAY",
                 "REACTOME_APOPTOSIS","HAMAI_APOPTOSIS_VIA_TRAIL_UP","ALCALA_APOPTOSIS","WP_FERROPTOSIS","REACTOME_PYROPTOSIS")
MyPathways <- c("GSE15659_RESTING_VS_ACTIVATED_TREG_DN",
                 "GSE15659_RESTING_VS_ACTIVATED_TREG_UP")
plotGseaTable(fgsea_sets[MyPathways], ranks, fgseaResTidy,
              gseaParam=0.5)