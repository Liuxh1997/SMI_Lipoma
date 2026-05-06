library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

obj = readRDS('obj.raw.rds')

# Figure 2 ####
DEGs.B_M = FindMarkers(obj,ident.1 = 'Lipoma',group.by = 'Tumor_subtype')|>
  filter(p_val_adj < 0.05)|>arrange(avg_log2FC)
avg_exp = AverageExpression(obj,assays = 'SCT',features = rownames(DEGs.B_M),
                            group.by = 'Tumor_subtype')
pdf('Figure2/ht_tumor_sig.pdf',width = 4,height = 10)
pheatmap(avg_exp$SCT,scale = 'row',
         treeheight_row = 10,treeheight_col = 10,
         clustering_method = 'ward.D2',
         cluster_cols = F,
         border_color = NA,
         color = rev(COL2('RdYlBu')),
         cellwidth = 24,cellheight = 8,
         cutree_rows = 5,
         fontsize_row = 8)
dev.off()

obj@meta.data|>group_by(cell_type,Tumor_subtype)|>select(AP:Prot)|>
  summarise(across(everything(),median))|>
  pivot_longer(AP:Prot)|>
  ggplot(aes(name,value))+
  # geom_polygon(aes(color=cell_type,
  #                  group=cell_type),alpha=0.2,fill=NA)+
  geom_line(aes(color=cell_type,
                group=cell_type))+
  facet_wrap(~Tumor_subtype)+
  scale_color_npg()+
  coord_radial()+
  labs(x='',y='median of scocre')+
  theme(legend.position = 'top')
ggsave('Figure2/radar_tumor_sig.pdf',width = 8,height = 6)


row_cluster = read.csv('tumor.heterogeneity.degs.csv',row.names = 1)
subtype_signature = row_cluster|>group_by(clusters)|>
  summarise(genes = list(genes))|>
  deframe()|>as.list()
BP.list = lapply(subtype_signature, FUN = function(x){
  runGSEA(x,universe = rownames(obj),category = 'C5',subcategory = 'BP',pval.thr = 0.05)})
RC.list = lapply(subtype_signature, FUN = function(x){
  runGSEA(x,universe = rownames(obj),category = 'C2',subcategory = 'REACTOME',pval.thr = 0.05)})
HK.list = lapply(subtype_signature, FUN = function(x){
  runGSEA(x,universe = rownames(obj),category = 'H',pval.thr = 0.05)})


library(purrr)
combine_gsea <- function(gsea_list, category) {
  imap(gsea_list, ~ {
    .x$Type <- .y
    .x$Category <- category
    .x
  }) %>% list_rbind()
}
final_df <- list(
  combine_gsea(BP.list, "BP"),
  combine_gsea(RC.list, "REACTOME"),
  combine_gsea(HK.list, "H")
) %>% list_rbind()


final_df|>group_by(Type,Category)|>
  arrange(pval)|>
  slice_head(n=2)|>
  ggplot(aes(-log10(pval),pathway))+
  geom_col(aes(fill=Category))+
  scale_fill_aaas(alpha=0.8)+
  facet_wrap(~Type,scales = 'free_y',ncol=1,strip.position = 'right')+
  theme_bw()+
  theme(legend.position = 'bottom')
ggsave('Figure2/bar_sig_anno.pdf',width = 10,height = 6)


fov.show = obj@meta.data|>group_by(Tumor_subtype)|>
  count(fov)|>arrange(desc(n))|>
  slice_head(n=1)|>pull(fov)
filter(obj@meta.data,fov %in% c(61,49,104,91,17))|>
  ggplot(aes(x_FOV_px,-y_FOV_px))+
  geom_point(aes(color=cell_type),size=0.5,shape=16)+
  scale_color_npg()+
  facet_wrap(~Tumor_subtype,nrow=1,scales = 'free')+
  theme_dark()+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave('Figure2/point_spatial_subtype.pdf',width = 12,height = 4)


obj@meta.data|>select(fov,Tumor_subtype,x_FOV_px,y_FOV_px,AP:Prot)|>
  mutate(AP = scale(AP),
         Fat = scale(Fat),
         ECM = scale(ECM))|>
  filter(fov %in% c(61,49,104,91,17))|>
  pivot_longer(cols = c('AP','Fat','ECM'))|>
  ggplot(aes(x_FOV_px,-y_FOV_px))+
  geom_point(aes(color=value),shape=16,size=0.8,alpha=0.5)+
  # scale_color_npg()+
  scale_color_viridis_c(option = 'A')+
  facet_wrap(name~Tumor_subtype,nrow=3,scales = 'free')+
  theme_dark()+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave('Figure2/point_spatial_sig.pdf',width = 12,height = 10)


Idents(obj)
DotPlot(obj,group.by = 'cell_type',features = subtype_signature$ECM,idents = c('Fibroblasts','Tumor'),
        assay = 'RNA')+
  coord_flip()+
  theme_bw(base_line_size = 0)+
  theme(aspect.ratio = 3)
ggsave('Figure2/dot_ECM_genes.pdf',width = 6,height = 6)


VlnPlot(obj,idents = 'Adipocytes',group.by = 'Tumor_subtype',pt.size = 0,
        features = subtype_signature$Fat,adjust = 4,stack = T)+
  geom_boxplot(fill='white',width=0.2,outliers = F)
ggsave('Figure2/vln_Fat_genes.pdf',width = 8,height = 6)
