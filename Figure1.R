library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

obj = readRDS('obj.raw.rds')

# Figure 1 ####
obj@meta.data$cell_subtype = case_when(
  obj@meta.data$cell_type == 'Tumor' ~ paste('Mal',obj@meta.data$cell_subtype,sep = '_'),
  .default = obj@meta.data$cell_subtype
)

cell_subtype = obj@meta.data|>
  group_by(cell_type)|>
  select(cell_subtype)|>
  distinct()|>
  ungroup()|>
  arrange(cell_type)|>
  pull(cell_subtype)

obj@meta.data$cell_subtype = factor(obj@meta.data$cell_subtype,levels = cell_subtype)


DimPlot(obj,group.by = c('cell_subtype'),pt.size = 1)+
  coord_fixed()+
  scale_color_igv()+
  theme_void()
ggsave('Figure1/umap_cell_subtype.pdf',width = 8,height = 6)


DimPlot(obj,group.by = 'cell_type',pt.size = 1,label = T,repel = T)+
  coord_fixed()+
  scale_color_npg()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')
ggsave('Figure1/umap_cell_type.pdf',width = 6,height = 6)


DotPlot(obj,group.by = 'cell_type',assay = 'SCT',
        features = c('FABP4','CD36','G0S2','ADIPOQ',
                     'COL1A1','DCN','FN1',
                     'PTPRCAP','CD3D','CD8B','IL32','CCL5',
                     'IGFBP7','PECAM1','ACTA2','RGS5','VWF',
                     'CD74','CD14','CD68','SPP1',
                     'IGKC','IGHG1','IGHA1','IGHM','MZB1',
                     'S100A6','CLU','IGF2','COL6A1','LGALS1'))+
  scale_color_viridis_c(option = 'H')+
  theme(aspect.ratio = 1/3,legend.direction = 'vertical',
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
ggsave('Figure1/dot_cell_marker.pdf',width = 10,height = 5)


ggplot(obj@meta.data,aes(Tumor_subtype))+
  geom_bar(aes(fill=cell_type),position = 'fill')+
  scale_fill_npg()+
  facet_wrap(~Tumor_subtype,scales = 'free')+
  coord_polar('y')+
  theme_void()
ggsave('Figure1/pie_tumor_subtype.pdf',width = 8,height = 6)


ggplot(obj@meta.data,aes(sample))+
  geom_bar(aes(fill=cell_type),position='fill')+
  facet_grid(~Tumor_subtype,space = 'free',scales = 'free')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
  scale_fill_npg()+
  labs(x='',y='proportion')
ggsave('Figure1/bar_sample_prop.pdf',width = 8,height = 6)


prop_data <- obj@meta.data |> 
  group_by(sample, Tumor_subtype) |>
  count(cell_type) |> 
  mutate(prop = n / sum(n)) |> 
  ungroup() |> 
  filter(Tumor_subtype %in% c('Lipoma','WDLPS','DDLPS')) |>
  select(cell_type, prop, sample) |>
  pivot_wider(names_from = 'cell_type', values_from = 'prop', values_fill = 0) |> 
  column_to_rownames('sample')

corr <- cor(prop_data, method = "pearson")
p.mat <- cor.mtest(prop_data, method = "pearson")$p
pdf('Figure1/cor_cell_type_LWD.pdf',width = 6,height = 6)
corrplot(
  corr,
  col = rev(COL2('RdBu')),
  tl.col = 'black',
  order = 'hclust',
  addCoef.col = 'black',
  addrect = 3,
  # p.mat = p.mat,
  sig.level = 0.05,
  insig = 'label_sig',
  pch.col = 'red',
  pch.cex = 1
)
dev.off()

sample_info = obj@meta.data|>select(sample,Tumor_subtype)|>distinct()
prop_data = prop_data|>rownames_to_column('sample')|>left_join(sample_info)

p1=ggscatter(data = prop_data,y='Adipocytes',x='CAFs',color = 'Tumor_subtype',
             add = 'reg.line',add.params = list(color = "blue", fill = "grey85",linetype=2),
             conf.int = TRUE,rug=TRUE,
             cor.coef = TRUE,
             cor.coeff.args = list(method = "pearson",label.sep = "\n"))+
  theme(aspect.ratio = 1)

p2=ggscatter(data = prop_data,y='Tumor',x='Myeloid cells',color = 'Tumor_subtype',
             add = 'reg.line',add.params = list(color = "blue", fill = "grey85",linetype=2),
             conf.int = TRUE,rug=TRUE,
             cor.coef = TRUE,
             cor.coeff.args = list(method = "pearson",label.sep = "\n"))+
  theme(aspect.ratio = 1)
p1+p2
ggsave('Figure1/cor_point_LWD.pdf',width = 8,height = 6)


saveRDS(obj,'obj.raw.rds')