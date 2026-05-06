library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

# Figure 4 ####
obj = readRDS('obj.raw.rds')
all = obj
# umap
obj = subset(obj,subset=cell_type=='Mesenchymal cells')
obj = RunPCA(obj,npcs = 30)
obj = RunHarmony(obj,group.by.vars='sample')
obj = RunUMAP(obj,reduction = 'harmony',dims = 1:10,
              umap.method = 'umap-learn',metric = 'correlation')

DimPlot(obj,group.by = 'cell_subtype',
        label = T,repel = T,label.size = 6,
        alpha = 0.5,
        raster = T,raster.dpi = c(1060,1060),
        pt.size = 4)+
  scale_color_frontiers()+
  coord_equal()
ggsave('Figure4/umap_mes_subtype.pdf',width = 8,height = 6)


DimPlot(obj,split.by = 'Tumor_subtype',group.by = 'cell_subtype',
        ncol=3,raster = T,pt.size = 3)+
  coord_fixed()+
  theme(legend.position = c(0.8,0.2))
ggsave('Supplementary/FigS4_subcluster_umap.pdf',width = 8,height = 6)


obj@meta.data|>
  group_by(Tumor_subtype)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  ggplot(aes(Tumor_subtype,prop))+
  geom_col(aes(fill=cell_subtype))+
  labs(title = 'Mesenchymal cells',y='proportion')
ggsave('Supplementary/FigS4_subcluster_prop.pdf',width = 8,height = 6)


# Dot
Idents(obj)<-'cell_subtype'
obj=JoinLayers(obj,'RNA')
DEGs = FindAllMarkers(obj,assay = 'RNA',slot = 'counts',
                      logfc.threshold = 0.5,min.diff.pct = 0.05,
                      only.pos = T,return.thresh = 0.05)
write.csv(DEGs,'Supplementary/TableS_Mes_subcluster_DEGs.csv',row.names = F)


features = DEGs|>group_by(cluster)|>
  slice_head(n=3)|>
  pull(gene)

DotPlot(obj,assay = 'RNA',features = unique(features))+
  coord_flip()+
  theme_bw(base_line_size = 0)+
  scale_color_viridis_c(option = 'H')+
  theme(aspect.ratio = 3,axis.text.x = element_text(angle =90,vjust = 0.5,hjust = 1))
ggsave('Figure4/dot_mes_genes.pdf',width = 6,height = 8)


# prop
all@meta.data|>
  filter(Tumor_subtype%in%c('DDLPS','PLPS','MLPS'))|>
  group_by(Tumor_subtype,fov)|>
  count(cell_subtype)|>
  mutate(prop = n/sum(n))|>
  filter(cell_subtype%in%c('Endo_PECAM1','Myo_ACTA2','Pericytes',
                           'Endo_IGFBP3','Myo_PFN1','Endo_CCL21'))|>
  ggplot(aes(Tumor_subtype,prop))+
  geom_boxplot(aes(fill=Tumor_subtype),
               staplewidth = 0.8,width=0.5,
               alpha=0.2,outliers = F)+
  stat_compare_means()+
  geom_point(aes(color=Tumor_subtype),
             position = position_jitter(width = 0.2),
             shape=16)+
  facet_wrap(~cell_subtype,scales = 'free')+
  theme_bw(base_size = 12)+
  theme(legend.position = 'top')+
  labs(x='',y='Cell density')
ggsave('Figure4/box_mes_density.pdf',width = 6,height = 6)


# cell num
obj@meta.data %>%
  group_by(Tumor_subtype) %>%
  summarise(
    n = n_distinct(fov),
    cells = n(),
    cells.fov = cells / n
  )|>
  ggplot(aes(Tumor_subtype,cells.fov))+
  geom_col(aes(fill=Tumor_subtype))+
  theme(aspect.ratio = 1.5)
ggsave('Figure4/col_fov_cells.pdf',width = 6,height = 8)


# fov show
df=all@meta.data|>
  filter(Tumor_subtype%in%c('DDLPS','PLPS','MLPS'))|>
  group_by(Tumor_subtype,fov)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  filter(cell_subtype%in%c('Endo_PECAM1','Myo_ACTA2','Pericytes'))|>
  arrange(-prop)|>
  group_by(Tumor_subtype,cell_subtype)|>
  slice_head(n=3)


all@meta.data|>
  filter(fov%in%c(101,147,21))|>
  mutate(cell_subtype=case_when(
    cell_subtype%in%c('Endo_PECAM1','Myo_ACTA2','Pericytes',
                      'Mal_Proliferation','Mal_Myogenesis','Mal_Stemness') ~ cell_subtype,
    .default = 'others'
  ))|>
  ggplot(aes(x_slide_mm,y_slide_mm))+
  geom_point(aes(color=cell_subtype),shape=16)+
  facet_wrap(~Tumor_subtype,scales = 'free')+
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave('Figure4/insitu_mes_fov.pdf',width = 10,height = 5)


# correlation
t.type = c('DDLPS','MLPS','PLPS')

pdf('Figure4/mes_cor_pl.pdf',width = 8,height = 8)
for (i in t.type) {
  all@meta.data|>
    filter(Tumor_subtype==i)|>
    group_by(fov)|>
    count(cell_subtype)|>
    mutate(prop=n/sum(n))|>
    filter(cell_subtype%in%c('Endo_PECAM1','Myo_ACTA2','Pericytes',
                             'Endo_IGFBP3','Myo_PFN1','Endo_CCL21','CAFs',
                             'Mal_Proliferation','Mal_Myogenesis','Mal_Stemness',
                             'Mal_EMT','Mal_Invasion','Mal_Interferon response'))|>
    ungroup()|>select(fov,cell_subtype,prop)|>
    pivot_wider(names_from = 'cell_subtype',values_from = 'prop',values_fill = 0)|>
    column_to_rownames('fov')|>
    cor()|>
    corrplot::corrplot(title = i,
                       tl.pos = 'l',
                       col = rev(COL2('RdBu', 200)))->p
  p
}
dev.off()


# signature
library(msigdbr)
library(GSVA)
library(UCell)

human = msigdbr(species = "Homo sapiens",category = 'H')
halkmark = human%>%select(gs_name,gene_symbol)%>%split(x = .$gene_symbol, f = .$gs_name)
obj = AddModuleScore_UCell(obj,features = halkmark,ncores = 4,name = '')


# the CV calculation
calculate_cv <- function(x) {
  if (!is.numeric(x)) return(NA)  
  if (all(is.na(x))) return(NA)   
  mean_val <- mean(x, na.rm = TRUE)
  if (mean_val == 0) return(NA)   
  sd(x, na.rm = TRUE) / mean_val
}

top5 = obj@meta.data|>
  group_by(cell_subtype)|>
  select(names(halkmark))|>
  summarise(across(everything(),median))|>
  summarise(across(where(is.numeric), ~ calculate_cv(.))) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "CV") %>%
  filter(!is.na(CV))|>
  arrange(-CV)|>
  head(5)

data_rd=obj@meta.data|>
  group_by(cell_subtype)|>
  select(top5)|>
  summarise(across(everything(),median))|>
  ungroup()

ggradar(data_rd,
        group.point.size = 0)
ggsave('Figure4/radar_mes_hallmark.pdf',width = 8,height = 6)


# ANGIOGENESIS
ANGIOGENESIS = list(ANGIOGENESIS=c(
  'CD93',
  'COL4A1',
  'COL4A2',
  'HSPG2',
  'LAMA4',
  'PXDN',
  'RGCC',
  'VWA1',
  'VWF'
))

obj = AddModuleScore_UCell(obj,features = ANGIOGENESIS,ncores = 4,name = '')
p2=VlnPlot(obj,pt.size = 0,features = 'ANGIOGENESIS',
           idents = c('Myo_PFN1'), 
           group.by = 'Tumor_subtype',adjust = 3)+
  geom_boxplot(width=0.2,fill='white',outliers = F)+
  stat_compare_means(comparisons = list(c('DDLPS','MLPS'),
                                        c('MLPS','PLPS'),
                                        c('DDLPS','PLPS')),
                     label.y = 0.7)

p2
ggsave('Figure4/vln_mes_sig.pdf',width = 6,height = 8)


# neighbor
D.nei = read.csv('squidpy/fov2_neighbor_zsocre.csv')|>
  filter(cell_subtype!='Mesenchymal cells')|>
  select(cell_subtype,Mesenchymal.cells)|>
  mutate(group='DDLPS',
         value=case_when(
           Mesenchymal.cells < 0 ~ -log1p(abs(Mesenchymal.cells)),
           .default = log1p(abs(Mesenchymal.cells))
         ))

M.nei = read.csv('squidpy/fov3_neighbor_zsocre.csv')|>
  filter(cell_subtype!='Mesenchymal cells')|>
  select(cell_subtype,Mesenchymal.cells)|>
  mutate(group='MLPS',
         value=case_when(
           Mesenchymal.cells < 0 ~ -log1p(abs(Mesenchymal.cells)),
           .default = log1p(abs(Mesenchymal.cells))
         ))

P.nei = read.csv('squidpy/fov4_neighbor_zsocre.csv')|>
  filter(cell_subtype!='Mesenchymal cells')|>
  select(cell_subtype,Mesenchymal.cells)|>
  mutate(group='PLPS',
         value=case_when(
           Mesenchymal.cells < 0 ~ -log1p(abs(Mesenchymal.cells)),
           .default = log1p(abs(Mesenchymal.cells))
         ))

df = rbind(D.nei,M.nei,P.nei)

ggplot(df,aes(value,cell_subtype))+
  geom_col()+
  facet_wrap(~group)+
  theme_minimal()+
  labs(x='Z-score',y='',title = 'Neighbor enrichment')
ggsave('Figure4/col_neighbor_mes.pdf',width = 8,height = 6)


# Mal_EMT 
Mes.score = obj@meta.data|>
  group_by(fov)|>
  summarise(across('ANGIOGENESIS',median))

all@meta.data|>
  group_by(Tumor_subtype,fov)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  filter(cell_subtype=='Mal_EMT'&Tumor_subtype%in%c('DDLPS','MLPS','PLPS'))|>
  left_join(Mes.score)|>
  mutate(ANGIOGENESIS=ifelse(is.na(ANGIOGENESIS),0,ANGIOGENESIS))|>
  ggplot(aes(prop,ANGIOGENESIS))+
  geom_point(aes(color=Tumor_subtype),shape=16)+
  geom_smooth(se = F,method = 'lm',formula = 'y~x',lty=2,color='gray')+
  stat_cor(method = 'spearman')+
  facet_wrap(~Tumor_subtype,scales = 'free')+
  theme_classic()+
  theme(aspect.ratio = 2)+
  labs(x='Proprotion of Mal_EMT cells')
ggsave('Figure4/cor_mes_emt.pdf',width = 8,height = 6)


# genes
avg.exp = AverageExpression(obj,assays = 'RNA',slot = 'counts',
                            features = c('PECAM1'),
                            group.by = 'fov')$RNA|>
  as.data.frame()|>t()|>as.data.frame()|>
  rownames_to_column('fov')|>
  mutate(fov=sub('g','',fov)|>as.numeric())

emt = subset(all,subset = cell_subtype == 'Mal_EMT')

col.exp = AverageExpression(emt,assays = 'RNA',slot = 'counts',
                            features = c('COL1A1'),
                            group.by = 'fov')$RNA|>
  as.data.frame()|>t()|>as.data.frame()|>
  rownames_to_column('fov')|>
  mutate(fov=sub('g','',fov)|>as.numeric())

all@meta.data|>
  select(Tumor_subtype,fov)|>
  distinct()|>
  left_join(avg.exp)|>
  mutate(PECAM1=ifelse(is.na(PECAM1),0,PECAM1))|>
  left_join(col.exp)|>
  mutate(COL1A1=ifelse(is.na(COL1A1),0,COL1A1))|>
  filter(Tumor_subtype %in% c('DDLPS','MLPS','PLPS'))|>
  ggplot(aes(COL1A1,PECAM1))+
  geom_point(aes(color=Tumor_subtype),shape=16)+
  geom_smooth(aes(group=Tumor_subtype,color=Tumor_subtype),
              se = F,method = 'lm',formula = 'y~x',lty=2)+
  stat_cor(aes(group=Tumor_subtype),method = 'spearman')+
  facet_wrap(~Tumor_subtype,scales = 'free')+
  theme_classic()+
  theme(aspect.ratio = 1)

# fov show
all@meta.data|>
  filter(fov%in%c(91,104,165))|>
  mutate(cell_class = case_when(
    cell_subtype%in%c('Endo_PECAM1','Mal_EMT')~cell_subtype,
    cell_type%in%c('Lymphocytes','Myeloid cells','Plasma cells')~'Immune',
    cell_type%in%c('Tumor','CAFs')~'Tumor',
    .default = 'Stromal'
  ))|>
  ggplot(aes(x_slide_mm,y_slide_mm))+
  geom_point(aes(color=cell_class),shape=16)+
  facet_wrap(~Tumor_subtype,scales = 'free')+
  theme(aspect.ratio = 1)
ggsave('Figure4/mes_fov_show.pdf',width = 12,height = 9)


all@meta.data|>
  group_by(Tumor_subtype)|>
  count(fov)|>
  filter(Tumor_subtype%in%c('DDLPS','MLPS','PLPS'))|>
  arrange(-n)|>
  group_by(Tumor_subtype)|>
  slice_head(n=5)

all@meta.data|>
  mutate(cell_class = case_when(
    cell_type%in%c('Lymphocytes','Myeloid cells','Plasma cells')~'Immune',
    cell_type%in%c('Tumor','CAFs')~'Tumor',
    .default = 'Stromal'
  ))|>
  ggplot(aes(Tumor_subtype))+
  geom_bar(aes(fill=cell_class),position = 'fill')+
  scale_fill_npg()+
  theme(aspect.ratio = 3)+
  labs(x='',y='porprotion')
ggsave('Figure4/bar_type_prop.pdf',width = 5,height = 6)
