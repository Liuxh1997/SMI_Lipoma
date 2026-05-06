library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

# Figure 5 ####
obj = readRDS('obj.raw.rds')
DimPlot(obj,group.by = 'cell_type')

IO.set = subset(obj, subset = cell_type %in% c("Plasma cells","Myeloid cells","Lymphocytes"))

Idents(IO.set)<-'cell_subtype'
IO.set = PrepSCTFindMarkers(IO.set)
DEGs = FindAllMarkers(IO.set,assay = 'SCT',slot = 'counts',
                      logfc.threshold = 0.5,min.pct = 0.1,
                      only.pos = T,return.thresh = 0.05)
write.csv(DEGs,'Supplementary/TableS_IO_DEGs.csv',row.names = F)

IO.set = RunPCA(IO.set,features = unique(DEGs$gene))
ElbowPlot(IO.set,50,'pca')
IO.set = RunHarmony(IO.set,group.by.vars='sample')
ElbowPlot(IO.set,50,'harmony')
IO.set = RunUMAP(IO.set,reduction = 'harmony',dims = 1:10,
                 umap.method = 'umap-learn',metric = 'correlation',
                 min.dist = 1e-2,n.neighbors = 35)

DimPlot(IO.set,label = T,raster = T,repel = T,pt.size = 1)+
  scale_color_d3(palette = 'category20c')
ggsave('Figure5/umap_io_celltype.pdf',width = 8,height = 6)


IO.set@meta.data|>
  group_by(Tumor_subtype)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  ggplot(aes(Tumor_subtype,prop))+
  geom_col(aes(fill=cell_subtype))+
  facet_wrap(~Tumor_subtype,scales = 'free',nrow = 1)+
  coord_polar('y')+
  scale_fill_d3(palette = 'category20')+
  theme(legend.position = 'top',
        axis.text.y = element_blank())+
  labs(y='')
ggsave('Supplementary/FigS5_subcluster_pie.pdf',width = 8,height = 6)


features = DEGs|>group_by(cluster)|>slice_head(n=3)|>pull(gene)|>unique()

DotPlot(IO.set,features = features[c(1:18,21,24:27)])+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        aspect.ratio = 1.5)+
  scale_color_viridis_c()
ggsave('Figure5/dot_io_marker.pdf',width = 6,height = 8)


IO.set@meta.data|>
  group_by(fov,Tumor_subtype)|>
  count(cell_type)|>
  mutate(prop=n/sum(n))|>
  ggplot(aes(Tumor_subtype,prop))+
  geom_point(aes(color=Tumor_subtype),alpha=0.5,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(fill=Tumor_subtype),outliers = F,alpha=0.2)+
  stat_compare_means(ref.group = 'Lipoma',
                     label = 'p.signif',
                     method = 't.test')+
  facet_wrap(~cell_type,ncol = 2)+
  theme_bw()+
  theme(aspect.ratio = 1)
ggsave('Figure5/denstiny_io_boxpl.pdf',width = 8,height = 6)


## sub_type density ####
IO.set@meta.data|>
  group_by(fov,Tumor_subtype)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  ggplot(aes(Tumor_subtype,prop))+
  geom_point(aes(color=Tumor_subtype),alpha=0.5,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(fill=Tumor_subtype),outliers = F,alpha=0.2)+
  stat_compare_means(ref.group = 'Lipoma',
                     label = 'p.signif',
                     method = 't.test')+
  facet_wrap(~cell_subtype,scales = 'free')+
  theme_bw()+
  theme(aspect.ratio = 1)

## correlation ####
obj@meta.data|>
  group_by(fov,Tumor_subtype)|>
  count(cell_subtype)|>
  mutate(prop=n/sum(n))|>
  select(-n)|>
  pivot_wider(names_from = 'cell_subtype',values_from = 'prop',values_fill = 0)|>
  filter(Tumor_subtype=='DDLPS')|>ungroup()|>
  select(-c(fov,Tumor_subtype))|>
  cor()|>corrplot::corrplot()

## signature ####
library(UCell)
IO.set = AddModuleScore_UCell(IO.set,features = signature_list,assay = 'RNA',
                              ncores = 4,name = '')

imm.sig = c(
  'Antigen Presentation',
  'Apoptosis',
  'BCR Signaling',
  'Cellular Stress',
  'Cytotoxicity',
  'Cytokines & Chemokines',
  'Inflammation',
  'Interferon Response Genes',
  'Pattern Recognition Receptors',
  'Proteases',
  'T Cell Checkpoints',
  'T Cell Exhaustion',
  'TCR Signaling',
  'TGF-beta Signaling',
  'TNF Signaling'
)

IO.set@meta.data = IO.set@meta.data%>%
  mutate(across(all_of(imm.sig),~ ifelse(.<quantile(.,0.25,na.rm = T),0,.)))

Idents(IO.set)<-'Tumor_subtype'

p.list=list()
for (i in unique(IO.set$Tumor_subtype)) {
  p=DotPlot(IO.set,features = imm.sig,idents = i,group.by = 'cell_subtype')+
    scale_color_gradient2(low = '#5766e6',high = '#dc2a40')+
    theme_bw()+
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
    labs(title = i)+
    scale_size(breaks = seq(10,90,20),
               range = c(1,5))
  p.list[[i]]=p
}

library(cowplot)

plot_grid(plotlist = p.list,nrow = 1,scale = T)
ggsave('Figure5/dot_io_sig.pdf',width = 30,height = 5,limitsize = F)

## vln plot ####
Idents(IO.set)<-'cell_subtype'
IO.set@meta.data|>
  filter(cell_subtype=='CD8T' & `T Cell Exhaustion` > 0.5)|>
  group_by(fov,Tumor_subtype)|>
  summarise(across(c('T Cell Exhaustion'),mean))|>
  ggplot(aes(Tumor_subtype,`T Cell Exhaustion`,color=Tumor_subtype))+
  geom_violin(adjust=2,scale = 'width',width=0.6)+
  geom_boxplot(width=0.2,outliers = F)+
  geom_point(position = position_jitter(width = 0.2))+
  stat_compare_means()+
  theme_classic()+
  theme(aspect.ratio = 1)
ggsave('Figure5/CD8T_stat_vln.pdf',width = 8,height = 6)

IO.set@meta.data|>
  filter(cell_subtype=='CD4T' & `T Cell Checkpoints` > 0.5)|>
  group_by(fov,Tumor_subtype)|>
  summarise(across(c('T Cell Checkpoints'),mean))|>
  ggplot(aes(Tumor_subtype,`T Cell Checkpoints`,color=Tumor_subtype))+
  geom_violin(adjust=2,scale = 'width',width=0.6)+
  geom_boxplot(width=0.2,outliers = F)+
  geom_point(position = position_jitter(width = 0.2))+
  stat_compare_means()+
  theme_classic()+
  theme(aspect.ratio = 1)
ggsave('Figure5/CD4T_stat_vln.pdf',width = 8,height = 6)

## cell num ####
library(dplyr)
library(purrr)

find_nearest_cells <- function(df, target_cell_type = "A", k = 10) {
  result <- df %>%
    group_by(sample, Tumor_subtype) %>%
    group_modify(~ {
      target_cells <- .x %>% filter(cell_subtype == target_cell_type)
      if (nrow(target_cells) == 0) {
        return(data.frame(
          cell_subtype = character(),
          count = integer()
        ))
      }
      distances <- map_dfr(1:nrow(target_cells), function(i) {
        target_x <- target_cells$x_slide_mm[i]
        target_y <- target_cells$y_slide_mm[i]
        .x %>%
          mutate(
            distance = sqrt((x_slide_mm - target_x)^2 + (y_slide_mm - target_y)^2),
            target_cell_id = i
          )
      })
      nearest_cells <- distances %>%
        group_by(target_cell_id) %>%
        arrange(distance) %>%
        slice_head(n = k) %>%
        ungroup()
      cell_counts <- nearest_cells %>%
        count(cell_subtype, name = "count")
      return(cell_counts)
    }) %>%
    ungroup()
  return(result)
}

CD8.n = find_nearest_cells(df = obj@meta.data,target_cell_type = 'CD8T',k = 10)
CD4.n = find_nearest_cells(df = obj@meta.data,target_cell_type = 'CD4T',k = 10)

p.list = list()
n=1
for (i in list(CD8.n,CD4.n)) {
  p=i|>group_by(Tumor_subtype,cell_subtype)|>
    summarise(across(count,sum))|>
    arrange(-count)|>
    slice_head(n=5)|>
    ggplot(aes(cell_subtype,Tumor_subtype))+
    geom_point(aes(size=log1p(count),color=log1p(count)),shape=15)+
    scale_size(range = c(1,8))+
    scale_color_viridis_c(option = 'A')+
    theme(aspect.ratio = 1/3,axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
    labs(x='',y='')
  p.list[[n]]=p
  n=n+1
}

plot_grid(plotlist = p.list,scale = T,nrow = 2)
ggsave('Figure5/dot_io_cd4_cd8.pdf',width = 8,height = 6)

## 

obj@meta.data|>
  group_by(Tumor_subtype,sample)|>
  summarise(
    fov_count = n_distinct(fov),
    mac_163_sum = sum(cell_subtype=='Mac_CD163'),
    avg_num = mac_163_sum/fov_count,
    .groups = 'drop'
  )|>group_by(Tumor_subtype)|>
  summarise(
    mean_mac_fov = mean(avg_num),
    var_mac_fov = sd(avg_num),
    .groups = 'drop'
  )|>
  ggplot(aes(Tumor_subtype,mean_mac_fov))+
  geom_col(aes(fill=Tumor_subtype),width = 0.7)+
  geom_errorbar(aes(ymin=mean_mac_fov,ymax = mean_mac_fov+var_mac_fov,
                    color=Tumor_subtype))+
  coord_flip()+
  theme(aspect.ratio = 1.5)
ggsave('Figure5/col_mac_163.pdf',width = 6,height = 6)


## DSP ####
library(readxl)
library(ComplexHeatmap)

DSP.exp <- read_excel("DSP/normalized_Q3.xlsx",sheet = "TargetCountMatrix")|>
  column_to_rownames('TargetName')
DSP.anno <- read_excel("DSP/normalized_Q3.xlsx",sheet = "SegmentProperties")|>
  column_to_rownames('SegmentDisplayName')|>
  mutate(group = paste(tissue_type,Morphology,sep = '_'))

Heatmap(t(scale(t(DSP.exp))),
        show_row_names = F,
        show_column_names = F,column_title = NA,
        column_split = DSP.anno$group,cluster_column_slices = F,
        top_annotation = columnAnnotation(tumor = DSP.anno$tissue_type,
                                          seg = DSP.anno$Morphology),
        use_raster = T)

## MDM2 CDK4 ####
tumor.id = DSP.anno|>filter(Morphology=='TUMOR')|>rownames()
tumor.anno = DSP.anno|>filter(Morphology=='TUMOR')|>
  rownames_to_column('id')

tumor.exp = DSP.exp|>select(tumor.id)|>
  rownames_to_column('genes')|>
  filter(genes %in% c('CDK4','MDM2'))|>
  column_to_rownames('genes')|>
  t()|>
  as.data.frame()|>
  rownames_to_column('id')

left_join(tumor.anno,tumor.exp)|>
  pivot_longer(c('CDK4','MDM2'))|>
  mutate(tissue_type = factor(tissue_type,levels=c('Lipoma','WDLPS','DDLPS','MLPS','PLPS')))|>
  ggplot(aes(tissue_type,log1p(value),color=tissue_type))+
  geom_boxplot(staplewidth = 0.6,width=0.6,outliers = F)+
  geom_point(position = position_jitter(width = 0.2),shape=16)+
  facet_wrap(~name,ncol = 1,strip.position = 'left')+
  stat_compare_means()+
  theme_gray()+
  theme(aspect.ratio = 1.5)
ggsave('Figure5/dsp_marker_pl.pdf',width = 12,height = 9)

## signature ECM/FAT ####
library(UCell)
DSP = CreateSeuratObject(counts = DSP.exp,meta.data = DSP.anno)
DSP$tissue_type = factor(DSP$tissue_type,levels=c('Lipoma','WDLPS','DDLPS','MLPS','PLPS'))

sig.list = list('ECM'=c('S100A6','MRC2','COL6A1','FN1','VIM','HSPB1','CD63'),
                'Fat'=c('FABP4','PPIA','ADIRF','NEAT1','ADIPOQ'))

DSP = AddModuleScore_UCell(DSP,sig.list,ncores = 4,name = '')

DSP@meta.data|>
  pivot_longer(c('ECM','Fat'))|>
  filter(Morphology=='TUMOR')|>
  ggplot(aes(tissue_type,value,color=tissue_type))+
  geom_boxplot(staplewidth = 0.6,width=0.6,outliers = F)+
  geom_point(position = position_jitter(width = 0.2),shape=16)+
  facet_wrap(~name,ncol = 1,strip.position = 'left')+
  stat_compare_means()+
  theme_gray()+
  theme(aspect.ratio = 1.5)
ggsave('Figure5/dsp_sig_pl.pdf',width = 12,height = 9)

## ImmuneAI cell ####
ImmCellAI_hsa_result <- read.csv("DSP/result/results/6.Signature/ImmCellAI/ImmCellAI_hsa_result.csv",
                                 row.names = 1)
rownames(ImmCellAI_hsa_result) = gsub('_',' | ',rownames(ImmCellAI_hsa_result))

CD45.id = DSP.anno|>filter(Morphology=='IM')|>rownames()
CD45.cell = ImmCellAI_hsa_result[CD45.id,]

DSP.anno|>filter(Morphology=='IM')|>
  mutate(tissue_type = factor(tissue_type,levels=c('WDLPS','DDLPS','MLPS','PLPS')))|>
  cbind(CD45.cell)|>
  ggplot(aes(tissue_type,Exhausted,color=tissue_type))+
  geom_boxplot(width=0.6,outliers = F)+
  geom_point(position = position_jitter(width = 0.2))+
  stat_compare_means()+
  theme(aspect.ratio = 1)
ggsave('Figure5/dsp_cd45_exh.pdf',width = 6,height = 6)

Idents(DSP)<-'Morphology'
VlnPlot(DSP,features = 'CD163',group.by = 'tissue_type',idents = 'IM',pt.size = 2,adjust = 3)+
  geom_boxplot(width=0.2,fill='white',outliers = F)+
  theme(aspect.ratio = 1)
ggsave('Figure5/dsp_cd163_exp.pdf',width = 6,height = 6)

DSP.anno|>
  mutate(CD163 = DSP.exp['CD163',]|>as.numeric(),
         Exhausted = ImmCellAI_hsa_result$Exhausted,
         tissue_type = factor(tissue_type,levels=c('Lipoma','WDLPS','DDLPS','MLPS','PLPS')))|>
  ggplot(aes(log1p(CD163),Exhausted,color=tissue_type))+
  geom_point()+
  stat_cor()+
  geom_smooth(method = 'lm',formula = 'y~x',se = F)+
  theme_classic()+
  theme(aspect.ratio = 1/2)
ggsave('Figure5/dsp_cor.pdf',width = 8,height = 6)
