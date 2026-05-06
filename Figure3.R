library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

# Figure 3 ####
obj.tumor = readRDS('obj.tumor.rds')

obj.tumor@meta.data$M.D. = obj.tumor@meta.data$ECM/obj.tumor@meta.data$Fat
p1=obj.tumor@meta.data|>group_by(cell_subtype)|>
  summarise(M.D. = median(M.D.))|>
  ggplot(aes(reorder(cell_subtype,M.D.),y=1))+
  geom_point(aes(color=M.D.),size=10)+
  scale_color_gradient2(low = 'green',mid = 'blue',high = 'red',midpoint = 1)+
  theme_bw()+
  theme(aspect.ratio = 1/5)

obj.tumor@meta.data$cell_subtype = paste('Mal',obj.tumor$cell_subtype,sep = '_')
p2=DimPlot(obj.tumor,group.by = 'cell_subtype',raster = T,pt.size = 1.5)+
  coord_fixed()
p2/p1
ggsave('Figure3/umap_tumor_set.pdf',width = 8,height = 6)


Idents(obj.tumor)<-'cell_subtype'
obj.tumor = JoinLayers(obj.tumor,'RNA')
DEGs = FindAllMarkers(obj.tumor,assay = 'SCT',slot = 'data',
                      logfc.threshold = 0.5,min.pct = 0.25,
                      only.pos = T,return.thresh = 0.05)
features = DEGs|>group_by(cluster)|>
  slice_head(n=3)|>
  pull(gene)
DotPlot(obj.tumor,assay = 'SCT',features = features)+
  coord_flip()+
  theme_bw(base_line_size = 0)+
  scale_color_viridis_c(option = 'H')+
  theme(aspect.ratio = 3,axis.text.x = element_text(angle =90,vjust = 0.5,hjust = 1))
ggsave('Figure3/dot_tumor_set.pdf',width = 6,height = 8)


## SMI signature ####
library(readxl)
SMI_1KPanel_Genes <- read_excel("SMI-1KPanel-Genes.xlsx", sheet = "Annotations", skip = 1)|>
  select(-c('Add-On','Core','Cell Type Associated','Cell Typing'))|>
  column_to_rownames('Gene')

signature_list <- list()
for (col in colnames(SMI_1KPanel_Genes)) {
  index = SMI_1KPanel_Genes[[col]]
  genes = rownames(SMI_1KPanel_Genes)[index == '+']
  signature_list[[col]] <- genes
}


library(UCell)
obj.tumor = AddModuleScore_UCell(obj.tumor,features = signature_list,ncores = 10,
                                 assay = 'RNA',slot = 'counts',name = '')
obj.tumor@meta.data|>select(cell_subtype,names(signature_list))|>
  group_by(cell_subtype)|>
  summarise(across(everything(),median))|>
  summarise(cell_subtype = cell_subtype,
            across(names(signature_list),scale))|>
  pivot_longer(-cell_subtype)|>
  arrange(desc(value))|>
  group_by(cell_subtype)|>
  slice_head(n=5)|>
  pull(name)|>
  unique()->ht_show
obj.tumor@meta.data|>select(cell_subtype,ht_show)|>
  group_by(cell_subtype)|>
  summarise(across(everything(),median))|>
  ungroup()|>column_to_rownames('cell_subtype') -> tumor_ms


pdf('Figure3/ht_tumor_set.pdf',width = 8,height = 6)
pheatmap(mat = tumor_ms,clustering_method = 'ward.D2',
         scale = 'column',cluster_rows = T,
         cellheight = 20,cellwidth = 10,
         treeheight_row = 15,treeheight_col = 15,
         display_numbers = F)
dev.off()


ggplot(obj.tumor@meta.data,aes(Tumor_subtype))+
  geom_bar(aes(fill=cell_subtype),position = 'fill',alpha=0.85)+
  facet_wrap(~Tumor_subtype,scales = 'free_x',nrow = 1)+
  coord_polar('y')+
  theme_void()
ggsave('Figure3/pie_tumor_set.pdf',width = 8,height = 4)


obj.tumor@meta.data|>
  group_by(Tumor_subtype)|>
  count(cell_subtype)|>
  mutate(prop = n/sum(n))|>
  arrange(desc(prop))|>
  slice_head(n=1)
obj.tumor@meta.data|>
  group_by(sample,Tumor_subtype)|>
  select(c('NF-kB Signaling','MAPK Targets','PI3K-Akt Signaling',
           'Insulin Signaling','Androgen Signaling'))|>
  summarise(across(where(is.numeric),mean))|>
  pivot_longer(-c('sample','Tumor_subtype'))|>
  ggplot(aes(name,value))+
  geom_boxplot(aes(fill=Tumor_subtype),staplewidth = 0.8,outlier.size = 0.5)+
  stat_compare_means(aes(group=Tumor_subtype),label = 'p.signif')+
  scale_fill_bmj(alpha=0.5)+
  theme_classic(base_size = 15)+
  labs(x='',y='Mean score of signaling')
ggsave('Figure3/box_tumor_signaling.pdf',width = 8,height = 6)


## genes vlnplot ####
Idents(obj.tumor)<-'cell_subtype'
obj.tumor=NormalizeData(obj.tumor)
p1=VlnPlot(obj.tumor,
           features = c('CDKN1A','H2AZ1','PCNA'),
           idents = 'Proliferation',
           split.by = 'Tumor_subtype',
           pt.size = 0,log = T,slot = 'data',
           assay = 'RNA',adjust = 5)&
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 0))&
  labs(x='',y='')&
  scale_fill_bmj(alpha = 0.5)
p2=VlnPlot(obj.tumor,
           features = c('LGALS1','FN1','ANXA1'),
           idents = 'Myogenesis',
           split.by = 'Tumor_subtype',
           pt.size = 0,log = T,slot = 'data',
           assay = 'RNA',adjust = 5)&
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 0))&
  labs(x='',y='')&
  scale_fill_bmj(alpha = 0.5)
p3=VlnPlot(obj.tumor,
           features = c('POU5F1','TWIST1','CTNNB1'),
           idents = 'Stemness',
           split.by = 'Tumor_subtype',
           pt.size = 0,log = T,slot = 'data',
           assay = 'RNA',adjust = 5)&
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 0))&
  labs(x='',y='')&
  scale_fill_bmj(alpha = 0.5)
p1/p2/p3
ggsave('Figure3/vln_tumor_set.pdf',width = 8,height = 8)


## insitu fov ####
obj.tumor@meta.data|>
  group_by(Tumor_subtype,fov)|>
  count(cell_subtype)|>
  arrange(desc(n))|>
  group_by(Tumor_subtype)|>
  slice_head(n=3)
obj@meta.data|>
  filter(fov %in% c(61,31,104,183,23))|>
  mutate(cell_subtype = ifelse(cell_type=='Tumor',cell_subtype,NA))|>
  ggplot(aes(x_slide_mm,y_slide_mm))+
  geom_point(aes(color=cell_subtype),shape=16,size=1)+
  facet_wrap(~Tumor_subtype,scales = 'free',nrow = 1)+
  theme_dark()+
  theme(aspect.ratio = 1)
ggsave('Figure3/point_spatial_tumor_set.pdf',width = 10,height = 5)


## scRNA ####
sc = readRDS('ref.rds.gz')
sc = subset(sc,subset = nFeature_RNA>200 & nCount_RNA>500)
sc = SCTransform(sc)|>
  RunPCA()
sc = RunHarmony(sc,group.by.vars='patient')
sc = RunUMAP(sc,reduction = 'harmony',dims=1:10,n.neighbors = 50)

Idents(sc)<-'celltype.l1'

DimPlot(sc,raster = T,label = F,pt.size = 1.5)+
  scale_color_d3()+
  coord_equal()
ggsave('Figure3/umap_ref_sc.pdf',width = 8,height = 6)


sig.list = list(Fat = c('FABP4','PPIA','ADIRF','NEAT1','ADIPOQ'),
                ECM = c('COL5A1','COL6A3','COL5A2','COL1A2','COL3A1','COL1A1'))
library(UCell)
sc = AddModuleScore_UCell(sc,features = sig.list,
                          assay = 'SCT',slot = 'data',
                          ncores = 2,name = '')
sc@meta.data$tumor_type = case_when(
  sc@meta.data$type %in% c('DDLPS-DD','DDLPS-DDK','DDLPS-WD','DDLPS-sWD') ~ 'DDLPS',
  sc@meta.data$type == 'lipoma' ~ 'Lipoma',
  sc@meta.data$type == 'nl' ~ 'Adjacent normal',
  .default = 'WDLPS'
)
sc$tumor_type = factor(sc$tumor_type,levels = c('Adjacent normal','Lipoma',
                                                'WDLPS','DDLPS'))
sc$EF = sc$ECM/sc$Fat

VlnPlot(sc,group.by = 'tumor_type',features = 'EF',pt.size = 0,adjust = 3,log = T)+
  geom_boxplot(fill='white',outliers = F,width=.2)+
  stat_compare_means(comparisons = list(c('Adjacent normal','Lipoma'),
                                        c('Lipoma','WDLPS'),
                                        c('WDLPS','DDLPS')))
ggsave('Figure3/vln_ref_sig.pdf',width = 5,height = 6)


## tumor sc ####
obj.tumor = readRDS('obj.tumor.rds')
Idents(obj.tumor)<-'cell_subtype'
obj.tumor = JoinLayers(obj.tumor,'RNA')
DEGs = FindAllMarkers(obj.tumor,assay = 'SCT',slot = 'data',
                      logfc.threshold = 0.5,min.pct = 0.25,
                      only.pos = T,return.thresh = 0.05)
features = DEGs|>group_by(cluster)|>
  slice_head(n=5)|>summarise(genes = list(gene))|>
  deframe()|>as.list()

sc = subset(sc,subset = celltype.l1 == 'Tumor cells')
sc = AddModuleScore_UCell(sc,features = features,
                          assay = 'SCT',slot = 'data',ncores = 2,
                          name = '')
P.genes = DEGs|>filter(cluster=='Proliferation')|>pull(gene)|>unique()

sc = subset(sc,subset= type %in% c("DDLPS-DD","nl","WDLPS","lipoma","DDLPS-DDK"))

Idents(sc) <- 'tumor_type'
DotPlot(sc,features = P.genes,group.by = 'celltype.l2',
        assay = 'SCT',idents = 'DDLPS')+
  coord_flip()+
  scale_color_viridis_c(option = 'H')+
  theme(aspect.ratio = 3,
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1))
ggsave('Figure3/dot_ref_marker.pdf',width = 6,height = 6)


pdf('Figure3/ht_ref_sm.pdf',width = 6,height = 6)
sc@meta.data|>
  filter(tumor_type=='DDLPS')|>
  group_by(celltype.l2)|>
  summarise(across(names(features),mean))|>
  column_to_rownames('celltype.l2')|>
  pheatmap(scale = 'column',
           color = colorRampPalette(c('#60a9e1','white','#d95e6e'))(100),
           cellwidth = 30,cellheight = 30)
dev.off()

sc@meta.data|>
  filter(tumor_type=='DDLPS')|>
  ggplot(aes(tumor_type))+
  geom_bar(aes(fill=max_col),position = 'fill')


DotPlot(sc,features = names(features),group.by = 'tumor_type')+
  coord_flip()

sc_type = sc@meta.data%>%
  select(names(features))%>%
  mutate(
    max_col = names(.)[max.col(.,'first')]
  )
sc$max_col = sc_type$max_col

sc@meta.data|>
  group_by(tumor_type)|>
  count(max_col)|>
  mutate(prop =n/sum(n))|>
  filter(tumor_type!='Adjacent normal')|>
  ggplot(aes(tumor_type,prop))+
  geom_col(aes(fill=max_col),position = 'dodge')+
  facet_wrap(~max_col,scales = 'free_y',nrow = 2)
ggsave('Figure3/col_ref_prop.pdf',width = 8,height = 6)

## ref degs ####
Idents(sc)<-'celltype.l2'
sc = JoinLayers(sc,'RNA')
sc = NormalizeData(sc,assay = 'RNA')
DEGs = FindMarkers(sc,ident.1 = 'WDLPS',ident.2 = 'DDLPS',
                   assay = 'RNA',slot='counts',
                   group.by = 'tumor_type',
                   subset.ident = 'Invasion')

library(ggrepel)
DEGs|>
  rownames_to_column('gene')|>
  mutate(signif = case_when(
    avg_log2FC>1 & pct.1>0.25 & (pct.1-pct.2)>0.2 ~ 'DDLPS',
    avg_log2FC< -1 & pct.2>0.25 & (pct.2-pct.1)>0.2 ~ 'WDLPS'
  ))|>
  mutate(gene = ifelse(is.na(signif),NA,gene))|>
  ggplot(aes(pct.1,pct.2))+
  geom_point(aes(color=signif),shape=16)+
  geom_text_repel(aes(label=gene))+
  theme_classic()+
  coord_equal()
ggsave('Figure3/point_ref_degs.pdf',width = 8,height = 8)

DDLPS.gene = DEGs|>
  filter(avg_log2FC>1 & pct.1>0.25 & (pct.1-pct.2)>0.2)|>
  rownames()


exp = GetAssayData(sc,'RNA','data')

sc = AddMetaData(sc,exp['MDM4',],'MDM4')
sc = AddMetaData(sc,exp['CDK4',],'CDK4')
sc = AddMetaData(sc,exp['NR4A1',],'NR4A1')
sc = AddMetaData(sc,exp['MCL1',],'MCL1')


p1=sc@meta.data|>
  filter(celltype.l2=='Invasion')|>
  ggplot(aes(MDM4,NR4A1))+
  geom_point(aes(color=tumor_type))+
  stat_cor()+
  theme(aspect.ratio = 1)

p2=sc@meta.data|>
  filter(celltype.l2=='Invasion')|>
  ggplot(aes(MDM4,MCL1))+
  geom_point(aes(color=tumor_type))+
  stat_cor()+
  theme(aspect.ratio = 1)

p1+p2
ggsave('Figure3/point_ref_MDM4.pdf',width = 10,height = 5)
