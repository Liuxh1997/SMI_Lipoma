library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(GeneNMF)
library(ggpubr)

# Figure 6 ####
cci.l = readRDS('CCI/Lipoma_cell_type_CCI.rds')
cci.w = readRDS('CCI/WDLPS_cell_type_CCI.rds')
cci.d = readRDS('CCI/DDLPS_cell_type_CCI.rds')
cci.m = readRDS('CCI/MLPS_cell_type_CCI.rds')
cci.p = readRDS('CCI/PLPS_cell_type_CCI.rds')

cellchat <- mergeCellChat(list(cci.w,cci.d), add.names = c('WDLPS','DDLPS'), cell.prefix = TRUE)
## heatmap
par(mfrow = c(1,1))
pdf('Figure6/cci_ht.pdf',width = 6,height = 6)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

## signal
rankNet(cellchat, mode = "comparison",
        stacked = T,do.stat = TRUE,
        measure = 'weight',slot.name = 'netP',
        sources.use = 'CAFs',targets.use = 'Tumor')
ggsave('Figure6/cci_bar.pdf',width = 6,height = 6)

## nichenet
devtools::install_github("saeyslab/nichenetr")
library(nichenetr)

ligand_target_matrix = readRDS('../YKKY0060_hejinhua_TNBC/nichenet/ligand_target_matrix.rds')
lr_network = readRDS('../YKKY0060_hejinhua_TNBC/nichenet/lr_network.rds')
weighted_networks = readRDS('../YKKY0060_hejinhua_TNBC/nichenet/weighted_networks.rds')

s.set = subset(obj,subset = Tumor_subtype %in% c('WDLPS','DDLPS'))

DefaultAssay(s.set)<-'RNA'
s.set = JoinLayers(s.set,assay = 'RNA')
Idents(s.set)<-'cell_type'
NE=nichenet_seuratobj_aggregate(seurat_obj = s.set,
                                receiver = 'Tumor',sender = 'CAFs',
                                condition_colname = 'Tumor_subtype',
                                condition_oi = 'DDLPS',condition_reference = 'WDLPS',
                                ligand_target_matrix = ligand_target_matrix,
                                lr_network = lr_network,
                                weighted_networks = weighted_networks)

NE$ligand_activity_target_heatmap
ggsave('Figure6/cci_lr_activity.pdf',width = 15,height = 5)

## upset
install.packages("UpSetR")
library(UpSetR)

listInput=list()
for (i in unique(obj$Tumor_subtype)) {
  tmp <- read_csv(paste0("CCI/",i,"_cell_type_LRs.csv"))
  tmp <- filter(tmp,source == 'Tumor' | target == 'Tumor')
  listInput[[i]] <- unique(tmp$pathway_name)
}

pdf('Figure6/cci_upset.pdf',width = 8,height = 6)
upset(fromList(listInput),nsets = 50, order.by = "freq")
dev.off()


max_length <- max(sapply(listInput, length))
filled_list <- lapply(listInput, function(x) {
  c(x, rep(NA, max_length - length(x)))
})
df <- as.data.frame(filled_list, stringsAsFactors = FALSE)
write.csv(df,'cci_tumor_subtype.csv',row.names = F)


pathway = c(
  'NRG',
  'BAFF',
  'KIT',
  'VTN',
  'PECAM1',
  'ncWNT',
  'GALECTIN',
  'NCAM'
)

pathway = c(
  'TGFb',
  'ACTIVIN',
  'PDGF',
  'CCL',
  'IL6',
  'IL1',
  'CSF',
  'LIGHT',
  'NPR2',
  'COLLAGEN',
  'CDH1',
  'EPHA',
  'EPHB',
  'MHC-II'
)

p.list = list()
for (i in c('WDLPS','DDLPS','MLPS','PLPS')) {
  tmp <- read_csv(paste0("CCI/",i,"_cell_type_LRs.csv"))
  p=tmp|>group_by(pathway_name)|>
    summarise(across(prob,sum))|>
    filter(pathway_name%in%pathway)|>
    arrange(-prob)|>
    ggplot(aes(log10(prob*1e12),reorder(pathway_name,prob)))+
    geom_col()+
    labs(title = i)
  p.list[[i]]=p
}

library(cowplot)
plot_grid(plotlist = p.list,nrow = 1)
ggsave('Figure6/cci_signal_bar.pdf',width = 8,height = 6)


## signal show 
par(mfrow=c(1,4))
netVisual_aggregate(cci.w, signaling = 'COLLAGEN', layout = "circle",signaling.name = 'WDLPS')
netVisual_aggregate(cci.d, signaling = 'COLLAGEN', layout = "circle",signaling.name = 'DDLPS')
netVisual_aggregate(cci.m, signaling = 'COLLAGEN', layout = "circle",signaling.name = 'MLPS')
netVisual_aggregate(cci.p, signaling = 'COLLAGEN', layout = "circle",signaling.name = 'PLPS')


## collagen
collagen = cci.d@LR[["LRsig"]]|>filter(pathway_name=='COLLAGEN')|>pull(ligand)|>unique()
recepter = c('ITGA1','ITGA2','ITGA3','ITGB1','ITGB8','ITGAV','CD44')

DSP.exp|>t()|>as.data.frame()|>select(any_of(c(collagen)))|>
  cbind(DSP.anno)|>filter(Morphology=='TUMOR')|>
  mutate(tissue_type = factor(tissue_type,levels = c('Lipoma','WDLPS','DDLPS','MLPS','PLPS')))|>
  pivot_longer(c(any_of(collagen)))|>
  ggplot(aes(name,log1p(value)))+
  geom_boxplot(aes(fill=tissue_type))+
  stat_compare_means(aes(group = tissue_type),label = 'p.signif')+
  theme_bw(base_size = 15)+
  labs(x='',y='log1p(Exp)',title = 'TUMOR Morphology')+
  theme(aspect.ratio = 1/2,
        axis.text.x = element_text(angle = 90,hjust = 0,vjust = 0.5))
ggsave('Figure6/dsp_collagen_boxpl.pdf',width = 8,height = 6)
