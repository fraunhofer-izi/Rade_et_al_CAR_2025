print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Figure 2")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "scCustomize", "ggh4x", "ggrepel", "anndata", "scico", "survival"
)
.bioc_packages = c("dittoSeq", "speckle", "limma")

# Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

if (any(!"SeuratWrappers" %in% installed.packages())) {
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SeuratData" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"ProjecTILs" %in% installed.packages())) {
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

if (any(!"SeuratDisk" %in% installed.packages())) {
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")
source("code/helper/adt_rna_gene_mapping.R")
theme_set(mytheme(base_size = 8))

cust.th = theme(
  legend.spacing.y = unit(2, 'mm'),
  legend.key.size = unit(4, "mm"),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

base.size = 15

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(
  paste0(manifest$meta_pub, "integration/05_seurat_harmony_all_new.Rds")
)

se.meta$PATIENT_ID_SHORT = gsub("Patient0", "P", se.meta$SAMPLE_ID)
se.meta$PATIENT_ID_SHORT = gsub("_1$", "", se.meta$PATIENT_ID_SHORT)
se.meta@meta.data = droplevels(se.meta@meta.data)

se.meta$celltype_short_3 = ifelse(
  as.character(se.meta$celltype_short_3) == "Progenitor",
  as.character(se.meta$celltype), as.character(se.meta$celltype_short_3)
)
se.meta$celltype_short_3 = factor(se.meta$celltype_short_3)

se.t = readRDS(
  paste0(manifest$meta_pub$work, "integration/06_seurat_harmony_t_all_new.Rds")
)
se.t@meta.data = droplevels(se.t@meta.data)

pd = se.meta@meta.data
min(table(pd$orig.ident))
max(table(pd$orig.ident))
median(table(pd$orig.ident))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Number of cell per cell type
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.meta@meta.data
df = reshape2::melt(table(pd$celltype, pd$TIMEPOINT))

ggsave2(
  filename="publication/extended_data_files/Fig_2_nbr_cells_per_celltypes.jpg",
  ggplot(df, aes(Var1, value)) +
    geom_bar(stat="identity", width=0.8) +
    scale_y_log10(
    ) +
    facet_wrap(~ Var2, nrow = 3) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
      legend.title = element_text(margin = margin(b = 3, unit = "pt"), face = "plain")
    ) +
    ylab("Number of cells") + xlab(NULL),
  width = 180, height = 65, dpi = 200,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ct.pl =
  dimreduc_celltypes(
    obj = se.meta,
    ncol3 = 1,
    ncol2 = 1,
    base.size = base.size,
    leg.size = 6,
    raster = T,
    raster.pt.size = 1,
    dim = "umap",
    col.cc = F
  ) +
  theme(
    legend.spacing.y = unit(0, 'mm'),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    axis.title = element_blank(),
    legend.title = element_text(margin = margin(b = 3, unit = "pt"), face = "plain"),
    legend.text = element_text(
      margin = margin(l = 2, t=3, b=3, unit = "pt"), size = rel(.9)
    ),
  )

###

pd = get_metadata(se.meta)
q = quantile(pd$CellCycle_SCORE_UCell, .9999)
pd$CellCycle_SCORE_UCell[pd$CellCycle_SCORE_UCell > q] = q

cc.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = CellCycle_SCORE_UCell)) +
  scattermore::geom_scattermore(pointsize = 6, color="black")+
  scattermore::geom_scattermore(pointsize = 4, color="white") +
  # scattermore::geom_scattermore(pointsize = 1.5, pixels = c(512,512)) +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  scale_color_gradientn(
    colors = c(scico::scico(30, palette = "navia", direction = -1)),
    na.value = "white", limits = c(0.05, max(pd$CellCycle_SCORE_UCell)),
    breaks = scales::pretty_breaks(3)
  ) +
  theme(
    legend.text = element_text(size = rel(.8)),
    legend.spacing.y = unit(2, 'mm'),
    legend.spacing.x = unit(5, 'mm'),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.margin = margin(t=-5)
  ) +
  guides(
    color = guide_colorbar(
      title = "Cycling score", title.vjust = 1.1,
      barwidth = unit(7.5, 'lines'), barheight = unit(.6, 'lines'),
      ticks.linewidth = 2/.pt, ticks = T, frame.colour="black",
      frame.linewidth = 0.5/.pt
    )
  )

###

pd = get_metadata(se.meta)
q = quantile(pd$Perc_of_mito_genes, .999)
pd$Perc_of_mito_genes[pd$Perc_of_mito_genes > q] = q
mt.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = Perc_of_mito_genes)) +
  scattermore::geom_scattermore(pointsize = 6, color="black")+
  scattermore::geom_scattermore(pointsize = 4, color="white") +
  # scattermore::geom_scattermore(pointsize = 1.5, pixels = c(512,512)) +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  scale_color_gradientn(
    colors = c(scico::scico(30, palette = "navia", direction = -1)),
    breaks = scales::pretty_breaks(3)
  ) +
  theme(
    legend.text = element_text(size = rel(.8)),
    legend.spacing.y = unit(2, 'mm'),
    legend.spacing.x = unit(5, 'mm'),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.margin = margin(t=-5)
  ) +
  guides(
    color = guide_colorbar(
      title = "% Mito", title.vjust = 1.1,
      barwidth = unit(7.5, 'lines'), barheight = unit(.6, 'lines'),
      ticks.linewidth = 2/.pt, ticks = T, frame.colour="black",
      frame.linewidth = 0.5/.pt
    )
  )


###

pd = get_metadata(se.meta)
pd = pd %>%
  mutate(
    TMP = case_when(
      VDJ_T_AVAIL == T ~ "TCR",
      VDJ_B_AVAIL == T ~ "BCR"
    )
  )
vdj.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = TMP)) +
  scattermore::geom_scattermore(pointsize = 6, color="black")+
  scattermore::geom_scattermore(pointsize = 4, color="white") +
  # scattermore::geom_scattermore(pointsize = 1.5, pixels = c(1024,1024)) +
  ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
  mytheme(base_size = base.size) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t=-5, b = 9),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, rnow = 1, override.aes = list(shape = 16, size = 6)
  )) +
  scale_color_manual(
    values = c("#CCBB44", "#CC6677"), breaks = c("TCR", "BCR"),
    na.value = "#FFFFFF"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc: RNA marker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs_to_plot=c(
  'CD4', 'CD8A', 'KLRF1', "CD34", "TNFRSF17", 'CD19', "LILRA4", 'CD14',
  'FCGR3A',	"FCER1A"
)
names(ftrs_to_plot) = ftrs_to_plot
names(ftrs_to_plot)[names(ftrs_to_plot) == "TNFRSF17"] = "BCMA"
names(ftrs_to_plot)[names(ftrs_to_plot) == "FCGR3A"] = "CD16"
# FeaturePlot_scCustom(se.meta, features = "KLRF1", na_cutoff = 0.25)

global.max = quantile(unlist(FetchData(se.meta, ftrs_to_plot)), .99)

ftrs.l = list()
pd = get_metadata(se.meta)

for (i in 1:length(ftrs_to_plot)) {

  pd.pl = pd
  # i = 8
  exprs = FetchData(se.meta, vars = ftrs_to_plot[i], layer = "data")
  pd.pl$EXPRS = exprs[[1]]
  pd.pl$EXPRS[pd.pl$EXPRS >= global.max] = global.max
  c = c(scico::scico(30, palette = "navia", direction = -1))

  if(ftrs_to_plot[i] == "TNFRSF17"){
    pd.pl = pd.pl[order(pd.pl$EXPRS, decreasing = F), ]
  }

  pl =
    ggplot(pd.pl, aes(x = umap_1, y = umap_2, color = EXPRS)) +
    scattermore::geom_scattermore(pointsize = 6.5, color="black")+
    scattermore::geom_scattermore(pointsize = 4, color="white") +
    ggrastr::geom_jitter_rast(shape = ".", raster.dpi = 300, scale = .5) +
    # geom_point(shape = ".") +
    # scattermore::geom_scattermore(pointsize = 1.5, pixels = c(512,512)) +
    scale_color_gradientn(
      colors =c[1:(length(c)-4)],
      na.value = "white", limits = c(0, global.max),
      breaks = scales::pretty_breaks(3)
    ) +
    mytheme(base_size = base.size) +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title  = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(.8)),
      # legend.box.margin=margin(0,0,0,-10),
      plot.title = element_text(hjust = 0.5, colour = "black", size = rel(1))
    ) +
    ggtitle(names(ftrs_to_plot[i])) +
    guides(
      color = guide_colorbar(
        title = "RNA Exprs ", title.vjust = 1.1,
        barwidth = unit(7.5, 'lines'), barheight = unit(.6, 'lines'),
        ticks.linewidth = 2/.pt, ticks = T, frame.colour="black",
        frame.linewidth = 0.5/.pt
      )
    )
  ftrs.l[[i]] = pl
}

legend = get_legend(ftrs.l[[1]])
ftrs.l = lapply(ftrs.l, function(x){
  x = x + theme(legend.position='none', panel.border = element_blank())
})

ftrs.pl = plot_grid(plotlist = ftrs.l, nrow = 2, scale = .98)
ftrs.pl = plot_grid(ftrs.pl, ggdraw(legend), rel_heights = c(1, .15), nrow = 2)
# ftrs.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Barplot: sample overview
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
tmp = naturalsort(rowSums(table(p.data$PATIENT_ID_SHORT, p.data$TIMEPOINT) > 0), decreasing = T)
p.data = p.data[!duplicated(p.data$PATIENT_ID_SHORT), ]
p.data$TMP = tmp[match(p.data$PATIENT_ID_SHORT, names(tmp))]
p.data = p.data %>% dplyr::mutate(
  BEST_RESPONSE = dplyr::case_when(
    BEST_RESPONSE == "VGPR" | BEST_RESPONSE == "PR" ~ 'VGPR/PR',
    BEST_RESPONSE == "SD" | BEST_RESPONSE == "PD" ~ 'SD/PD',
    TRUE ~ BEST_RESPONSE
  )
)
p.data$BEST_RESPONSE = factor(p.data$BEST_RESPONSE, levels = c("CR", "VGPR/PR", "SD/PD"))

lp = unique(se.meta@meta.data[se.meta@meta.data$TIMEPOINT == "LP", ]$PATIENT_ID_SHORT)
l = unique(se.meta@meta.data[se.meta@meta.data$TIMEPOINT == "Late", ]$PATIENT_ID_SHORT)
vl = unique(se.meta@meta.data[se.meta@meta.data$TIMEPOINT == "Very Late", ]$PATIENT_ID_SHORT)
p.data = p.data %>% dplyr::mutate(
  TP_LP = dplyr::case_when(PATIENT_ID_SHORT %in% lp ~ "yes"),
  TP_LATE = dplyr::case_when(PATIENT_ID_SHORT %in% l ~ "yes"),
  TP_V_LATE = dplyr::case_when(PATIENT_ID_SHORT %in% vl ~ "yes")
)

p.data = p.data %>% dplyr::arrange(BEST_RESPONSE, PRODUCT, CRS_GRADE, SEX, desc(TMP))
# dplyr::select(TMP, PATIENT_ID_SHORT, PRODUCT, RESPONSE_OVERALL)

p.data.sub = p.data %>%
  dplyr::select(PATIENT_ID_SHORT, SEX, PRODUCT, CRS_GRADE, BEST_RESPONSE, TP_LP, TP_LATE, TP_V_LATE) %>%
  dplyr::mutate_if(is.factor, as.character)
lvls = p.data.sub$PATIENT_ID_SHORT
p.data.sub = reshape2::melt(p.data.sub, id = "PATIENT_ID_SHORT")
p.data.sub$PATIENT_ID_SHORT = factor(p.data.sub$PATIENT_ID_SHORT, levels = rev(lvls))

x.labels = setNames(
  c("Response", "Product", "CRS Gr", "LP", "Late", "Very Late"),
  c("BEST_RESPONSE", "PRODUCT", "CRS_GRADE", "TP_LP", "TP_LATE", "TP_V_LATE")
)

# Set levels for response
p.data.sub.resp = subset(p.data.sub, variable == "BEST_RESPONSE")
p.data.sub.resp$value = factor(p.data.sub.resp$value, levels = c("CR", "VGPR/PR", "SD/PD"))

leg.size = 5
p.data.geomtile.pl =
  ggplot() +
  geom_tile(
    data = p.data.sub.resp,
    aes(variable, PATIENT_ID_SHORT , fill = value),
    lwd = .5, linetype = 1, color = "white"
  ) +
  guides(fill = guide_legend(title = "Response", order = 1, override.aes = list(size = leg.size))) +
  scale_fill_manual(values = c(CR = "#7B9AB6", `VGPR/PR` = "#E9C54E" , `SD/PD`="#9B740A")) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = subset(p.data.sub, variable == "PRODUCT"),
    aes(variable, PATIENT_ID_SHORT , fill = value),
    lwd = .5, linetype = 1, color = "white"
  ) +
  guides(fill = guide_legend(title = "Product", order = 2, override.aes = list(size = leg.size))) +
  scale_fill_manual(values = c("#006EAE", "#44AA99")) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = subset(p.data.sub, variable == "CRS_GRADE"),
    aes(variable, PATIENT_ID_SHORT , fill = value),
    lwd = .5, linetype = 1, color = "white"
  ) +
  guides(fill = guide_legend(title = "CRS Grade", order = 3, override.aes = list(size = leg.size))) +
  scale_fill_manual(values = c("#6699CC", "#d98f9e", "#994455")) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = subset(p.data.sub, variable == "SEX"),
    aes(variable, PATIENT_ID_SHORT , fill = value),
    lwd = .5, linetype = 1, color = "white"
  ) +
  guides(fill = guide_legend(title = "Sex", order = 4, override.aes = list(size = leg.size))) +
  scale_fill_manual(values = c("#888888", "#DDDDDD")) +
  ggnewscale::new_scale_fill() +
  scale_x_discrete(labels = x.labels) +
  mytheme(base_size = base.size) +
  theme(
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.1),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(vjust = + 2),
    axis.ticks.y = element_blank(),
    legend.spacing.y = unit(0, 'mm'),
    legend.title = element_text(margin = margin(b = 3, unit = "pt"), face = "plain")
  ) +
  ylab("Patient ID")

sample.cells = as.data.frame(table(se.meta@meta.data$orig.ident))
pd.tp = se.meta@meta.data %>% dplyr::select(orig.ident, PATIENT_ID_SHORT, TIMEPOINT)
pd.tp = pd.tp[!duplicated(pd.tp$orig.ident), ]
pd.tp$CELLS = sample.cells$Freq[match(pd.tp$orig.ident, sample.cells$Var1)]
pd.tp$PATIENT_ID_SHORT = factor(pd.tp$PATIENT_ID_SHORT, levels = levels(p.data.sub$PATIENT_ID_SHORT))
pd.tp = pd.tp %>% dplyr::select(-orig.ident)
pd.tp = reshape2::melt(pd.tp, id.vars = c("PATIENT_ID_SHORT", "TIMEPOINT"))

nbr.cells = ggplot() +
  geom_tile(
    data = pd.tp, aes(TIMEPOINT, PATIENT_ID_SHORT , fill = value),
    lwd = .5, linetype = 1, color = "white"
  ) +
  scale_fill_gradientn(
    colors = c(scico::scico(20, palette = "navia", direction = -1), rep("#021326", 1)),
    breaks = scales::pretty_breaks(3)
  ) +
  mytheme(base_size = base.size) +
  theme(
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.1),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.spacing.y = unit(0, 'mm'),
    legend.title = element_text(margin = margin(b = 3, unit = "pt"), face = "plain")
  ) +
  guides(
    fill = guide_colourbar(
      title = "Nbr. of Cells", title.vjust = 1.1,
      barwidth = unit(.6, 'lines'), barheight = unit(7.5, 'lines'),
      ticks.linewidth = 2/.pt, ticks = T, frame.colour="black",
      frame.linewidth = 0.5/.pt
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Barplot: % of celltypes per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
p.data.ct = p.data %>%
  group_by(orig.ident, celltype_short_3) %>%
  dplyr::summarise(nbr.cells = n()) %>%
  dplyr::mutate(freq = nbr.cells / sum(nbr.cells))
p.data.ct$TIMEPOINT = se.meta$TIMEPOINT[match(p.data.ct$orig.ident, se.meta$orig.ident)]
p.data.ct$PATIENT = se.meta$PATIENT_ID_SHORT[match(p.data.ct$orig.ident, se.meta$orig.ident)]
p.data.ct$PRODUCT = se.meta$PRODUCT[match(p.data.ct$orig.ident, se.meta$orig.ident)]
p.data.ct$RESPONSE = se.meta$BEST_RESPONSE_CONSENSUS[match(p.data.ct$orig.ident, se.meta$orig.ident)]
lvls = names(ct.col)[names(ct.col) %in% p.data.ct$celltype_short_3]
p.data.ct$celltype_short_3 = factor(p.data.ct$celltype_short_3, levels = lvls)

p.data.ct$PATIENT = factor(p.data.ct$PATIENT, levels = levels(p.data.sub$PATIENT_ID_SHORT))

nbr.celltypes.pl =
  ggplot(p.data.ct, aes(freq, PATIENT, fill = celltype_short_3)) +
  geom_bar(stat="identity", width = .8) +
  facet_nested(RESPONSE + PRODUCT ~ TIMEPOINT, scales = "free", space = "free", nest_line = element_line(linetype = 1)) +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  scale_fill_manual(values = ct.col) +
  xlab("Cell type fraction") +
  mytheme(base_size = base.size) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(.25, "lines"),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(vjust = + 2),
    axis.ticks.y = element_blank()
  ) +
  ylab("Patient ID")

pd = se.meta@meta.data
mean(pd[pd$TIMEPOINT == "Late", ]$TIME_CAR_DAY_30)
min(pd[pd$TIMEPOINT == "Late", ]$TIME_CAR_DAY_30)
max(pd[pd$TIMEPOINT == "Late", ]$TIME_CAR_DAY_30)
mean(pd[pd$TIMEPOINT == "Very Late", ]$TIME_CAR_DAY_100)
min(pd[pd$TIMEPOINT == "Very Late", ]$TIME_CAR_DAY_100)
max(pd[pd$TIMEPOINT == "Very Late", ]$TIME_CAR_DAY_100)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc: T
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dim = "umap"
ct.to.pl = "celltype"
leg.size = 6

pd = get_metadata(se.t)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]
pd.o = pd

pd.cc = subset(pd, CellCycle == T)
pd.cc[[ct.to.pl]] = "Cycling"
pd = droplevels(pd[!pd$cell %in% pd.cc$cell, ])

pd.cd4 = pd[grepl("CD4", pd[[ct.to.pl]]), ]
pd.cd4[[ct.to.pl]] = factor(
  pd.cd4[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.cd4[[ct.to.pl]]])
)

pd.cd8 = pd[grepl("CD8", pd[[ct.to.pl]]), ]
pd.cd8[[ct.to.pl]] = factor(
  pd.cd8[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.cd8[[ct.to.pl]]])
)

pd.other = pd[!pd[[ct.to.pl]] %in% as.character(c(unique(pd.cd4[[ct.to.pl]]), unique(pd.cd8[[ct.to.pl]]))), ]
pd.other = rbind(pd.other, pd.cc)
pd.other[[ct.to.pl]] = factor(
  pd.other[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.other[[ct.to.pl]]])
)

stopifnot(
  length(colnames(se.t)) == ( nrow(pd.cd4) + nrow(pd.cd8) + nrow(pd.other))
)

ct.t.pl =
  ggplot() +
  scattermore::geom_scattermore(data = pd.o, aes(x = DIM1, y = DIM2), pointsize = 6, color="black")+
  scattermore::geom_scattermore(data = pd.o, aes(x = DIM1, y = DIM2), pointsize = 5, color="white") +
  # geom_point(data = pd.cd4, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  ggrastr::geom_jitter_rast(data = pd.cd4, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = 1) +
  guides(colour = guide_legend(
    title = "CD4", title.position = "top", ncol = 2, order = 1, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#BBBBBB", "Cycling")), na.value = "green") +
  ggnewscale::new_scale_color() +
  # geom_point(data = pd.cd8, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  ggrastr::geom_jitter_rast(data = pd.cd8, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = 1) +
  guides(colour = guide_legend(
    title = "CD8", title.position = "top", ncol = 2, order = 2, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#FFB92D", "Cycling")), na.value = "green") +
  ggnewscale::new_scale_color() +
  # geom_point(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  ggrastr::geom_jitter_rast(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = 1) +
  guides(colour = guide_legend(
    title = "Other", title.position = "top", ncol = 1, order = 3, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#BBBBBB", "green")), na.value = "green") +
  mytheme(base_size = base.size) +
  theme(
    legend.text = element_text(size = rel(.9)),
    legend.spacing.y = unit(1, 'mm'),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 10, hjust = 1.1, face = "plain")
  ) +
  coord_flip() + scale_y_reverse()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FACS vs scRNA-Seq: T-Cell subtypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata = readRDS("publication/clinicial_data/clinical_table_DF_2024_10_28.Rds")

t.s = pdata$pdata.t.s
t.s = t.s[t.s$DAY %in% c("Leukapheresis", "Day 30", "Day 100"), ]

# START ### Supps
t.s.supps = t.s
t.s.supps$SAMPLE_ID = NULL
t.s.supps$PATIENT_ID = gsub("Patient0", "P", t.s.supps$PATIENT_ID)
t.s.supps$DAY[t.s.supps$DAY == "Leukapheresis"] = "LP"

sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info_new/table_fc_t_subsets.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  t.s.supps,
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)
# END ### Supps

t.s = t.s[, grepl("_PERC|DAY|SAMPLE_ID", colnames(t.s))]

t.s = t.s %>% dplyr::mutate(
  DAY = dplyr::case_when(
    DAY == "Leukapheresis" ~ "LP",
    DAY == "Day 30" ~ "Late",
    DAY == "Day 100" ~ "Very Late"
  )
)

t.s$CD4_CENTRAL_MEMORY_PERC = t.s$CD4_CENTRAL_MEMORY_PERC + t.s$NAIVE_CD4_PERC
t.s$NAIVE_CD4_PERC = NULL
t.s = reshape::melt(t.s, id.vars = c("SAMPLE_ID", "DAY"))
colnames(t.s)[3:4] = c("FACS_CT", "FACS_PERC")
naturalsort(as.character(unique(t.s$FACS_CT)))
t.s = t.s %>% dplyr::mutate(
  FACS_CT_2 = dplyr::case_when(
    FACS_CT == "CD4_CENTRAL_MEMORY_PERC" | FACS_CT == "NAIVE_CD4_PERC" ~ "CD4 NaiveLike",
    FACS_CT == "TREGS_PERC" ~ "CD4 Treg",
    FACS_CT == "CD8_CENTRAL_MEMORY_PERC" ~ "CD8 CM",
    FACS_CT == "NAIVE_CD8_PERC" ~ "CD8 NaiveLike",
    FACS_CT == "CD8_EFF-MEMORY_PERC"  ~ "CD8 EM",
    FACS_CT == "CD8_EFFEKTORZELLEN_PERC" ~ "CD8 EMRA"
  )
)
t.s = t.s[!is.na(t.s$FACS_CT_2), ]
t.s$ID = paste0(t.s$SAMPLE_ID, "_", t.s$DAY, "_", t.s$FACS_CT_2)

pd = se.t@meta.data
pd = pd[grepl("^CD4|^CD8", pd$celltype_short_3), ]
pd$celltype = gsub("CD8 EMRA 1", "CD8 EMRA", pd$celltype)
pd$celltype = gsub("CD8 EMRA 2", "CD8 EMRA", pd$celltype)
pd = droplevels(pd)

pd.ct = pd %>%
  group_by(SAMPLE_ID, TIMEPOINT, celltype_short_3, celltype) %>%
  dplyr::summarise(nbr.cells = n()) %>%
  dplyr::mutate(freq = nbr.cells / sum(nbr.cells)) %>%
  dplyr::filter(celltype != "CD4 Treg") %>%
  data.frame()

pd.ct.t = pd %>%
  dplyr::mutate(celltype_short_3 = "T") %>%
  group_by(SAMPLE_ID, TIMEPOINT, celltype_short_3, celltype) %>%
  dplyr::summarise(nbr.cells = n()) %>%
  dplyr::mutate(freq = nbr.cells / sum(nbr.cells)) %>%
  dplyr::filter(celltype == "CD4 Treg") %>%
  data.frame()

pd.ct = rbind(pd.ct, pd.ct.t)

pd.ct$ID = paste0(pd.ct$SAMPLE_ID, "_", pd.ct$TIMEPOINT, "_", pd.ct$celltype)
t.s$RNA_PERC = pd.ct$freq[match(t.s$ID, pd.ct$ID)]
t.s$DAY = factor(t.s$DAY, levels = c("LP", "Late", "Very Late"))
t.s = t.s[paste0(t.s$SAMPLE_ID, t.s$DAY) %in% intersect(paste0(t.s$SAMPLE_ID, t.s$DAY), paste0(pd.ct$SAMPLE_ID, pd.ct$TIMEPOINT)), ]

fact.vs.sc.lp.t =
  ggplot(t.s, aes(RNA_PERC * 100, FACS_PERC * 100, color = DAY)) +
  geom_point(size = 1, pch = 21, stroke = 1.4, alpha = .75) +
  # scale_color_manual(values = c("#004488", "#DDAA33", "#BB5566")) +
  scale_color_manual(values = c("#88CCEE", "#DDCC77", "#CC6677")) +
  geom_abline(slope = 1, linewidth = .15, linetype = 2) +
  ggpubr::stat_cor(
    method = "spearman", label.y = sqrt(95), size = 4, color = "black"
  ) +
  mytheme(base_size = base.size) +
  theme(
    legend.position = "right",
    legend.box.margin=margin(0,0,0,-15),
    legend.text = element_text(margin = margin(l = 1), size = 13),
    axis.text = element_text(size = 13),
    panel.spacing = unit(.5, "lines")
  ) +
  scale_y_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  scale_x_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  xlab("% of cell type (by scRNA-Seq)") +
  ylab("% of cell type (by FC)") +
  facet_grid2(~ FACS_CT_2, independent = F) +
  guides(
    color = guide_legend(
      title = NULL, override.aes = list(shape = 16, size = 3)
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FACS vs scRNA-Seq: CAR+
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata.facs = pdata$pdata.imm.s
pdata.facs = pdata.facs[pdata.facs$DAY == "Day 30" | pdata.facs$DAY == "Day 100", ]
pdata.facs$DAY = ifelse(pdata.facs$DAY == "Day 30", "Late", "Very Late")

pdata.facs$ID = paste0(pdata.facs$SAMPLE_ID, "_", pdata.facs$DAY)

pd = se.t@meta.data
pd = pd[grepl("^CD4|^CD8", pd$celltype_short_3), ]
pd = droplevels(pd)
pd$CAR_BY_EXPRS = as.logical(pd$CAR_BY_EXPRS)

pd.ct =
  pd %>%
  dplyr::filter(TIMEPOINT %in% c("Late", "Very Late")) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(nbr.cells = n()) %>%
  dplyr::mutate(car_pos = sum(CAR_BY_EXPRS)) %>%
  dplyr::mutate(CARpos_PERC_RNA = (car_pos / nbr.cells)*100) %>%
  dplyr::distinct(orig.ident, .keep_all = T) %>%
  data.frame()
pd.ct$ID = paste0(pd.ct$SAMPLE_ID, "_", pd.ct$TIMEPOINT)
pd.ct$CARpos_PERC_FACS = pdata.facs$CD3_CAR_PERC[match(pd.ct$ID, pdata.facs$ID)]
pd.ct = pd.ct[!is.na(pd.ct$CARpos_PERC_FACS), ]

pd.ct$PRODUCT = ifelse(pd.ct$PRODUCT == "ide", "Ide-cel", "Cilta-cel")

fact.vs.sc.car =
  ggplot(pd.ct, aes(CARpos_PERC_RNA, CARpos_PERC_FACS, fill = PRODUCT)) +
  geom_point(size = 2.5, pch = 21, stroke = .4, colour = "black", alpha = .45) +
  ggpubr::stat_cor(method = "spearman", label.y = sqrt(95), size = 3.5, color = "black") +
  scale_y_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  scale_x_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  scale_fill_manual(values = c("#004488", "#009988")) +
  geom_abline(slope = 1, linewidth = .3, linetype = "dashed") +
  xlab("% of CAR+ (by scRNA-Seq)") +
  ylab("% of CAR+ (by FC)") +
  mytheme(base_size = base.size) +
  theme(
    # aspect.ratio = 1,
    legend.position = "none",
    axis.text = element_text(size = 13),
    legend.title = element_blank(),
    panel.spacing = unit(.5, "lines")
  ) +
  facet_grid2(TIMEPOINT ~ PRODUCT, independent = F)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Survival analysis with cell type fractions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd.tp = se.meta@meta.data

pd.tp$celltype = pd.tp$celltype
pd.tp = droplevels(pd.tp)

keep.ct = pd.tp %>% dplyr::group_by(celltype, TIMEPOINT) %>%
  dplyr::count() %>%
  dplyr::filter(n > 200) %>%
  dplyr::group_by(celltype) %>%
  dplyr::count() %>%
  dplyr::filter(n >= 1) %>%
  dplyr::pull(celltype) %>%
  as.character()

pd.tp = pd.tp[pd.tp$celltype %in% keep.ct, ]
pd.tp = split(pd.tp, pd.tp$TIMEPOINT)

# file = pd.tp$Late
cox.ct.res = parallel::mclapply(pd.tp, function(file) {

  TP = as.character(file$TIMEPOINT[1])
  pd = droplevels(file)

  prop = table(pd$celltype, pd$orig.ident)
  prop = prop[rowSums(prop) > 200, ]
  cell.nbr = rowSums(prop)
  prop = prop.table(prop, margin = 2)
  prop = as.data.frame.matrix(prop)
  prop.ori = prop
  prop = as.matrix(prop)
  prop = standardize(prop)

  pheno = pd
  pheno = pheno[!duplicated(pheno$orig.ident), ]
  rownames(pheno) = pheno$orig.ident
  pheno = pheno[base::intersect(rownames(pheno), colnames(prop)), ]
  prop = prop[, base::intersect(rownames(pheno), colnames(prop))]
  stopifnot(identical(rownames(pheno), colnames(prop)))

  # Survival object
  survivalDOD = Surv(pheno$PFS, pheno$PROGRESSION)
  survModel = survfit(survivalDOD ~ 1, data = pheno)

  run_cox = function(obj) {

    exprs.mat.spl = split(as.matrix(obj), rownames(obj))

    coxModel.l = parallel::mclapply(exprs.mat.spl, function(celltype) {
      res = summary(survival::coxph(survivalDOD ~ celltype))
      log.rank =  tryCatch(
        survival::survdiff(
          survivalDOD ~ ifelse(celltype > median(celltype), 2, 1)
        )$pvalue, error=function(e) "error"
      )
      if(log.rank == "error") {
        log.rank = 1
      }

      res = c(
        "HR" = res$coef[1, , drop = F][2],
        "logHR" = res$coef[1, , drop = F][1],
        "SE_logHR" = res$coef[1, , drop = F][3],
        "L95CI" = res$conf.int[1, , drop = F][,"lower .95"],
        "U95CI" = res$conf.int[1, , drop = F][,"upper .95"],
        "Rsquare" = res$rsq[[1]],
        "LR_Pval" = unname(res$logtest[3]),
        "LogRank_Pval" = log.rank,
        "Pval" = res$coef[1, , drop = F][5],
        "PH" =  unname(
          cox.zph(survival::coxph(survivalDOD ~ celltype))$table[, 3][2]
        )
      )
    }, mc.cores = 1)

    coxModel = data.frame(do.call("rbind", coxModel.l))
    coxModel
  }

  # Cox survival model fitting
  coxModel = run_cox(obj = prop)
  coxModel$TIMEPOINT = TP
  coxModel = coxModel[order(coxModel$Pval, decreasing = F), ]
  coxModel$CELLNBR = cell.nbr[match(rownames(coxModel), names(cell.nbr))]
  coxModel$LogRank_Pval_adj = p.adjust(coxModel$LogRank_Pval, method = "BH")
  coxModel$Pval_adj = p.adjust(coxModel$Pval, method = "BH")

  prop = prop[rownames(coxModel), ]

  list(coxModel = coxModel, prop = prop, pheno = pheno, prop.ori = prop.ori)

}, mc.cores = 1)

res = lapply(names(cox.ct.res), function(x){
  cox = cox.ct.res[[x]]$coxModel
  cox$CELLTYPES = rownames(cox)
  cox
})
df.surv = do.call("rbind", res)
df.surv$TIMEPOINT = gsub("Apheresis", "LP", df.surv$TIMEPOINT)
df.surv$TIMEPOINT = factor(
  df.surv$TIMEPOINT, levels = c("LP", "Late", "Very Late")
)
df.surv$logHR_DIR = ifelse(df.surv$logHR > 0, ">0", "<0")

logrank.hr.pl =
ggplot(df.surv, aes(x = TIMEPOINT, y = CELLTYPES, label = CELLNBR)) +
  scale_x_discrete() +
  scale_y_discrete() +
  geom_point(
    aes(shape = ifelse(LogRank_Pval >= 0.1, NA, "alpha.s")),
    size = 3, stroke = 2.5, color = "#228833"
  ) +
  geom_point(
    aes(shape = ifelse(LogRank_Pval_adj >= 0.1, NA, "fdr.s")),
    size = 3, stroke = 2.5, color = "red"
  ) +
  # scale_shape_manual(
  #   values = c(22), name = paste0("Log Rank p < 0.1"), labels = c("yes", "")
  # ) +
  scale_shape_manual(
    values = c(22, 22),
    name = NULL,
    labels = c(paste0("p < ", 0.1), paste0("FDR < ", 0.1), "")
  ) +
  geom_point(aes(color = logHR_DIR), size = 4, shape = 15) +
  scale_color_manual(values = c(">0" = "#997700", "<0" = "#6699CC")) +
  geom_text(nudge_x = .15, hjust = 0, size = 3) +
  mytheme() +
  theme(
    panel.grid.major = element_line(
      colour = "grey80", linetype = 2, linewidth = .3
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = rel(1)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(
      margin = margin(l = 1,  unit = "pt"), size = rel(1)
    ),
    legend.margin = margin(t=-3, l = 5),
  ) +
  guides(
    shape = guide_legend(
      order = 3, title.position = "top",
      override.aes = list(size = 3, stroke = 1)
    ),
    color = guide_legend(
      title = "Log hazard ratio", order = 1,
      title.hjust = 0, title.position = "top",
      override.aes = list(shape = 16, size = 4)
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | cell type correlation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# x = "Very Late"
cor.l = lapply(names(cox.ct.res), function(x){
  df = cox.ct.res[[x]]$prop.ori
  df = df[naturalorder(rownames(df)), ]
  corr = cor(t(df),method="spearman",use="pairwise.complete.obs")
  p.mat <- ggcorrplot::cor_pmat(t(df),method="spearman",exact=FALSE)
  pt = ggcorrplot::ggcorrplot(
    corr = corr,
    lab=T,
    show.legend = F,
    type = "lower",
    colors = c("#4477AA","#FFFFFF","#CC6677"),
    outline.color = "#AAAAAA",
    lab_size = 1,
    insig = "pch",
    pch.cex = 4,
    pch.col = "#AAAAAA",
    tl.cex = 7,
    p.mat = p.mat,
    hc.order = F
  ) +
    ggtitle(x)
  pt
})

ggsave2(
  filename="publication/extended_data_files/Fig_02_celltypes_cor.png",
  plot_grid(plotlist = cor.l, nrow = 1, scale = .95),
  width = 180, height = 60, dpi = 300,
  bg = "white", units = "mm", scale = 2.5,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SUPP | DA | grouped by outcome
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
celltype = "celltype"
p.data = se.meta@meta.data

l = list()
for (i in as.character(unique(p.data$BEST_RESPONSE_CONSENSUS))) {

  df = droplevels(p.data[p.data$BEST_RESPONSE_CONSENSUS == i, ])
  # df = droplevels(df[grepl("B", df$celltype), ])

  prop.ct = speckle::getTransformedProps(
    clusters = df[[celltype]],
    sample = df$orig.ident
  )$Proportions

  prop.ct = reshape2::melt(prop.ct)
  prop.ct$TIMEPOINT = df$TIMEPOINT[match(prop.ct$sample, df$orig.ident)]
  prop.ct$OUTCOME = i
  l[[i]] = prop.ct

}

prop.ct = rbind(l$`non-CR`, l$`CR`)

p.data.w = p.data
p.data.w$TIMEPOINT = gsub(" ", "", p.data.w$TIMEPOINT)
p.data.rspns = split(p.data.w, p.data.w$BEST_RESPONSE_CONSENSUS)
# For some reason, variables that are passed to the limma::makeContrasts
# function must be defined globally
grp1 = "LP"; grp2 = "Late"; da.late.lp = da_paired(p.data.rspns, transform = "logit")
grp1 = "Late"; grp2 = "VeryLate"; da.vlate.late = da_paired(obj = p.data.rspns, transform = "logit")
grp1 = "LP"; grp2 = "VeryLate"; da.vlate.lp = da_paired(p.data.rspns, transform = "logit")

da.tps = rbind(da.late.lp, da.vlate.late, da.vlate.lp)
da.tps$CTRST = gsub("VeryLate", "Very Late", da.tps$CTRST)
da.tps$group1 = gsub("_vs.+", "", da.tps$CTRST)
da.tps$group2 = gsub(".+_vs_", "", da.tps$CTRST)
da.tps = da.tps[da.tps$Present.Grp1 > 2 | da.tps$Present.Grp2 > 2, ]
da.tps$ID = paste0(
  da.tps$CLUSTER, da.tps$group1, da.tps$group2, da.tps$GROUP
)

stat.test <- prop.ct %>%
  dplyr::group_by(OUTCOME, clusters) %>%
  rstatix::t_test(value ~ TIMEPOINT)

stat.test <- stat.test %>%
  rstatix::add_xy_position(x = "OUTCOME", dodge = 0.8, step.increase = .25)

stat.test$y.position = stat.test$y.position + 3.4

stat.test$ID = paste0(
  stat.test$clusters, stat.test$group1, stat.test$group2, stat.test$OUTCOME
)
stat.test$p = da.tps$P.Value[match(stat.test$ID, da.tps$ID)]
stat.test$p.adj = da.tps$FDR[match(stat.test$ID, da.tps$ID)]

stat.test = stat.test %>%  rstatix::add_significance(
  "p.adj",
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***",  "**", "*", "ns")
)
stat.test$p.adj.signif[stat.test$p.adj.signif == ""] = "ns"
stat.test$p.adj = round(stat.test$p.adj, 5)
stat.test$p.adj[stat.test$p.adj.signif == "ns"] = "ns"

k = table(stat.test$clusters, stat.test$p.adj < .1)[, 2]
k = names(k[k > 0])

prop.ct = prop.ct[prop.ct$clusters %in% k, ]

pl =
  ggplot(prop.ct, aes(OUTCOME, value, color = TIMEPOINT)) +
  geom_boxplot(
    outliers = F, fatten = 1.5, alpha = 1, lwd = 0.5
  ) +
  facet_wrap(~ clusters, nrow = 3, scales = "free_x") +
  ggpubr::stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", hide.ns = TRUE# , size = 2.1
  ) +
  geom_point(
    aes(group=TIMEPOINT),
    position=position_dodge(width=0.75), size = .1
  ) +
  scale_color_manual(values = c("#555555", "#6699CC", "#004488")) +
  xlab(NULL) + labs(fill = NULL) +
  ylab("Proportion") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.3)),
    axis.title.x = element_blank(),
    # axis.text.x = element_text(angle=45, vjust=1, hjust=1),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(.25, "lines"),
    strip.text = element_text(size = rel(1), face = "plain"),
    strip.background = element_blank()
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10),
    breaks=c(0, 0.001, 0.01, 0.1, 1) , limits = c(0, 22)
  ) +
  labs(color = NULL)

ggsave2(
  filename="publication/extended_data_files/Fig_02_celltypes_over_time.png",
  pl,
  width = 180, height = 80, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Survival curves
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
x.breaks = 150
do.trunc = T
trunc.cut = 2
tp = "Late"

pheno = cox.ct.res[[tp]]$pheno
coxModel = cox.ct.res[[tp]]$coxModel
prop = cox.ct.res[[tp]]$prop

prop = prop[rownames(coxModel[coxModel$LogRank_Pval < 0.1, ]), ]

a = coxModel[coxModel$logHR > 0, ]
a = rownames(a[order(a$LogRank_Pval, decreasing = F), ])
a = a[!grepl("gdT", a)]

b = rbind(
  coxModel[coxModel$logHR < 0, ],
  coxModel[grepl("gdT", rownames(coxModel)), ]
)
b = rownames(b[order(b$LogRank_Pval, decreasing = F), ])
prop = prop[c(a, b)[c(a, b) %in% rownames(prop)], ]

# ct = "Plasma cell"
kapl.l = list()
for (ct in rownames(prop)) {

  ftr = prop[ct, ]
  pheno$GROUP = NULL
  pheno$GROUP = ifelse(ftr > median(ftr), 1, 0)
  survivalDOD = Surv(pheno$PFS, pheno$PROGRESSION)
  survModel = survfit(survivalDOD ~ GROUP, data = pheno)

  # ungrouped KM estimate to determine cutoff
  survModel.ungrouped = survminer::surv_fit(survivalDOD ~ 1, data = pheno)

  if (do.trunc == T) {
    cutoff  = min(
      survModel.ungrouped$time[survModel.ungrouped$n.risk <= trunc.cut]
    )
  }

  cust.th = theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
    plot.margin = unit(c(1,17,7,7), "points")
  )

  log.rank.pval = coxModel[ct, ]$LogRank_Pval
  if (log.rank.pval <= 0.001) {
    log.rank.pval = formatC(log.rank.pval, format = "e", digits= 2)
  } else if (log.rank.pval > 0.001 & log.rank.pval <= 0.01) {
    log.rank.pval = format(round(log.rank.pval, digits=5), nsmall = 5)
  } else {
    log.rank.pval = format(round(log.rank.pval, digits=4), nsmall = 3)
  }

  log.rank.pval.adj = coxModel[ct, ]$LogRank_Pval_adj
  if (log.rank.pval.adj <= 0.0001) {
    log.rank.pval.adj = formatC(log.rank.pval.adj, format = "e", digits= 2)
  } else if (log.rank.pval.adj > 0.0001 & log.rank.pval.adj <= 0.001) {
    log.rank.pval.adj = format(round(log.rank.pval.adj, digits=5), nsmall = 5)
  } else {
    log.rank.pval.adj = format(round(log.rank.pval.adj, digits=4), nsmall = 3)
  }

  surv.pl =
    survminer::ggsurvplot(
      survModel,
      conf.int = F,
      pval = F,
      break.time.by = x.breaks,
      legend = "none",
      xlab = NULL,
      ylab = NULL,
      censor.size = 3,
      censor.shape = 124,
      size = 1,
      ggtheme = mytheme(base_size = base.size) + cust.th,
      tables.theme =  survminer::theme_cleantable(),
      risk.table.fontsize = 3.5,
      risk.table.title = "",
      risk.table.pos = "in",
      risk.table = "absolute",
      palette = c("#6699CC", "#BB5566"),
      title = paste0(
        ct, " | Late\np: ", log.rank.pval, "; FDR: ", log.rank.pval.adj
      )
    )

  surv.pl$plot = surv.pl$plot +
    theme(
      axis.text = element_text(size = rel(.8))
    )

  suppressMessages(
    if (do.trunc == T) {
      surv.pl$plot = surv.pl$plot +
        coord_cartesian(xlim=c(0, cutoff)) +
        scale_x_continuous(
          limits = c(0, cutoff), breaks = seq(0, cutoff, by = x.breaks)
        )
    })

  surv.pl$table = surv.pl$table +
    theme(
      panel.border = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank(),
      plot.margin =  margin(-10, 0, 0, 0, unit = "pt")
    )

  suppressMessages(
    if (do.trunc == T) {
      surv.pl$table = surv.pl$table +
        coord_cartesian(xlim=c(0, cutoff)) +
        scale_x_continuous(
          limits = c(0, cutoff), breaks = seq(0, cutoff, by = x.breaks)
        )
    })

  surv.pl =
    surv.pl$plot +
    patchwork::inset_element(surv.pl$table, left = 0, bottom = .01, right = 1, top =.4)
  kapl.l[[ct]] = surv.pl
}

kapl = plot_grid(plotlist = kapl.l, ncol = 4, scale = 1)
kapl =
  annotate_figure(
    kapl,
    left = text_grob("Progression-free survival", rot = 90, size = base.size),
    bottom  = text_grob("Days", rot = 0, size = base.size, vjust = -.5)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# forest, logHR + CI
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = lapply(names(cox.ct.res), function(x){
  cox = cox.ct.res[[x]]$coxModel
  # cox$L95CI = cox$logHR - cox$SE_logHR * 1.96
  # cox$U95CI = cox$logHR + cox$SE_logHR * 1.96
  cox$L95CI = exp(cox$logHR - cox$SE_logHR * 1.65)
  cox$U95CI = exp(cox$logHR + cox$SE_logHR * 1.65)
  cox$CELLTYPE = rownames(cox)
  cox
})
df = do.call("rbind", res)

keep.ct = reshape2::melt(table(df$CELLTYPE, df$Pval < 0.1)) %>%
  dplyr::filter(Var2 == T) %>%
  dplyr::filter(value > 0) %>%
  dplyr::pull(Var1) %>%
  as.character()
df = df[df$CELLTYPE %in% keep.ct, ]

order.ct = df %>%
  dplyr::group_by(CELLTYPE) %>%
  dplyr::summarise(ave = mean(logHR)) %>%
  dplyr::arrange(ave) %>%
  dplyr::pull(CELLTYPE)
df$CELLTYPE = factor(df$CELLTYPE, levels = rev(order.ct))

df$Cox_Pval_Sign = ifelse(df$Pval < 0.1, "#549C74", "#555555")

dotCOLS = c("#9c9c9c", "#6699CC", "#004488")
barCOLS = c("#9c9c9c", "#6699CC", "#004488")

df$Cox_Pval
l = list()
for (i in 1:length(df$Pval)) {
  if (df$Pval[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval[i] > 0.0001 & df$Pval[i] <= 0.001) {
    l[[i]] = format(round(df$Pval[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_pl = unlist(l)
df$Cox_Pval_pl = ifelse(
  df$Pval < 0.1, paste0(df$Cox_Pval_pl, " *"), df$Cox_Pval_pl
)

l = list()
for (i in 1:length(df$Pval_adj)) {
  if (df$Pval_adj[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval_adj[i] > 0.0001 & df$Pval_adj[i] <= 0.001) {
    l[[i]] = format(round(df$Pval_adj[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval_adj[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_adj_pl = unlist(l)
df$Cox_Pval_adj_pl = ifelse(
  df$Pval_adj < 0.1, paste0(df$Cox_Pval_adj_pl, " *"), df$Cox_Pval_adj_pl
)

df$TIMEPOINT = factor(df$TIMEPOINT, levels = c("LP", "Late", "Very Late"))

hr.pl =
ggplot(
    df,
    aes(
      x = CELLTYPE, y = HR, ymin = L95CI, ymax = U95CI,
      col = TIMEPOINT, fill=TIMEPOINT)
  ) +
  geom_hline(yintercept=1, lty=2, lwd = .5, col = "#555555") +
  geom_linerange(linewidth=1.5, position=position_dodge(width = 0.8)) +
  facet_wrap(~ CELLTYPE, ncol = 1, scales = "free_y") +
  geom_point(
    size=3, shape=21, colour="white", stroke = .5,
    position=position_dodge(width = 0.8),
  ) +
  geom_text(
    data = df[!df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = U95CI + 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = -0.1, vjust = 0.5, size = rel(3.25)
    # hjust = -0.1, vjust = 0.5, size = rel(2.5)
  ) +
  geom_text(
    data = df[df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = L95CI - 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = 1.1, vjust = 0.5, size = rel(3.25)
    # hjust = 1.1, vjust = 0.5, size = rel(2.5)
  ) +
  coord_flip() +
  scale_colour_manual(values = alpha(barCOLS, .7), drop=F) +
  scale_fill_manual(values = dotCOLS) +
  scale_y_log10(
    name= "Hazard ratio (90% CI)",
    limits = c(min(df$L95CI), max(df$U95CI) + 0)
  ) +
  xlab(NULL) +
  mytheme(base_size = base.size) +
  # mytheme() +
  theme(
    strip.text.x  = element_blank(),
    panel.spacing  = unit(.4, "lines"),
    legend.title=element_blank(),
    legend.box.margin=margin(-5,0,5,0),
    legend.position = "right",
    axis.ticks.y = element_blank(),
    legend.text = element_text(margin = margin(l = 3, unit = "pt"), size = rel(1))
  ) +
  scale_x_discrete(expand = c(0.5, 0))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supps
# logrank and uni cox model
# multivariate cox model
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("publication/figure_scripts/fig_02_suppl_cox_multivar_source.R")
ggsave2(
  filename="publication/extended_data_files/Fig_2_multi_cox.jpg",
  plot_grid(
    logrank.hr.pl, NULL, hr.multi.2.pl,
    nrow = 1, rel_widths = c(1, .15, 1),
    align = "vh", labels = c("e", "", "f"),
    label_fontface = "bold", label_size = 11
  ),
  width = 180, height = 130, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# ggsave2(
#   filename="publication/extended_data_files/Fig_2_multi_cox_rebuttal.jpg",
#   plot_grid(
#     hr.pl, NULL, hr.multi.1.pl, NULL, hr.multi.2.pl,
#     nrow = 1, rel_widths = c(1, .1, 1, .1, 1),
#     align = "vh", labels = c("a", "", "b", "", "c"),
#     label_fontface = "bold", label_size = 11
#   ),
#   width = 180, height = 130, dpi = 300,
#   bg = "white", units = "mm", scale = 1.6,
#   device = png, type = "cairo"
# )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fin
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mtcars$cyl = as.factor(mtcars$cyl)
mtcars = mtcars[mtcars$cyl != 4, ]
score.leg =
  get_legend(
    ggplot(mtcars, aes(x=wt, y=mpg, color=cyl)) +
      geom_point(shape="|") +
      geom_smooth(method=lm, aes(fill=cyl), se=FALSE, linewidth = 2) +
      mytheme(base_size = 15) +
      theme(
        legend.position = "right",
        legend.text = element_text(margin = margin(l = -1, r = 5, unit = "pt"), size = rel(1)),
        legend.margin = margin(t=-25)
      ) +
      guides(fill = guide_legend(
        title.position="top",
        title.vjust = 0.5,
        #title.hjust = 1,
        keywidth = 1.5,
        keyheight = .8,
        override.aes = list(size = .7, alpha = .4))
      ) +
      scale_fill_manual(
        name = "Cell type\nproportion",
        labels = c("Low", "High"),
        # labels = c(expression("" <= 0), expression("" > 0)),
        values = c("#6699CC", "#BB5566")) +
      scale_colour_manual(
        name = "Cell type\nproportion",
        labels = c("Low", "High"),
        # labels = c(expression("" <= 0), expression("" > 0)),
        values = c("#6699CC", "#BB5566"))
  )
# plot_grid(score.leg)

ggsave2(
  filename="publication/figures_main/fig_02.png",
  plot_grid(
    plot_grid(
      plot_grid(NULL, plot_grid(
        NULL,
        plot_grid(p.data.geomtile.pl, nbr.cells, rel_widths = c(1, .8)),
        nrow = 2, rel_heights = c(.05, 1)), ncol = 3, rel_widths = c(0, .5, .0)
      ),
      plot_grid(
        NULL,
        ct.pl + theme(legend.position = "none"),
        plot_grid(NULL, plot_grid(cowplot::get_legend(ct.pl)), rel_widths = c(.7, 1)),
        nrow = 3, rel_heights = c(.05, 0.5, 1)
      ),
      NULL,
      NULL,
      plot_grid(
        NULL,
        ftrs.pl,
        NULL,
        plot_grid(
          cc.pl, NULL, vdj.pl, NULL, mt.pl, nrow = 1, rel_widths = c(1, .05, 1, .05, 1),
          labels = c("D", "", "E", "", "F"), label_fontface = "bold", label_size = 24, hjust = 0.1
        ),
        nrow = 4, rel_heights = c(.075, 1.5, .1, 1)
      ),
      ncol = 5, rel_widths = c(1.5, .5, .1, .05, 1.75),
      labels = c("A", "B", "", "C", ""), label_fontface = "bold", label_size = 24,  hjust = -.4
    ),
    # NULL,
    NULL,
    plot_grid(
      nbr.celltypes.pl,
      NULL,
      plot_grid(
        plot_grid(
          ct.t.pl, NULL, fact.vs.sc.car, ncol = 3, rel_widths = c(1.8, .1, 1),
          labels = c("H", "", "J"), label_fontface = "bold", label_size = 24, vjust = .5
        ),
        NULL,
        plot_grid(fact.vs.sc.lp.t),
        nrow = 3, rel_heights = c(1.35, .1, 1),
        labels = c("", "", "I"), label_fontface = "bold", label_size = 24, vjust = .5
      ),
      ncol = 3, rel_widths = c(.3, .025, .618),
      labels = c("G", "", "", "", ""), label_fontface = "bold", label_size = 24, vjust = .5
    ),
    NULL,
    plot_grid(
      plot_grid(kapl, score.leg, ncol = 2, rel_widths = c(1, .1)),
      NULL,
      plot_grid(NULL, hr.pl, NULL, nrow = 3, rel_heights = c(.0, 1, .005)),
      ncol = 3, rel_widths = c(.618, .015,.382),
      labels = c(" K", "", "I"), label_fontface = "bold", label_size = 24, vjust = .4
    ),
    # NULL,
    nrow = 5, rel_heights = c(1.05, .075, 1, .075, 1.35)
  ),
  width = 165, height = 205, dpi = 300, bg = "white", units = "mm", scale = 3,
  device = png, type = "cairo"
)


ggsave2(
  filename="publication/figures_main/fig_02.pdf",
  plot_grid(
    plot_grid(
      plot_grid(NULL, plot_grid(
        NULL,
        plot_grid(p.data.geomtile.pl, nbr.cells, rel_widths = c(1, .8)),
        nrow = 2, rel_heights = c(.05, 1)), ncol = 3, rel_widths = c(0, .5, .0)
      ),
      plot_grid(
        NULL,
        ct.pl + theme(legend.position = "none"),
        plot_grid(NULL, plot_grid(cowplot::get_legend(ct.pl)), rel_widths = c(.7, 1)),
        nrow = 3, rel_heights = c(.05, 0.5, 1)
      ),
      NULL,
      NULL,
      plot_grid(
        NULL,
        ftrs.pl,
        NULL,
        plot_grid(
          cc.pl, NULL, vdj.pl, NULL, mt.pl, nrow = 1, rel_widths = c(1, .05, 1, .05, 1),
          labels = c("D", "", "E", "", "F"), label_fontface = "bold", label_size = 24, hjust = 0.1
        ),
        nrow = 4, rel_heights = c(.075, 1.5, .1, 1)
      ),
      ncol = 5, rel_widths = c(1.5, .5, .1, .05, 1.75),
      labels = c("A", "B", "", "C", ""), label_fontface = "bold", label_size = 24,  hjust = -.4
    ),
    # NULL,
    NULL,
    plot_grid(
      nbr.celltypes.pl,
      NULL,
      plot_grid(
        plot_grid(
          ct.t.pl, NULL, fact.vs.sc.car, ncol = 3, rel_widths = c(1.8, .1, 1),
          labels = c("H", "", "J"), label_fontface = "bold", label_size = 24, vjust = .5
        ),
        NULL,
        plot_grid(fact.vs.sc.lp.t),
        nrow = 3, rel_heights = c(1.35, .1, 1),
        labels = c("", "", "I"), label_fontface = "bold", label_size = 24, vjust = .5
      ),
      ncol = 3, rel_widths = c(.3, .025, .618),
      labels = c("G", "", "", "", ""), label_fontface = "bold", label_size = 24, vjust = .5
    ),
    NULL,
    plot_grid(
      plot_grid(kapl, score.leg, ncol = 2, rel_widths = c(1, .1)),
      NULL,
      plot_grid(NULL, hr.pl, NULL, nrow = 3, rel_heights = c(.0, 1, .005)),
      ncol = 3, rel_widths = c(.618, .015,.382),
      labels = c(" K", "", "L"), label_fontface = "bold", label_size = 24, vjust = .4
    ),
    # NULL,
    nrow = 5, rel_heights = c(1.05, .075, 1, .075, 1.35)
  ),
  width = 165, height = 205, dpi = 300, bg = "white", units = "mm", scale = 3
  # device = cairo_pdf
)

