print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Figure 3")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggh4x", "ggrepel", "scico", "DescTools",
  "scCustomize", "ggpubr", "immunarch", "Ckmeans.1d.dp", "PieGlyph", "ggstats", "survival",
  "openxlsx"
)
.bioc_packages = c(
  "dittoSeq", "SummarizedExperiment", "slingshot", "destiny"
)

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

if (any(!"Startrac" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("Japrin/sscVis")
  devtools::install_github("ncborcherding/Startrac")
}
library(Startrac)

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)

# library(escape)

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/trajectory.R")
source("code/helper/ora.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.t = readRDS(
  paste0(manifest$meta_pub$work, "integration/06_seurat_harmony_t_all_new.Rds")
)
vdj.t = readRDS(paste0(manifest$meta_pub$work, "integration/05_vdj_t_new.Rds"))
se.t = AddMetaData(se.t, vdj.t)

mir = readRDS("data/signatures/mir_gene_set_collection.Rds")

# Clone Pseudo IDs
pd = se.t@meta.data
pd$barcode = rownames(pd)
pd = pd[!is.na(pd$CTstrict), ]
df = pd %>%
  dplyr::group_by(CTstrict) %>%
  dplyr::count(name = "cloneFreqAll") %>%
  dplyr::arrange(desc(cloneFreqAll))
df$pseudo_id = paste0("Clone_", seq(1:nrow(df)))
se.t$CLONE_PSEUDO_ID = df$pseudo_id[match(se.t$CTstrict, df$CTstrict)]

dm.cd4 = readRDS(paste0(manifest$meta_pub, "diffusion/dm_cd4_new.Rds"))
dm.cd8 = readRDS(paste0(manifest$meta_pub, "diffusion/dm_cd8_new.Rds"))

dm.cd8$se$celltype = gsub("CD8 EMRA KLRC2\\+", "CD8 EMRA 2", dm.cd8$se$celltype)
dm.cd8$se$celltype = gsub("CD8 EMRA$", "CD8 EMRA 1", dm.cd8$se$celltype)
dm.cd8$se$celltype = factor(dm.cd8$se$celltype)

dm.cd8$sce$celltype = gsub("CD8 EMRA KLRC2\\+", "CD8 EMRA 2", dm.cd8$sce$celltype)
dm.cd8$sce$celltype = gsub("CD8 EMRA$", "CD8 EMRA 1", dm.cd8$sce$celltype)
dm.cd8$sce$celltype = factor(dm.cd8$sce$celltype)

dm.cd4$se = cluster_nearest_milestones(dm.cd4$sce, dm.cd4$se)
dm.cd8$se = cluster_nearest_milestones(dm.cd8$sce, dm.cd8$se)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | DGEA | non-Cr vs. CR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dg.outcome = dgea_gex(
  obj = se.t, group = "BEST_RESPONSE_CONSENSUS",
  ctrs.grp1 = "non-CR", ctrs.grp2 = "CR",
  latent.vars = c("STUDY", "nFeature_RNA"),
  split.by.tp = T, logfc.threshold = log2(1.1),
  subsample = T, threads = 20, min.de.genes = 10
)

dg = dg.outcome
dg$res.dgea.sign = dg$res.dgea.sign[
  abs(dg$res.dgea.sign$avg_log2FC) > log2(1.25),
]
table(dg$res.dgea.sign$cluster)


dg$res.dgea[grepl("BCMA", dg$res.dgea$feature), ]

# START Supp Table
supps.de.sign = dg$res.dgea.sign %>% dplyr::select(
  Gene_name = feature, log2FC = avg_log2FC, Pvalue = p_val, FDR = p_val_adj,
  Time_point = timepoint, Cell_identity = celltype
)
supps.de.sign = supps.de.sign %>% dplyr::arrange(Time_point, Cell_identity, FDR)

sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info/table_de_all_t_noncr_vs_cr.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  supps.de.sign,
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)
# END Supp Table

ftrs = c("CAR-BCMA", bm.ftrs)
ftrs = unique(ftrs[ftrs != ""])

k = sort(table(dg$res.dgea.sign$cluster)) >= 10
k = names(k[k])
dg$res.dgea.sign = dg$res.dgea.sign[
  dg$res.dgea.sign$cluster %in% k,
]
table(dg$res.dgea.sign$cluster)


max(dg$res.dgea[dg$res.dgea$timepoint == "LP", ]$avg_log2FC)
min(dg$res.dgea[dg$res.dgea$timepoint == "LP", ]$avg_log2FC)
max(dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ]$avg_log2FC)
min(dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ]$avg_log2FC)

dgea.pl.lp =
  dgea_plot(
    dgea.res = dg$res.dgea[dg$res.dgea$timepoint == "LP", ],
    dgea.res.sign = dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "LP", ], nbr.tops = 5,
    box.padding = .1, text.repel.size = 2.2, text.de.nbr.size = 2.2, cluster = "celltype", dge.nbr.up = 0,
    # ylim.extend.up = 1, ylim.extend.dn = 4,
    subset.tops = ftrs, axis.max = 3, axis.min = -3, legend.margin.t = -25
  ) +
  ggtitle("LP | DE genes comparing non-CR with CR\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))

dgea.pl.l =
dgea_plot(
    dgea.res = dg$res.dgea[dg$res.dgea$timepoint == "Late", ],
    dgea.res.sign = dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ], nbr.tops  = 5,
    box.padding = .125, text.repel.size = 2, text.de.nbr.size = 2, cluster = "celltype", dge.nbr.up = 0,
    # ylim.extend.up = 4, ylim.extend.dn = 2,
    subset.tops = ftrs,
    axis.max = 2.7, axis.min = -5, legend.margin.t = -15
  )  +
  ggtitle("Late | DE genes comparing non-CR with CR\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))

dgea.pl.vl =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Very Late", ], nbr.tops = 5,
    box.padding = .1, text.repel.size = 2.2, cluster = "celltype", dge.nbr.up = 0,
    subset.tops = ftrs, axis.max = 6.5, axis.min = -4.5, legend.margin.t = -25
  ) +
  ggtitle("Very Lates | DE genes comparing non-CR with CR\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | ORA | non-Cr vs. CR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge = dg$res.dgea.sign
dge$cluster = gsub("CD4\\.", "CD4\n", dge$cluster)
dge$cluster = gsub("CD8\\.", "CD8\n", dge$cluster)

ftrs.l = dge
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
ftrs.l.1 = ftrs.l[grepl("CD4|CD8", names(ftrs.l))]
lengths(ftrs.l.1)

ora.go = parallel::mclapply(ftrs.l.1, function(x){
  run_nmf_ora(
    genes = x, universe = rownames(se.t),
    category = "C5", subcategory = "BP"
  )
}, mc.cores = length(ftrs.l.1))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | DGEA | ORA | non-Cr vs. CR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora.go.immu = lapply(ora.go, function(x){
  x[x$pathway %in% go.immu.pathways, ]
})

ora.go.immu = lapply(ora.go.immu, function(x){
  x = x[x$overlap >= 4, ]
  if(nrow(x) > 1) {x}
})

ora.go.immu = ora.go.immu[grepl("_Late", names(ora.go.immu))]
names(ora.go.immu) = gsub("_Late", "", names(ora.go.immu))
names(ora.go.immu) = gsub("\n", " ", names(ora.go.immu))

ftrs.l.sub = ftrs.l.1[grepl("_Late", names(ftrs.l.1))]
names(ftrs.l.sub) = gsub("_Late", "", names(ftrs.l.sub))
names(ftrs.l.sub) = gsub("\n", " ", names(ftrs.l.sub))

ora.main =
  ora_bubble(
    gsea.res = ora.go.immu,
    ftrs.list = ftrs.l.sub,
    nbr.tops = 7,
    min.genes = 4,
    term.length = 60,
    sort.by.padj = F,
    facet.spit = F
  ) &
  theme(
    panel.grid.major.y = element_line(colour = "grey80", linetype = 2, linewidth = .2),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    legend.box.margin = margin(0, 160, 0, 0),
    legend.ticks.length = unit(0.05, 'cm'),
    legend.title = element_text(margin = margin(r = 5, l = 1, unit = "pt"), size = rel(1)),
    legend.margin = margin(t=-4, l = 5),
    legend.text = element_text(margin = margin(l = 1, t = 2, b = 2)),
    plot.title = element_text(margin = margin(0, 180, 4, 4), face = "plain"),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 6)
  ) &
  guides(
    fill = guide_colorbar(
      title = "Pathway\ndirection",
      title.vjust = .56,
      title.position = "left",
      barwidth = unit(5, 'lines'),
      barheight = unit(.35, 'lines'),
      order = 1,
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  ) &
  ggtitle("Late | Enrichment test for DE genes comparing non-CR with CR")


main.pl.a = plot_grid(
  NULL,
  plot_grid(
    plot_grid(
      dgea.pl.l + theme(
        axis.title.y = element_text(vjust = + 0),
        axis.text.y = element_text(size = rel(.9))
      ),
      NULL, nrow = 2, rel_heights = c(1, .01),
      labels = c("A", ""), label_fontface = "bold", label_size = 12, vjust = .8
    ),
    NULL,
    ora.main,
    ncol = 3, rel_widths = c(.658, .02, .32),
    labels = c("", "", "B"), label_fontface = "bold", label_size = 12,
    vjust = .8, hjust = 2.5
  ),
  nrow = 2, rel_heights = c(.01, 1)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DGEA | ORA | non-Cr vs. CR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora.go.immu = lapply(ora.go, function(x){
  x[x$pathway %in% unique(go.immu.pathways), ]
})

ora.go.immu = lapply(ora.go.immu, function(x){
  x = x[x$overlap >= 3, ]
  if(nrow(x) > 1) {x}
})

names(ora.go.immu) = gsub("CD4 ", "CD4\n", names(ora.go.immu))
names(ora.go.immu) = gsub("CD8 ", "CD8\n", names(ora.go.immu))
names(ftrs.l.1) = gsub("CD4 ", "CD4\n", names(ftrs.l.1))
names(ftrs.l.1) = gsub("CD8 ", "CD8\n", names(ftrs.l.1))

ora.supps =
  ora_bubble(
    gsea.res = ora.go.immu,
    ftrs.list = ftrs.l.1,
    nbr.tops = 5,
    min.genes = 4,
    term.length = 80,
    sort.by.padj = F,
    facet.spit = T,
    font.size = 7
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    strip.text.x = element_text(size = rel(1))
  ) &
  guides(
    fill = guide_colorbar(
      title = "Pathway\ndirection",
      title.vjust = 1,
      barwidth = unit(4, 'lines'),
      barheight = unit(.35, 'lines'),
      order = 1,
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )

ggsave2(
  filename="publication/extended_data_files/Fig_03_dgea_ora_cr_noncr.png",
  plot_grid(
    plot_grid(
      dgea.pl.lp, NULL, dgea.pl.vl, nrow = 3, rel_heights = c(1, .025, 1),
      labels = c("a", "", "b"), label_fontface = "bold", label_size = 12, vjust = 1.1
    ),
    NULL,
    ora.supps,
    nrow = 3, rel_heights = c(1.2, .1, 1.5),
    labels = c("", "", "c"), label_fontface = "bold", label_size = 12, vjust = 1.1
  ),
  width = 180, height = 200, dpi = 500, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Milestone graphs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ti_plots = function(
    se.dm = dm.cd8$se,
    sce.dm = dm.cd8$sce,
    res = "cluster",
    max.value = 1.65,
    swtch.nodes = NULL
){

  se.dm$cloneSize = as.character(se.dm$cloneSize)
  se.dm$cloneSize = as.factor(ifelse(
    grepl("Hyper", se.dm$cloneSize), "Hyperexpanded (100 > X)", se.dm$cloneSize
  ))
  se.dm$cloneSize = fct_relevel(se.dm$cloneSize, "Single (0 < X <= 1)", after = Inf)
  sce.dm$cloneSize = as.character(sce.dm$cloneSize)
  sce.dm$cloneSize = as.factor(ifelse(
    grepl("Hyper", sce.dm$cloneSize), "Hyperexpanded (100 > X)", sce.dm$cloneSize
  ))
  sce.dm$cloneSize = fct_relevel(sce.dm$cloneSize, "Single (0 < X <= 1)", after = Inf)

  dm.pd = cbind(
    se.dm@meta.data,
    Embeddings(se.dm, reduction = "DS")[, c(1:5)],
    slingPseudotime(sce.dm)
  )

  mst <- slingMST(sce.dm, as.df = TRUE)
  mst = mst[mst$Order == 1, ][1, ]
  dm.pl.1 =
    ggplot(dm.pd, aes(x = DS_1, y = DS_2, color = celltype)) +
    # scattermore::geom_scattermore(mapping=aes(x=DS_1, y=DS_2), pointsize = 5, color="black")+
    # scattermore::geom_scattermore(mapping=aes(x=DS_1, y=DS_2), pointsize = 4, color="white")+
    geom_point(shape = ".") +
    xlab("DC 1") +
    ylab("DC 2") +
    scale_color_manual(values = til.col) +
    guides(colour = guide_legend(
      title = NULL, override.aes = list(shape = 16, size = 2.75)
    )) +
    geom_path(
      data = slingCurves(
        sce.dm, as.df = TRUE) %>%
        dplyr::arrange(Order) %>%
        dplyr::rename("DS_1" = DS_1, "DS_2" = DS_2),
      aes(group = Lineage), col = "black", size = .3
    ) +
    # geom_path(data = mst, col = "black", size = .5)  +
    geom_point(
      data = mst, aes(col = "black"), color = "black", shape = 21,
      fill = "grey",  size = 2
    ) +
    theme(
      legend.position = c(0.85, 0.75),
      legend.key.size = unit(3, "mm"),
      legend.title = element_text( size=6),
      legend.text=element_text(size=6),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )


  clono.col = setNames(
    c("#003285", "#5AB2FF", "#A0DEFF", "#F08A5D", "#B83B5E", "#BBBBBB"),
    c(rev(levels(se.dm$cloneSize)), "NA")
  )
  tmp = dm.pd[order(dm.pd$clonalFrequency, decreasing = F, na.last = F), ]
  dm.pl.2 =
    ggplot(tmp, aes(x = DS_1, y = DS_2, color = cloneSize)) +
    geom_point(shape = ".") +
    xlab("DC 1") +
    ylab("DC 2") +
    scale_color_manual(values = clono.col) +
    guides(colour = guide_legend(
      title = "Clonotype group",
      override.aes = list(shape = 16, size = 2.75)
    )) +
    theme(
      legend.position = c(0.85, 0.79),
      legend.key.size = unit(3, "mm"),
      legend.title = element_text( size=6, vjust = -1),
      legend.text=element_text(size=6),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )

  ti.net = ti_networks(sce.dm, clust.res = res)
  nearest.cl = group_onto_nearest_cluster(ti.net)
  stopifnot(identical(names(nearest.cl), rownames(colData(sce.dm))))
  sce.dm$cluster = nearest.cl

  sce.dm.lp = sce.dm[, sce.dm$TIMEPOINT == "LP"]
  ti.net.lp = ti_networks(sds = sce.dm.lp, clust.res = res)
  top.ratio.lp = ti_topology_ratio(
    sds = sce.dm.lp, ti = ti.net.lp, target = "BEST_RESPONSE_CONSENSUS",
    max.value = max.value, permuations = 1000
  )

  sce.dm.l = sce.dm[, sce.dm$TIMEPOINT == "Late"]
  ti.net.l = ti_networks(sds = sce.dm.l, clust.res = res)
  top.ratio.l = ti_topology_ratio(
    sds = sce.dm.l, ti = ti.net.l, target = "BEST_RESPONSE_CONSENSUS",
    max.value = max.value, permuations = 1000
  )

  sce.dm.vl = sce.dm[, sce.dm$TIMEPOINT == "Very Late"]
  ti.net.vl = ti_networks(sds = sce.dm.vl, clust.res = res)
  top.ratio.vl = ti_topology_ratio(
    sds = sce.dm.vl, ti = ti.net.vl, target = "BEST_RESPONSE_CONSENSUS",
    max.value = max.value, permuations = 1000
  )

  ti.net = ti_networks(sce.dm, clust.res = res)

  top.ct.all = ti_topology_ct_comp(
    sds = sce.dm, ti = ti.net, target = "celltype", swtch.nodes = swtch.nodes
  )

  top.cl.all = ti_topology_ct_comp(
    sds = sce.dm, ti = ti.net, target = "cloneSize", fill.col = clono.col,
    swtch.nodes = swtch.nodes
  )

  list(
    dm.pl.1 = dm.pl.1,
    dm.pl.2 = dm.pl.2,
    top.ratio.lp = top.ratio.lp,
    top.ratio.l = top.ratio.l,
    top.ratio.vl = top.ratio.vl,
    top.ct.all = top.ct.all,
    top.cl.all = top.cl.all
  )
}

ti_cd4_plots = ti_plots(
  se.dm = dm.cd4$se, sce.dm = dm.cd4$sce, res = "cluster", max.value = 1.4,
  swtch.nodes = c(7,6)
)

ti_cd8_plots = ti_plots(
  se.dm = dm.cd8$se, sce.dm = dm.cd8$sce, res = "cluster", max.value = .88,
  swtch.nodes = c(6,5)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Plot Part | I
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c.th = theme(
  plot.margin = unit(c(0,0,10,0), "pt"),
  plot.title = element_text(hjust = 0.5, margin=margin(5,0,0,0))
)

clono.leg =
  plot_grid(
    get_legend(
      ti_cd4_plots$top.cl.all + labs(fill = "Clonotype group") +
        guides(fill = guide_legend(ncol = 2, title.position = "left", title.vjust = .95)) +
        theme(
          legend.key.size = unit(3, "mm"),
          legend.title = element_text( size=6),
          legend.text=element_text(size=6, margin = margin(l = 2)),
          #legend.position = c(0.45, 1.5)
        )
    )
  )

ct.leg.cd4 =
  plot_grid(
    get_legend(
      ti_cd4_plots$top.ct.all + labs(fill = NULL) +
        #guides(fill = guide_legend(ncol = 1, title.position = "left")) +
        theme(
          legend.key.size = unit(3, "mm"),
          legend.title = element_text( size=6),
          legend.text=element_text(size=6, margin = margin(l = 2))
        )
    )
  )

ct.leg.cd8 =
  plot_grid(
    get_legend(
      ti_cd8_plots$top.ct.all + labs(fill = NULL) +
        #guides(fill = guide_legend(ncol = 1, title.position = "left")) +
        theme(
          legend.key.size = unit(3, "mm"),
          legend.title = element_text( size=6),
          legend.text=element_text(size=6, margin = margin(l = 2))
        )
    )
  )

cd4.lim = c(
  max(ti_cd4_plots$top.ratio.l$data$y), min(ti_cd4_plots$top.ratio.l$data$y) - .11
)

ti_cd4_plots$top.ratio.vl + theme(legend.position = "none") +
  ggtitle("Very Late") + c.th + scale_y_reverse(limits = cd4.lim)

dm.cd4.pl =
  plot_grid(
    plot_grid(
      NULL,
      plot_grid(
        NULL,
        NULL,
        # ti_cd4_plots$dm.pl.1 +
        #   theme(
        #     legend.position = c(1.3, 0.25),
        #     plot.margin = unit(c(4,4,-10,4), "pt")
        #   ) + scale_y_reverse(),
        NULL,
        ti_cd4_plots$top.ct.all + theme(legend.position = "none") + scale_y_reverse(),
        ti_cd4_plots$top.cl.all + theme(legend.position = "none") + scale_y_reverse(),
        ncol = 5, rel_widths = c(0, 0, .1, 1, 1),
        labels = c("C", "", "", "", "D"), label_fontface = "bold", label_size = 11, vjust = 0.3
      ),
      NULL, nrow = 3, rel_heights = c(.1, 1, .1)
    ),
    NULL,
    plot_grid(
      plot_grid(
        ti_cd4_plots$top.ratio.lp + theme(legend.position = "none") +
          ggtitle("LP") + c.th + scale_y_reverse(limits = cd4.lim),
        ti_cd4_plots$top.ratio.l + theme(legend.position = "none") +
          ggtitle("Late") + c.th + scale_y_reverse(limits = cd4.lim),
        ti_cd4_plots$top.ratio.vl + theme(legend.position = "none") +
          ggtitle("Very Late") + c.th + scale_y_reverse(limits = cd4.lim),
        nrow = 1, scale = 1
      ),
      NULL,
      plot_grid(get_legend(ti_cd4_plots$top.ratio.lp), NULL, ncol = 2, rel_widths = c(2, -.4)),
      nrow = 3, rel_heights = c(1, .0, .1),
      labels = c("E", "", ""), label_fontface = "bold", label_size = 11
    ),
    NULL,
    ncol = 4, rel_widths = c(1.1, .05, 1.25, .01)
  )

cd8.lim = c(
  max(ti_cd8_plots$top.ratio.l$data$x), min(ti_cd8_plots$top.ratio.l$data$x) - .125
)

dm.cd8.pl =
  plot_grid(
    plot_grid(
      NULL,
      plot_grid(
        NULL,
        NULL,
        # ti_cd8_plots$dm.pl.1 +
        #   theme(
        #     legend.position = c(1.3, 0.25),
        #     plot.margin = unit(c(4,4,-10,4), "pt")
        #   ),
        NULL,
        ti_cd8_plots$top.ct.all + theme(legend.position = "none", plot.margin = unit(c(10,4,-15,4), "pt")) + coord_flip()  + scale_x_reverse(),
        ti_cd8_plots$top.cl.all + theme(legend.position = "none", plot.margin = unit(c(10,4,-15,4), "pt")) + coord_flip()  + scale_x_reverse(),
        ncol = 5, rel_widths = c(0, 0, .1, 1, 1)
      ),
      NULL, nrow = 3, rel_heights = c(.0, 1, .25)
    ),
    NULL,
    plot_grid(
      plot_grid(
        ti_cd8_plots$top.ratio.lp + theme(legend.position = "none") +
          ggtitle("LP") + c.th + coord_flip()  + scale_x_reverse(limits = cd8.lim),
        ti_cd8_plots$top.ratio.l + theme(legend.position = "none") +
          ggtitle("Late") + c.th + coord_flip()  + scale_x_reverse(limits = cd8.lim),
        ti_cd8_plots$top.ratio.vl + theme(legend.position = "none") +
          ggtitle("Very Late") + c.th + coord_flip()  + scale_x_reverse(limits = cd8.lim),
        nrow = 1, scale = 1
      ),
      NULL,
      plot_grid(get_legend(ti_cd8_plots$top.ratio.lp), NULL, ncol = 2, rel_widths = c(2, -.4)),
      nrow = 3, rel_heights = c(1, .0, .1)
    ),
    NULL,
    ncol = 4, rel_widths = c(1.1, .05, 1.25, .01)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Singletons vs Outcome
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.t@meta.data
pd$BEST_RESPONSE_CONSENSUS = as.character(pd$BEST_RESPONSE_CONSENSUS)
pd = pd[!is.na(pd$clonalFrequency), ]
singlletons = names(table(pd$CLONE_PSEUDO_ID)[table(pd$CLONE_PSEUDO_ID) == 1])
pd$clonalFrequency = ifelse(pd$CLONE_PSEUDO_ID %in% singlletons, "singletons", "expanded")
pd$celltype = pd$celltype_short_3

# Grundwert jeweils CD4 und CD8
pd.cl = lapply(as.character(unique(pd$celltype)), function(x){
  pd.sub = pd[pd$celltype == x, ]
  pd.sub = prop.table(table(pd.sub$orig.ident, pd.sub$clonalFrequency), margin = 1)
  pd.sub = as.data.frame.matrix(pd.sub)
  pd.sub = pd.sub[rowSums(is.na(pd.sub)) == 0, ]
  pd.sub$celltype = x
  pd.sub$orig.ident = rownames(pd.sub)
  pd.sub
}) %>% data.table::rbindlist()
pd.cl$singletons = pd.cl$singletons * 100

pd.cl.pl = merge(
  pd.cl, pd[!duplicated(pd$orig.ident), ], by = c("orig.ident"), all.x = T
)

boxpl = function(df = NULL, ctrst = "RESPONSE_CONSENSUS_2", leg.t = NULL) {

  ggplot(df, aes(TIMEPOINT, singletons, fill = .data[[ctrst]])) +
    geom_boxplot(
      outlier.colour = NA, position=position_dodge(0.8), fatten = 1,
      alpha = .7, lwd = 0.3
    ) +
    geom_dotplot(
      binaxis='y', stackdir='center', position=position_dodge(0.8),
      stroke = .5, dotsize = .6
    ) +
    ylim(0, max(pd.cl$singletons + 15)) +
    stat_compare_means(
      aes(group = .data[[ctrst]], label = ..p.format..), label.x = 2,
      hjust = .5, vjust = 1, size = 2.2
    ) +
    scale_fill_manual(
      values = setNames(c("#997700", "#6699CC"), c("non-CR", "CR"))
    ) +
    xlab(NULL) + ylab("% of singletons") +
    theme(
      panel.spacing = unit(.5, "lines"),
      axis.title.y = element_text(vjust = + 2),
      legend.margin=margin(t=5),
      legend.box.spacing = unit(0, "pt"),
      legend.position = "bottom"
    ) +
    guides(
      fill = guide_legend(title = NULL, title.position="left", title.hjust = 0.5)
    ) +
    facet_wrap( ~ celltype.x) +
    scale_y_continuous(
      breaks=c(0, 25, 50, 75, 100),
      limits = c(0, 120)
    )
}

boxpl(pd.cl.pl, "BEST_RESPONSE_CONSENSUS", leg.t = "Response")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | STARtac
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
startrac = function(obj, tp = "Late", celltype = "celltype_short_3"){

  pd = obj@meta.data

  pd = pd[pd$TIMEPOINT == tp, ]
  pd$barcode = rownames(pd)
  pd$celltype = pd[[celltype]]
  pd = droplevels(pd)
  pd = pd %>% dplyr::select(
    "Cell_Name" = barcode, "clone.id" = CLONE_PSEUDO_ID, "patient" = orig.ident,
    "majorCluster" = celltype, "group" = BEST_RESPONSE_CONSENSUS
  )
  pd$loc = "T"
  pd$patient = as.character(pd$patient)

  suppressWarnings({
    suppressMessages({
      out <- Startrac.run(pd, proj="bcma", cores=NULL,verbose=F)
    })
  })

  # plot(out,index.type="pairwise.tran",byPatient=T)
  plot(out,index.type="cluster.all",byPatient=T)

  df = out@cluster.data
  df = df[df$aid != "bcma", ]
  df$NCells = NULL; df$migr = NULL
  df = reshape2::melt(df, keys = c("aid", "majorCluster"))
  df = df[df$variable == "expa", ]
  df = df[!is.na(df$value), ]
  df
}

star.lp = startrac(obj = se.t, tp = "LP", celltype = "celltype_short_3")
star.l = startrac(obj = se.t, tp = "Late", celltype = "celltype_short_3")
star.vl = startrac(se.t, tp = "Very Late", celltype = "celltype_short_3")

star = rbind(star.lp, star.l, star.vl)
star = merge(
  star, se.t@meta.data[!duplicated(se.t$orig.ident), ],
  by.x = "aid", by.y = "orig.ident"
)
star$value[star$value < 0] = 0

shannon.pl =
  ggplot(star, aes(TIMEPOINT, value, fill = BEST_RESPONSE_CONSENSUS)) +
  geom_boxplot(outlier.colour = NA, position=position_dodge(0.8), fatten = 1, alpha = .7, lwd = 0.3) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), stroke = .5, dotsize = .6) +
  stat_compare_means(
    aes(group = BEST_RESPONSE_CONSENSUS, label = ..p.format..), label.x = 2, hjust = .5, vjust = 1,
    size = 2.2
  ) +
  scale_fill_manual(values = c("CR" = "#6699CC", "non-CR" = "#997700")) +
  # ylim(0, max(pd.cl$freq + .1)) +
  scale_y_sqrt(
    #breaks = c(0, .05, .2, .45),
    limits = c(0, round(max(star$value), digits = 1) + .1)
  ) +
  xlab(NULL) + ylab("Clonality") +
  theme(
    panel.spacing = unit(.5, "lines"),
    axis.title.y = element_text(vjust = + 2),
    legend.margin=margin(t=-5),
    legend.box.spacing = unit(0, "pt"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "", title.position="top", title.hjust = 0.5)) +
  facet_wrap( ~ majorCluster)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Final plot | II
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cde.png = ggdraw() +
  cowplot::draw_image(
    "publication/figures_main/fig_03_cde.png", scale = 1.15)

main.pl.b = plot_grid(
  cde.png + theme(plot.margin = unit(c(0,0,-27, 33), "pt")),
  # plot_grid(dm.cd4.pl, NULL, dm.cd8.pl, nrow = 3, rel_heights = c(1, .15, 1)),
  NULL,
  plot_grid(
    NULL,
    boxpl(pd.cl.pl, "BEST_RESPONSE_CONSENSUS", leg.t = "Response"),
    NULL,
    shannon.pl,
    nrow = 4, rel_heights = c(0, 1, .05, 1),
    labels = c("F", "", "G"), label_fontface = "bold", label_size = 12
  ),
  ncol = 3, rel_widths = c(.618, 0.05, .382)
)

ggsave2(
  filename= "publication/figures_main/fig_03_2_legend.png",
  plot_grid(ct.leg.cd4, clono.leg, ct.leg.cd8, nrow = 3),
  width = 180, height = 60, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | CD4/CD8 gene enrichment for milestones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# CD4
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cust.gene.sets = mir[grepl("Chu", mir$TermID), ]
cust.gene.sets = cust.gene.sets[grepl("CD4", cust.gene.sets$TermID), ]
cust.gene.sets$TermID = gsub("CD4-", "", cust.gene.sets$TermID)
cust.gene.sets$TermID = gsub(" .+", "", cust.gene.sets$TermID)
cust.gene.sets$GeneID = gsub("MT\\.", "MT-", cust.gene.sets$GeneID)
cust.gene.sets = split(cust.gene.sets, cust.gene.sets$TermID)
cust.gene.sets = lapply(cust.gene.sets, function(x){x$GeneID})

se.gs.cd4 = ucell_enrich(se.w = dm.cd4$se, cust.gene.sets)

pd = se.gs.cd4@meta.data
bc = lapply(levels(pd$TIMEPOINT), function(i){
  print(i)
  pd.s = pd[pd$TIMEPOINT == i, ]
  print(table(pd.s$BEST_RESPONSE_CONSENSUS, pd.s$nearest.cl))
  set.seed(1234)
  bc  = pd.s %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::group_by(nearest.cl, BEST_RESPONSE_CONSENSUS) %>%
    dplyr::slice_sample(n = 500) %>%
    dplyr::pull("barcode")
  bc
})

se.gs.cd4@meta.data = se.gs.cd4@meta.data %>%
  dplyr::mutate(
    nearest.cl.2 = dplyr::case_when(
      nearest.cl == "7" ~ "6",
      nearest.cl == "6" ~ "7",
      TRUE ~ nearest.cl
    )
  )

se.gs.cd4$nearest.cl.2 = factor(
  se.gs.cd4$nearest.cl.2,
  levels = naturalsort(unique(as.character(se.gs.cd4$nearest.cl.2)))
)

hm.cd4 =
  enrich_heatmap_milestones(
    se.w = subset(se.gs.cd4, cells = unlist(bc)),
    meta.vars = c("TIMEPOINT", "nearest.cl.2")
  )
hm.cd4 = hm.cd4 + ggtitle("CD4")

ti.net.cd4 = ti_networks(sds = dm.cd4$sce, clust.res = "cluster")
milestones.cd4 = ti_topology_ct_comp(
  sds = dm.cd4$sce, ti = ti.net.cd4, target = "celltype", swtch.nodes = c(7,6)
) + theme(legend.position = "bottom") +  scale_y_reverse()

milestones.cd4 = milestones.cd4 +
  guides(fill = guide_legend(title = NULL, ncol = 3)) +
  theme(
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(t = -5),
    legend.key.spacing.y = unit(1, 'pt'),
    legend.key.spacing.x = unit(7, 'pt'),
    legend.text = element_text(margin = margin(l = 2))
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CD8
cust.gene.sets = mir[grepl("Chu", mir$TermID), ]
cust.gene.sets = cust.gene.sets[grepl("CD8", cust.gene.sets$TermID), ]
cust.gene.sets$TermID = gsub("CD8-", "", cust.gene.sets$TermID)
cust.gene.sets$TermID = gsub(" .+", "", cust.gene.sets$TermID)
cust.gene.sets$GeneID = gsub("MT\\.", "MT-", cust.gene.sets$GeneID)
cust.gene.sets = split(cust.gene.sets, cust.gene.sets$TermID)
cust.gene.sets = lapply(cust.gene.sets, function(x){x$GeneID})

se.gs.cd8 = ucell_enrich(se.w = dm.cd8$se, cust.gene.sets)

pd = se.gs.cd8@meta.data
bc = lapply(levels(pd$TIMEPOINT), function(i){
  print(i)
  pd.s = pd[pd$TIMEPOINT == i, ]
  print(table(pd.s$BEST_RESPONSE_CONSENSUS, pd.s$nearest.cl))
  set.seed(1234)
  bc  = pd.s %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::group_by(nearest.cl, BEST_RESPONSE_CONSENSUS) %>%
    dplyr::slice_sample(n = 2000) %>%
    dplyr::pull("barcode")
  bc
})

se.gs.cd8@meta.data = se.gs.cd8@meta.data %>%
  dplyr::mutate(
    nearest.cl.2 = dplyr::case_when(
      nearest.cl == "6" ~ "5",
      nearest.cl == "5" ~ "6",
      TRUE ~ nearest.cl
    )
  )


hm.cd8 =
  enrich_heatmap_milestones(
    se.w = subset(se.gs.cd8, cells = unlist(bc)),
    meta.vars = c("TIMEPOINT", "nearest.cl.2")
  )
hm.cd8 = hm.cd8 + ggtitle("CD8")

ti.net.cd8 = ti_networks(dm.cd8$sce, clust.res = "cluster")
milestones.cd8 = ti_topology_ct_comp(
  sds = dm.cd8$sce, ti = ti.net.cd8, target = "celltype", swtch.nodes = c(6,5)
) + labs(fill = NULL) + theme(legend.position = "bottom") + coord_flip()  + scale_x_reverse()
milestones.cd8 =
  milestones.cd8 +
  guides(fill = guide_legend(title = NULL, ncol = 4)) +
  theme(
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(l = 25, t = -5),
    legend.key.spacing.y = unit(1, 'pt'),
    legend.key.spacing.x = unit(7, 'pt'),
    legend.text = element_text(margin = margin(l = 2))
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supp | DGEA and ORA (GO terms) for milestones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_ora_cd4 = dgea_gex(
  obj = dm.cd4$se, group = "BEST_RESPONSE_CONSENSUS", ctrs.grp1 = "non-CR", ctrs.grp2 = "CR",
  .target = "nearest.cl",
  latent.vars = c("STUDY", "nFeature_RNA"),
  split.by.tp = T, logfc.threshold = log2(1.25), min.pct = .2,
  subsample = T, threads = 15, min.de.genes = 1
)

dgea_ora_cd8 = dgea_gex(
  obj = dm.cd8$se, group = "BEST_RESPONSE_CONSENSUS", ctrs.grp1 = "non-CR", ctrs.grp2 = "CR",
  .target = "nearest.cl",
  latent.vars = c("STUDY", "nFeature_RNA"),
  split.by.tp = T, logfc.threshold = log2(1.25), min.pct = .2,
  subsample = T, threads = 15, min.de.genes = 1
)

dg = dgea_ora_cd4$res.dgea.sign
ftrs.l = dg
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
names(ftrs.l)
ftrs.l.cd4 = ftrs.l

ora.cd4 = parallel::mclapply(ftrs.l.cd4, function(x){
  run_nmf_ora(genes = x, universe = rownames(se.t), category = "CD4")
}, mc.cores = 1)

dg = dgea_ora_cd8$res.dgea.sign
ftrs.l = dg
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
names(ftrs.l)
ftrs.l.cd8 = ftrs.l

ora.cd8 = parallel::mclapply(ftrs.l.cd8, function(x){
  run_nmf_ora(genes = x, universe = rownames(se.t), category = "CD8")
}, mc.cores = 1)

a =
  ora_bubble(
    gsea.res = ora.cd4,
    ftrs.list = ftrs.l.cd4,
    nbr.tops = 6,
    min.genes = 3,
    # max.value = .75,
    # quantile.filter = T
    term.length = 60,
    sort.by.padj = F,
    facet.spit = T
  ) +
  theme(
    axis.text.y = element_text(size = 8, lineheight = .75),
    panel.grid.major = element_line(colour = "grey80", linetype = 2, linewidth = .2)
  ) +
  ggtitle("CD4: non-CR vs. CR") +
  guides(
    fill = guide_colorbar(
      title =  "Pathway\ndirection",
      barwidth = unit(.5, 'lines'),
      barheight = unit(3, 'lines'),
      order = 1, ticks.linewidth = .75/.pt,  frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )


b =
  ora_bubble(
    gsea.res = ora.cd8,
    ftrs.list = ftrs.l.cd8,
    nbr.tops = 5,
    min.genes = 3,
    # max.value = .75,
    # quantile.filter = T,
    term.length = 60,
    sort.by.padj = F,
    facet.spit = T
  ) +
  theme(
    axis.text.y = element_text(size = 8, lineheight = .75),
    panel.grid.major = element_line(colour = "grey80", linetype = 2, linewidth = .2)
  ) +
  ggtitle("CD8: non-CR vs. CR") +
  guides(
    fill = guide_colorbar(
      title =  "Pathway\ndirection",
      barwidth = unit(.5, 'lines'),
      barheight = unit(3, 'lines'),
      order = 1, ticks.linewidth = .75/.pt,  frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot

ggsave2(
  filename="publication/extended_data_files/Fig_03_milestones.png",
  plot_grid(
    plot_grid(
      plot_grid(
        NULL,
        plot_grid(NULL, milestones.cd4 , nrow = 2, rel_heights = c(.075, 1)),
        NULL, hm.cd4, nrow = 1, rel_widths = c(.05, 1, .1, 3)
      ),
      NULL,
      plot_grid(
        NULL,
        plot_grid(NULL, milestones.cd8 , nrow = 2, rel_heights = c(.075, 1)),
        NULL, hm.cd8, nrow = 1, rel_widths = c(.05, 1, .1, 3)
      ),
      nrow = 3, rel_heights = c(1, .1, 1),
      labels = c("a", "", "b"), label_fontface = "bold", label_size = 11
    ),
    NULL,
    plot_grid(
      a, NULL, b, ncol = 3, rel_widths  = c(1, .1, 1),
      labels = c("c", "", "d"), label_fontface = "bold", label_size = 11
    ),
    nrow = 3, rel_heights = c(2.5, .1, 1)
  ),
  width = 180, height = 130, dpi = 500,
  bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Potential tumor reactive T cells (pTRTs)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
reactive_geneset = list("sigMM"=c(
  "GNLY","ZNF683","GZMH","FGFBP2","GZMB","NKG7","CCL5","HOPX","KLRD1","EFHD2",
  "CD8A","CTSW","CST7","ITGB1","BHLHE40"
))

se.tcr = ucell_enrich(
  se.w = se.t,
  reactive_geneset, ncores = 30
)

DefaultAssay(se.tcr) = "custom_UCell_score"
se.tcr = AddMetaData(se.tcr, FetchData(se.tcr, vars = "sigMM"))

###

ptrt.es.violin =
  ggplot(se.tcr@meta.data, aes(x = celltype, y = sigMM, fill = celltype)) +
  ggrastr::geom_jitter_rast(
    #position = position_dodge(width = .8),
    width = .08,
    alpha = 1, show.legend = F, color='#555555', shape = ".", raster.dpi = 200, scale = .25
  ) +
  geom_violin(linewidth = 0.2, alpha = .8, scale = "width") +
  stat_summary(
    fun= "median",geom = "crossbar", width = 0.2, show.legend = F,
    lwd = .3, position = position_dodge(width = 0.8)
  )+
  scale_fill_manual(values = til.col) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05, size = rel(.9)),
    axis.ticks.x  =  element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = rel(1)),
    axis.title.y = element_text(vjust = + 2),
    plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1)),
  ) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(-0.0001, 1)) +
  ylab("Enrichment\nfor pTRTs") +
  labs(fill = NULL) +
  ggtitle("Potential Tumor Reactive T Cells (pTRTs)")

####################

pd = se.tcr@meta.data
pd$pTRT_DICHO = as.factor(ifelse(
  as.numeric(scale(pd$sigMM)) > median(as.numeric(scale(pd$sigMM))),
  "pTRTs", "non-pTRTs"
))
# pd = pd[!is.na(pd$clonalFrequency), ]
pd$celltype = pd$celltype

pd.cl = lapply(as.character(unique(pd$celltype)), function(x){
  pd.sub = pd[pd$celltype == x, ]
  pd.sub = prop.table(table(pd.sub$orig.ident, pd.sub$pTRT_DICHO), margin = 1)
  pd.sub = as.data.frame.matrix(pd.sub)
  pd.sub[is.na(pd.sub)] = 0
  pd.sub$celltype = x
  pd.sub$orig.ident = rownames(pd.sub)
  pd.sub
}) %>% data.table::rbindlist()

ct.fltr = pd.cl %>% dplyr::group_by(celltype) %>%
  dplyr::summarise(counts = sum(pTRTs > 0.5, na.rm = TRUE)) %>%
  dplyr::filter(counts > 10) %>%
  dplyr::pull(celltype)

pd.cl = pd.cl[pd.cl$celltype %in% ct.fltr, ]
pd.cl = merge(pd.cl, pd[!duplicated(pd$orig.ident), ], by = c("orig.ident"), all.x = T)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Survival analysis with cell type fractions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd.tp = split(pd.cl, pd.cl$TIMEPOINT)

# file = pd.tp$Late

cox.ct.res = parallel::mclapply(pd.tp, function(file) {

  TP = as.character(file$TIMEPOINT[1])
  pd = droplevels(file)

  prop = pd %>% dplyr::select(orig.ident, celltype.x, pTRTs) %>%
    reshape2::dcast(celltype.x  ~ orig.ident)
  rownames(prop) = prop$celltype.x; prop$celltype.x = NULL
  prop[is.na(prop)] = 0
  prop = as.matrix(prop)
  prop = prop[matrixStats::rowSums2(prop, na.rm = T) != 0, ]

  pheno = pd
  pheno = pheno[!duplicated(pheno$orig.ident), ] %>% data.frame()
  rownames(pheno) = pheno$orig.ident
  pheno = pheno[base::intersect(rownames(pheno), colnames(prop)), ]
  prop = prop[, base::intersect(rownames(pheno), colnames(prop))]
  stopifnot(identical(rownames(pheno), colnames(prop)))

  # Survival object
  survivalDOD = Surv(pheno$PFS, pheno$PROGRESSION)
  survModel = survfit(survivalDOD ~ 1, data = pheno)

  run_cox = function(obj) {

    exprs.mat.spl = split(as.matrix(obj), rownames(obj))

    # celltype = exprs.mat.spl$CD8.EMRA.2
    coxModel.l = parallel::mclapply(exprs.mat.spl, function(celltype) {
      res = summary(survival::coxph(survivalDOD ~ celltype + STUDY, data = pheno))
      res = c(
        "HR" = res$coef[1, , drop = F][2],
        "logHR" = res$coef[1, , drop = F][1],
        "SE_logHR" = res$coef[1, , drop = F][3],
        "L95CI" = res$conf.int[1, , drop = F][,"lower .95"],
        "U95CI" = res$conf.int[1, , drop = F][,"upper .95"],
        "Rsquare" = res$rsq[[1]],
        "Pval" = res$coef[1, , drop = F][5],
        "PH" =  unname(cox.zph(survival::coxph(survivalDOD ~ celltype))$table[, 3][2])
      )
    }, mc.cores = 1)

    coxModel = data.frame(do.call("rbind", coxModel.l))
    coxModel$Pval_adj = p.adjust(coxModel$Pval, method = "BH")
    coxModel
  }

  # Cox survival model fitting
  coxModel = run_cox(obj = prop)
  coxModel$TIMEPOINT = TP
  coxModel = coxModel[order(coxModel$Pval, decreasing = F), ]
  prop = prop[rownames(coxModel), ]

  list(coxModel = coxModel, prop = prop, pheno = pheno)

}, mc.cores = 1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | forest, logHR + CI
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = lapply(names(cox.ct.res), function(x){
  cox = cox.ct.res[[x]]$coxModel
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
df$Cox_Pval_pl = ifelse(df$Pval < 0.1, paste0(df$Cox_Pval_pl, " *"), df$Cox_Pval_pl)

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
ggplot(df, aes(x = CELLTYPE, y = HR, ymin = L95CI, ymax = U95CI, col = TIMEPOINT, fill=TIMEPOINT)) +
  geom_hline(yintercept=1, lty=2, lwd = .2, col = "#555555") +
  geom_linerange(linewidth=1.5, position=position_dodge(width = 0.8)) +
  facet_wrap(~ CELLTYPE, ncol = 1, scales = "free_y") +
  geom_point(size=2.5, shape=21, colour="white", stroke = 1, position=position_dodge(width = 0.8)) +
  geom_text(
    aes(y = U95CI + 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)),
    color = "black",
    position=position_dodge(width = .8),
    hjust = -0.1, vjust = 0.5, size = rel(2.25)
  ) +
  coord_flip() +
  scale_colour_manual(values = alpha(barCOLS, .7))+
  scale_fill_manual(values = dotCOLS)+
  scale_y_log10(
    name= "Hazard ratio (90% CI)",
    limits = c(min(df$L95CI), max(df$U95CI) + 2000),
    breaks = c(0.001, 0.1, 1, 10),
    labels = c(0.001, 0.1, 1, 10)
    # labels = function(x) sprintf("%g", x)
  ) +
  xlab(NULL) +
  theme(
    strip.text.x  = element_blank(),
    # panel.spacing  = unit(.25, "lines"),
    legend.title=element_blank(),
    legend.position = "right",
    axis.ticks.y = element_blank(),
    legend.box.spacing = unit(0, "pt"),
    legend.text = element_text(margin = margin(l = 2, unit = "pt"), size = rel(1))
  ) +
  scale_x_discrete(expand = c(0.5, 0))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | Barplot: cell type composition for singletons and expanded clones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.t@meta.data
pd = pd[!is.na(pd$clonalFrequency), ]
singlletons = names(table(pd$CLONE_PSEUDO_ID)[table(pd$CLONE_PSEUDO_ID) == 1])
pd$clonalFrequency = ifelse(pd$CLONE_PSEUDO_ID %in% singlletons, "singletons", "expanded")
# pd$clonalFrequency = ifelse(pd$clonalFrequency > 1, "expanded", "singletons")

pd.ct = pd %>%
  group_by(clonalFrequency, CAR_BY_EXPRS, TIMEPOINT, celltype_short_3, celltype) %>%
  dplyr::summarise(nbr.cells = n()) %>%
  dplyr::mutate(freq = nbr.cells / sum(nbr.cells)) %>%
  data.frame()

pd.ct$clonalFrequency = gsub("singletons", "Singletons", pd.ct$clonalFrequency)
pd.ct$clonalFrequency = gsub("expanded", "Expanded clones", pd.ct$clonalFrequency)
pd.ct$freq = pd.ct$freq * 100
pd.ct$clonalFrequency = factor(pd.ct$clonalFrequency, levels = c("Singletons", "Expanded clones"))
pd.ct$CAR_BY_EXPRS = ifelse(pd.ct$CAR_BY_EXPRS == TRUE, "CAR+", "CAR-")
pd.ct = droplevels(pd.ct)

barpl = function(df = NULL, lin = "CD4 T-Cell", title = NULL) {

  lbl = unique(df$celltype)
  lbl = setNames(gsub("CD4.|CD8.", "", lbl), lbl)

  if(is.null(title)) {
    title = lin
  }

  ggplot(df[df$celltype_short_3 == lin, ], aes(x = TIMEPOINT, y = freq, fill = celltype)) +
    geom_bar(stat="identity") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(margin = margin(l = 1.5, unit = "pt"), size = 6),
      legend.position = "bottom",
      legend.margin = margin(t=-5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1.1, size = rel(1)),
      ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
      panel.spacing = unit(.5, "lines")
    ) +
    facet_nested(~ clonalFrequency+ CAR_BY_EXPRS, scales = "free", space = "free") +
    scale_fill_manual(values = til.col, labels = lbl) +
    ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
    ggtitle(title) +
    geom_text(
      aes(label = ifelse(freq > 10, paste0(round(freq,1)), '')),
      position = position_stack(vjust = 0.5), size = 2
    ) +
    guides(fill = guide_legend(title = NULL, nrow = 1))
}

# barpl(df = pd.ct) / barpl(pd.ct, lin = "CD8 T-Cell", title = "Composition of CD8 subtypes")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main | DGEA | CAR+ Late vs LP; CAR+ vs. CAR-
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_cars = function(
    se.w = se.t,
    group = c("Late_TRUE", "LP_FALSE"),
    contrast.group = "TIMEPOINT",
    contrast = c("Late", "LP"),
    min.pct = .2,
    logfc.threshold = log2(1.25),
    .min.cells = 5,
    celltype = "celltype_short_3"
){

  se.w = se.w[ ,grepl("CD4|CD8", se.w$celltype_short_3)]

  se.w$TIMEPOINT_CAR = paste0(se.w$TIMEPOINT, "_", se.w$CAR_BY_EXPRS)

  se.w = se.w[, se.w$TIMEPOINT_CAR %in% group]
  se.w@meta.data[[celltype]] = factor(gsub(" T-Cell", "", se.w@meta.data[[celltype]]))
  se.w@meta.data = droplevels(se.w@meta.data)
  se.w@meta.data$TMP = paste0(se.w$orig.ident, "_", se.w@meta.data[[celltype]], "_", se.w$CAR_BY_EXPRS)

  k = table(paste0(se.w$PATIENT_ID, "_", se.w@meta.data[[celltype]], "_", se.w$CAR_BY_EXPRS)) < 10
  k = k[k]
  se.w = se.w[, !paste0(se.w$PATIENT_ID, "_", se.w@meta.data[[celltype]]) %in% gsub("_TRUE|_FALSE", "", names(k))]

  k = table(paste0(se.w$PATIENT_ID, "_", se.w@meta.data[[celltype]], "_", se.w$CAR_BY_EXPRS))
  k = table(gsub("_TRUE|_FALSE", "", names(k))) == 2
  k = k[k]
  se.w = se.w[, paste0(se.w$PATIENT_ID, "_", se.w@meta.data[[celltype]]) %in% names(k)]
  se.w@meta.data = droplevels(se.w@meta.data)

  # se.w = se.w[, se.w$PATIENT_ID %in% rownames(df[df$Late_TRUE >= 50, ])]
  Idents(se.w) = "TMP"
  se.w = subset(se.w, downsample = 100)
  Idents(se.w) = "orig.ident"
  se.w@meta.data = droplevels(se.w@meta.data)

  pd = se.w@meta.data
  pd = droplevels(pd[!duplicated(pd$orig.ident), ])
  df = table(se.w$PATIENT_ID, paste0(se.w$TIMEPOINT, "_", se.w$celltype_short_3))
  df = data.frame(rbind(df))
  df$RESP = pd[["BEST_RESPONSE_CONSENSUS"]][match(rownames(df), pd$PATIENT_ID)]
  df = df[order(df$RESP), ]
  print(df)

  res = run_wilx(
    obj = se.w,
    target = celltype, min.cells = .min.cells, rm.var.chains = T,
    lfc.thresh = logfc.threshold, min.pct.thres = min.pct, assay = "RNA",
    contrast.group = contrast.group, contrast = contrast,
    test.method = "MAST", latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA")
  )
  res
}
ctrst.1 = dgea_cars()
ctrst.2 = dgea_cars(
  contrast.group = "CAR_BY_EXPRS", group = c("Late_TRUE", "Late_FALSE"),
  contrast = c("TRUE", "FALSE")
)
# table(ctrst.1$cluster, ctrst.1$significant)
# table(ctrst.2$cluster, ctrst.2$significant)

ctrst.1.sign = subset(ctrst.1, significant == T)
ctrst.2.sign = subset(ctrst.2, significant == T)

ctrst.sign = rbind(ctrst.1.sign, ctrst.2.sign)

###

stem = c("HNRNPLL", "CCR7", "SELL", "LEF1", "IL7R", "TCF7", "KLF3", "BTG2")
act = c("CD28", "CD27", "GZMK", "GZMA", "PRF1", "CCL4", "GZMM", "CCL5", "NKG7", "GZMB", "GNLY", "LYAR", "CXCR3", "CXCR4", "TXNIP")
pre_exh = c("TIGIT", "ENTPD1", "HAVCR2", "CTLA4", "LAG3", "PDCD1", "EOMES", "TOX", "PMCH", "BATF", "PRDM1", "CX3CR1", "LYST")
res = c("MCM5", "TNFRSF9", "ITGAE", "CD69", "NR4A2")
df.ftrs = rbind(
  data.frame(feature = stem, type = "Stem-like,\nmemory"),
  data.frame(feature = act, type = "Activation, cytotoxicity,\neffector function"),
  data.frame(feature = pre_exh, type = "Exhaustion")
  # data.frame(feature = res, type = "Resident")
)
res.wlx.pl = ctrst.sign[ctrst.sign$feature %in% df.ftrs$feature, ]
res.wlx.pl$Type = df.ftrs$type[match(res.wlx.pl$feature, df.ftrs$feature)]
res.wlx.pl$group.1 = factor(
  ifelse(res.wlx.pl$group.1 == "Late", "CAR+ vs. LP", "CAR+ vs. CAR-"),
  levels = c("CAR+ vs. LP", "CAR+ vs. CAR-")
)

top.de =
  ggplot(res.wlx.pl, aes(avg_log2FC, feature, fill = cluster)) +
  geom_bar(stat="identity", width = .8) +
  facet_grid2(
    Type ~ group.1,
    scales = "free_y", independent = F, space = "free"
  ) +
  geom_vline(xintercept = 0, linewidth = .2, linetype = "dashed") +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(3, "mm"),
    panel.spacing = unit(.5, "lines"),
    legend.margin = margin(t=-2.5),
    legend.position = "bottom",
    legend.title = element_blank(),

    strip.text.x = element_text(margin = margin(0, 4.4, 4.4, 4.4, "pt"), face = "plain"),
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1))
  ) +
  labs(x = "log2 fold change") +
  scale_fill_manual(values = c("#555555", "#BBBBBB"))
# ggtitle("Late | CAR+ vs. LP")

# Supp Table
sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info/table_de_car_vs_lp.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  ctrst.1.sign %>% dplyr::select(
    Gene_name = feature, log2FC = avg_log2FC, Pvalue = p_val, FDR = p_val_adj,
    Time_point = group.1, Cell_identity = cluster
  ) %>%
    dplyr::arrange(Time_point, Cell_identity, FDR),
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)

# Supp Table
sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info/table_de_car_vs_noncar.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  ctrst.2.sign %>% dplyr::select(
    Gene_name = feature, log2FC = avg_log2FC, Pvalue = p_val, FDR = p_val_adj,
    Time_point = group.1, Cell_identity = cluster
  ) %>%
    dplyr::arrange(Time_point, Cell_identity, FDR),
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final Plot | Part III
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
main.pl.c = plot_grid(
  NULL,
  plot_grid(
    NULL,
    plot_grid(
      barpl(pd.ct, title = "Composition of CD4 subtypes"),
      NULL,
      barpl(pd.ct, lin = "CD8 T-Cell", title = "Composition of CD8 subtypes"),
      nrow = 3, rel_heights = c(1, .1, 1)
    ),
    NULL,
    plot_grid(NULL, top.de, nrow = 2, rel_heights = c(.0, 1)),
    NULL,
    plot_grid(
      plot_grid(
        plot_grid(NULL, ptrt.es.violin, nrow = 2, rel_heights = c(.0, 1)),
        NULL, rel_widths = c(1, .1)
      ),
      NULL,
      hr.pl, nrow = 3, rel_heights = c(1.2, .05, 2),
      labels = c("", "", "K"), label_fontface = "bold", label_size = 12, vjust = .5
    ),
    ncol = 6, rel_widths = c(.01, 1, .15, 1, .1, 1.25),
    labels = c("H", "", "", "I", "", "J"), label_fontface = "bold",
    label_size = 12, vjust = .6
  ),
  nrow = 2, rel_heights = c(.015, 1)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="publication/figures_main/fig_03.png",
  plot_grid(
    main.pl.a, NULL, main.pl.b, NULL, main.pl.c,
    nrow = 5, rel_heights = c(.9, .05, .9, .15, 1)
  ),
  width = 165, height = 200, dpi = 500, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

set.seed(1243)
ggsave2(
  filename="publication/figures_main/fig_03.pdf",
  plot_grid(
    main.pl.a, NULL, main.pl.b, NULL, main.pl.c,
    nrow = 5, rel_heights = c(.9, .05, .9, .15, 1)
  ),
  width = 165, height = 200, dpi = 500, bg = "white", units = "mm", scale = 1.6,
  device = cairo_pdf
)



