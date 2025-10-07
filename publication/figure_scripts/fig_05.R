print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Figure 5")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")


Sys.setenv( TZ="Etc/GMT+1" )

.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggh4x", "ggrepel", "scico", "DescTools",
  "scCustomize", "ggpubr", "immunarch", "Ckmeans.1d.dp", "PieGlyph", "ggstats"
)
.bioc_packages = c(
  "dittoSeq", "SummarizedExperiment", "limma", "impute"
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

if (any(!"speckle" %in% installed.packages())) {
  remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE,  dependencies = "Suggest")
}
library(speckle)

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/ora.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.t = readRDS(
  paste0(manifest$meta_pub$work, "integration/06_seurat_harmony_t_all_new.Rds")
)
vdj.t = readRDS(
  paste0(manifest$meta_pub$work, "integration/05_vdj_t_new.Rds")
)

se.t = AddMetaData(se.t, vdj.t)

mir = readRDS(
  "data/signatures/mir_gene_set_collection.Rds"
)

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

se.meta = readRDS(
  paste0(manifest$meta_pub, "integration/05_seurat_harmony_all_new.Rds")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FACs: CAR+ and CD4/CD8 Ratio
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata = readRDS("publication/clinicial_data/clinical_table_DF_2024_10_28.Rds")
pdata.clin = pdata$pdata.clin

pdata.imm = pdata$pdata.imm.s
pdata.imm$CAR_RATIO = pdata.imm$CD4_CAR_PERC / pdata.imm$CD8_CAR_PERC

plot_expansion = function(
    obj.imm = imm.s,
    obj.pd = pdata.clin,
    group = "PRODUCT",
    palette = c(cilta = "#004488", ide = "#44AA99"),
    labels = c("cilta-cel", "ide-cel"),
    leg.pos = "bottom"
) {

  pdata.exp = obj.pd
  pdata.exp$GROUP = pdata.exp[[group]]

  imm.exp = obj.imm %>%
    dplyr::select(SAMPLE_ID, DAY, CD3_CAR_PERC, CAR_RATIO) %>%
    dplyr::filter(DAY != "Leukapheresis") %>%
    dplyr::left_join(pdata.exp[,c("SAMPLE_ID", "GROUP")]) %>%
    dplyr::mutate(CAR_RATIO = ifelse(DAY == "Day 0", 0, CAR_RATIO)) %>%
    tidyr::pivot_longer(c("CD3_CAR_PERC", "CAR_RATIO")) %>%
    dplyr::mutate(
      DAY = factor(DAY, levels = c("Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))
    )

  p.wlx = imm.exp %>%
    dplyr::filter(DAY != "Day 0") %>%
    dplyr::group_by(DAY, name) %>%
    rstatix::wilcox_test(value ~ GROUP) %>%
    rstatix::add_significance(p.col = "p", output.col = "p_signif")
  p.wlx$p_signif[p.wlx$p_signif == "ns"] = ""

  l = list()
  for (i in 1:length(p.wlx$p)) {
    if (p.wlx$p[i] <= 0.0001) {
      l[[i]] = formatC(p.wlx$p[i], format = "e", digits= 2)
    } else if (p.wlx$p[i] > 0.0001 & p.wlx$p[i] <= 0.001) {
      l[[i]] = format(round(p.wlx$p[i], digits=5), nsmall = 5)
    } else {
      l[[i]] = format(round(p.wlx$p[i], digits=3), nsmall = 3)
    }
  }
  p.wlx$p_label = unlist(l)
  # p.wlx$p_label = paste0("p=", p.wlx$p_label)
  p.wlx$p_label = ifelse(p.wlx$p < .1, p.wlx$p_label, "")

  a = imm.exp %>% filter(name == "CD3_CAR_PERC", !is.na(value))
  tmp = ggplot(a, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess")
  max.v = max(ggplot_build(tmp)[[1]][[1]]$ymax)
  min.v = min(ggplot_build(tmp)[[1]][[1]]$ymin)

  pl.perc = ggplot(a, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_rect(aes(xmin=1, xmax=2, ymin=-Inf, ymax = Inf), alpha=0.025, fill="#DDDDDD", color = NA) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess") +
    coord_cartesian(ylim=c(min.v, max.v+ 5), xlim = c(1.35, 4.75)) +
    geom_text(
      data = subset(p.wlx, name == "CD3_CAR_PERC"),
      aes(x=DAY, y=Inf, label = p_signif), inherit.aes=F, size=3, vjust=1.5, hjust=.5
    ) +
    scale_color_manual(name = " ", values = palette, labels = labels) +
    scale_fill_manual(name = " ", values = palette, labels = labels) +
    ylab("% CD3+CAR+") +
    theme(
      # aspect.ratio = 1,
      legend.position = "bottom",
      axis.title.y = element_text(vjust = + 2),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      plot.title = element_text(hjust = 0.5, size = rel(1), face = "plain"),
      legend.key.spacing.x = unit(10, "pt")
    ) +
    guides(colour = guide_legend(keyheight = .7)) +
    ggtitle("CD3+CAR+")

  ###

  b = imm.exp %>% filter(name == "CAR_RATIO", !is.na(value))
  tmp = ggplot(b, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess")
  max.v = max(ggplot_build(tmp)[[1]][[1]]$ymax)
  min.v = min(ggplot_build(tmp)[[1]][[1]]$ymin)

  pl.ratio =
    ggplot(b, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_rect(aes(xmin=1, xmax=2, ymin=-Inf, ymax = Inf), alpha=0.025, fill="#DDDDDD", color = NA) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess") +
    coord_cartesian(ylim=c(min.v, max.v + .1), xlim = c(1.35, 4.75)) +
    geom_text(
      data = subset(p.wlx, name == "CAR_RATIO"),
      aes(x=DAY, y=Inf, label = p_signif), inherit.aes=F, size=3, vjust=1.5, hjust=.5
    ) +
    scale_color_manual(name = " ", values = palette, labels = labels) +
    scale_fill_manual(name = " ", values = palette, labels = labels) +
    ylab("CD4+CAR+/CD8+CAR+") +
    theme(
      # aspect.ratio = 1,
      legend.position = "bottom",
      axis.title.y = element_text(vjust = + 2),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      plot.title = element_text(hjust = 0.5, size = rel(1), face = "plain"),
      legend.key.spacing.x = unit(10, "pt")
    ) +
    ggtitle("CAR+ ratio")

  plot_grid(
    plot_grid(
      pl.perc + theme(legend.position = "none"),
      NULL,
      pl.ratio + theme(legend.position = "none"),
      ncol = 3, rel_widths = c(1, .1, 1)
    ),
    plot_grid(ggpubr::get_legend(pl.perc)),
    nrow = 2, rel_heights = c(1, .1)
  )
}

pl.exp.product = plot_expansion(
  obj.pd = pdata.clin,
  obj.imm = pdata.imm,
  group = "PRODUCT",
  palette = c(cilta = "#004488", ide = "#44AA99"),
  labels = c("cilta" = "Cilta-cel", "ide" = "Ide-cel"),
  leg.pos="bottom"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Differential abundance analysis (DA)
# Pre-filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
f = sort(rowSums(table(p.data$orig.ident, p.data$celltype) > 5))
p.data = droplevels(p.data[p.data$orig.ident %in% names(f[f > 10]), ])

ct.fltr = p.data %>%
  dplyr::group_by(celltype, TIMEPOINT, PRODUCT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(thres = n > 200) %>%
  dplyr::summarise(n = sum(thres)) %>%
  dplyr::filter(n >= 1) %>%
  dplyr::mutate(ID = paste0(celltype, TIMEPOINT)) %>%
  data.frame()
p.data = droplevels(p.data[paste0(p.data$celltype, p.data$TIMEPOINT) %in% ct.fltr$ID, ])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DA
# Base: all cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
celltype = "celltype"
grp1 = "cilta"
grp2 = "ide"
ctrst = "PRODUCT"
res.speckle = lapply(unique(as.character(p.data$TIMEPOINT)), function(tp){
  pd = droplevels(p.data[p.data$TIMEPOINT == tp, ])
  pd = droplevels(pd[pd[[ctrst]] %in% c(grp1, grp2), ])
  print(tp)
  print(table(pd[!duplicated(pd$orig.ident), ][[ctrst]]))
  pd[[ctrst]] = factor(pd[[ctrst]], levels = c(grp1, grp2))
  res = speckle::propeller(
    clusters = pd[[celltype]], sample = pd$orig.ident,
    group = pd[[ctrst]], transform = "logit"
  ) %>% dplyr::mutate(Group = tp)
}) %>% do.call("rbind", .)

# Log fold change (LFC)
res.speckle$FC = log2(
  res.speckle[[paste0("PropMean.", grp1)]] / res.speckle[[paste0("PropMean.", grp2)]]
)

# if LFC is NA or -Inf, set to max/min LFC value
if(any(is.na(res.speckle$FC) | is.infinite(res.speckle$FC))){
  print("warning: one group has 0 cell for a celltype. LFC ist set to max value. Check it")
  print(res.speckle[(is.na(res.speckle$FC) | is.infinite(res.speckle$FC)), ])
  fc = res.speckle[!is.infinite(res.speckle$FC), ]$FC
  res.speckle = res.speckle %>% dplyr::mutate(
    FC = dplyr::case_when(
      is.infinite(FC) ~ min(fc[fc < 0]),
      is.na(FC) ~ max(FC[FC > 0], na.rm = T),
      TRUE ~ FC
    )
  )
}

res.speckle$Group = factor(res.speckle$Group, levels = rev(c("LP", "Late", "Very Late")))
lvls = naturalsort(unique(as.character(res.speckle$BaselineProp.clusters)))
res.speckle$BaselineProp.clusters = factor(res.speckle$BaselineProp.clusters, levels = lvls)

f = as.character(res.speckle[res.speckle$P.Value < .1, ]$BaselineProp.clusters)
res.speckle.sign = res.speckle[res.speckle$BaselineProp.clusters %in% f, ]

# Thresholf min/max value
v.max = 2.5
res.speckle.sign$FC[res.speckle.sign$FC > v.max & !is.infinite(res.speckle.sign$FC)] = v.max
res.speckle.sign$FC[res.speckle.sign$FC < -v.max & !is.infinite(res.speckle.sign$FC)] = -v.max

da.pl = da_tile_pl(
  df = res.speckle.sign,
  grp1 = "Cilta",
  grp2 = "Ide",
  pl.title = "Comparision of cell type composition between patients\ntreated with Cilta and Ide",
  tile.size = 3,
  stroke.size = 1
) + theme(plot.title = element_text(hjust = 0.5, face = "plain"))

# sc_ct_sample_fraction(
#   inpMeta = p.data[p.data$celltype %in% unique(res.speckle.sign$BaselineProp.clusters), ],
#   label = "celltype",
#   group.facet = "TIMEPOINT",
#   group.color = "PRODUCT",
#   nbr.cell.cut = 10
# ) +
#   scale_y_continuous(
#     trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10),
#     breaks=c(0, 0.001, 0.01, 0.1, 1)
#   )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DA
# Late vs. LP, Very late vs. LP | paired desing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data.w = se.meta@meta.data
p.data.w$TIMEPOINT = gsub(" ", "", p.data.w$TIMEPOINT)

p.data.pr = split(p.data.w, p.data.w$PRODUCT)
p.data.tp = split(p.data.w, p.data.w$TIMEPOINT)

# For some reason, variables that are passed to the limma::makeContrasts
# function must be defined globally
grp1 = "Late"; grp2 = "LP"; da.late.lp = da_paired(p.data.pr, only.paired = F)
grp1 = "VeryLate"; grp2 = "Late"; da.vlate.late = da_paired(p.data.pr, only.paired = F)

da.tps = rbind(da.late.lp, da.vlate.late)
da.tps = da.tps[da.tps$Present.Grp1 > 2 | da.tps$Present.Grp2 > 2, ]

alpha = .1
f = da.tps %>% group_by(CLUSTER) %>%
  filter(P.Value < alpha) %>%
  pull(CLUSTER)
da.tps = da.tps[da.tps$CLUSTER %in% f, ]

da.tps$LogFC_DIR = ifelse(da.tps$LogFC_cont > 0, "up", "down")

# tp = p.data.tp$Late
freq.l = parallel::mclapply(p.data.tp, function(tp) {
  pd = droplevels(tp)
  prop = speckle::getTransformedProps(pd$celltype, pd$orig.ident, "asin")
  prop = as.data.frame(prop$TransformedProps)
  # prop = as.data.frame(prop$Proportions)
  prop$sample = pd$PATIENT_ID[match(prop$sample, pd$orig.ident)]
  prop
}, mc.cores = 1)

freq.diff = freq.l$LP
freq.diff$Freq_2 = freq.l$Late$Freq[
  match(paste0(freq.diff$clusters, freq.diff$sample), paste0(freq.l$Late$clusters, freq.l$Late$sample))
]
freq.diff$Freq_3 = freq.l$VeryLate$Freq[
  match(paste0(freq.diff$clusters, freq.diff$sample), paste0(freq.l$VeryLate$clusters, freq.l$VeryLate$sample))
]

freq.diff = freq.diff[rowSums(is.na(freq.diff)) != 2, ]

freq.diff$diff_2_1 = freq.diff$Freq_2 - freq.diff$Freq
freq.diff$diff_3_2 = freq.diff$Freq_3 - freq.diff$Freq_2

freq.diff = rbind(
  freq.diff %>% select(clusters, sample, diff = diff_2_1) %>% mutate(TP = "T12"),
  freq.diff %>% select(clusters, sample, diff = diff_3_2) %>% mutate(TP = "T23")
)

pd = se.meta@meta.data
pd = pd[!duplicated(pd$PATIENT_ID), ]
freq.diff$GROUP = pd$PRODUCT[match(freq.diff$sample, pd$PATIENT_ID)]

p.val.df = da.tps %>% dplyr::select(
  clusters = CLUSTER, GROUP = GROUP, group1 = CTRST, group2 = CTRST, p = FDR
)
p.val.df$p = ifelse(p.val.df$p < .1, p.val.df$p, 1)
pval = p.val.df$p
pval = lapply(pval, function(pval){
  if (pval <= 0.0001) {
    pval = formatC(pval, format = "e", digits= 2)
  } else if (pval > 0.0001 & pval <= 0.001) {
    pval = format(round(pval, digits=5), nsmall = 5)
  } else {
    pval = format(round(pval, digits=3), nsmall = 3)
  }

})
p.val.df$p = gsub("1.000", "", unlist(pval))

p.val.df = p.val.df %>% dplyr::mutate(.y. = "diff")
p.val.df$group1 = gsub("Late_vs_LP", "T12", p.val.df$group1)
p.val.df$group1 = gsub("VeryLate_vs_Late", "T23", p.val.df$group1)
p.val.df$group2 = gsub("Late_vs_LP", "T12", p.val.df$group2)
p.val.df$group2 = gsub("VeryLate_vs_Late", "T23", p.val.df$group2)
p.val.df$TP = p.val.df$group2

max.y = freq.diff %>%
  dplyr::group_by(clusters) %>%
  dplyr::summarise(max = max(diff, na.rm = T))
p.val.df$yPos = max.y$max[match(p.val.df$clusters, max.y$clusters)]

check.singl = table(paste0(p.val.df$clusters, p.val.df$GROUP))
check.singl = lapply(names(check.singl), function(i){
  if(check.singl[i] == 1) {
    r = p.val.df[paste0(p.val.df$clusters, p.val.df$GROUP) %in% i, , drop = F]
    r$p = "-"
    r$TP = ifelse(r$TP == "T12", "T23", "T12")
    return(r)
  }
})
check.singl = check.singl[lengths(check.singl) != 0] %>% do.call("rbind", .)
p.val.df = rbind(p.val.df, check.singl)

p.val.pl.all =
  ggplot(freq.diff, aes(GROUP, diff, fill = TP)) +
  geom_boxplot(outlier.size = .1, outliers = F, alpha = 1, fatten = 1.25) +
  geom_point(
    position=position_jitterdodge(jitter.width = .1), size = .01,
    color = "black", alpha = .5
  ) +
  facet_wrap(~ clusters, scales = "free_y", ncol = 6) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_fill_manual(
    values = c("#33BBEE", "#0077BB"),
    labels = c("T12" = "%Late - %LP", "T23" = "%Very Late -%Late")
  ) +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(.5, "lines")
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("% diff") + xlab(NULL) + labs(fill = "% Diff") +
  stat_pvalue_manual(
    p.val.df, x = "GROUP",
    y.position = p.val.df$yPos -(p.val.df$yPos / 100 * 15), color = "black",
    label = "p", size = 2,
    position = position_dodge(0.8)
  )

ggsave2(
  filename="publication/extended_data_files/Fig_5_supps_da_paired.png",
  p.val.pl.all,
  width = 180, height = 180, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

q = "Mono CD16$|EOMES|TEX|EM$"
freq.diff.m = freq.diff[grepl(q, freq.diff$clusters, ignore.case = T), ]
p.val.df.m = p.val.df[grepl(q, p.val.df$clusters, ignore.case = T), ]

p.val.pl.main =
  ggplot(freq.diff.m, aes(GROUP, diff, fill = TP)) +
  geom_boxplot(outlier.size = .01, outliers = T, alpha = 1, outlier.alpha = .5,
               fatten = 1.25, linewidth = .3
  ) +
  facet_wrap(clusters ~., scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0, linewidth = .3, linetype = "dashed") +
  scale_fill_manual(
    values = c("#33BBEE", "#0077BB"),
    labels = c("T12" = "%Late - %LP", "T23" = "%Very Late -%Late")
  ) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t=-6),
    panel.spacing = unit(.25, "lines"),
    legend.key.spacing.y = unit(-2, "pt"),
    strip.text = element_text(hjust = 0.5, margin = margin(2,0,4.5,0))
  ) +
  guides(fill=guide_legend(nrow=2, title = NULL)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("% Difference of cell type proportions") + xlab(NULL) + labs(fill = NULL) +
  stat_pvalue_manual(
    p.val.df.m, x = "GROUP", y.position = p.val.df.m$yPos, color = "black",
    #vjust = -.75,
    vjust = 1,
    label = "p", size = 1.5,
    position = position_dodge(0.8)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Cilta vs. Ide | All T-cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.w = se.t
se.w$celltype_short_3 = as.factor(gsub("CTL_", "CTL.", se.w$celltype))

dg.cilta.l = dgea_gex(
  obj = se.w[, se.w$PRODUCT == "cilta"], group = "TIMEPOINT", ctrs.grp1 = "Late", ctrs.grp2 = "LP",
  .target = "celltype",
  latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA"),
  split.by.tp = F, logfc.threshold = log2(1.25), min.pct = .20,
  subsample = T, threads = 15
)
dg.cilta.vl = dgea_gex(
  obj = se.w[, se.w$PRODUCT == "cilta"], group = "TIMEPOINT", ctrs.grp1 = "Very Late", ctrs.grp2 = "Late",
  .target = "celltype",
  latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA"),
  split.by.tp = F, logfc.threshold = log2(1.25), min.pct = .20,
  subsample = T, threads = 15
)

dg.ide.l = dgea_gex(
  obj = se.w[, se.w$PRODUCT == "ide"], group = "TIMEPOINT", ctrs.grp1 = "Late", ctrs.grp2 = "LP",
  .target = "celltype",
  latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA"),
  split.by.tp = F, logfc.threshold = log2(1.25), min.pct = .20,
  subsample = T, threads = 15
)
dg.ide.vl = dgea_gex(
  obj = se.w[, se.w$PRODUCT == "ide"], group = "TIMEPOINT", ctrs.grp1 = "Very Late", ctrs.grp2 = "Late",
  .target = "celltype",
  latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA"),
  split.by.tp = F, logfc.threshold = log2(1.25), min.pct = .20,
  subsample = T, threads = 15
)

dg.cilta.l$res.dgea$cluster = paste0(dg.cilta.l$res.dgea$celltype, "_", "Late.Cilta")
dg.cilta.vl$res.dgea$cluster = paste0(dg.cilta.vl$res.dgea$celltype, "_", "VeryLate.Cilta")
dg.ide.l$res.dgea$cluster = paste0(dg.ide.l$res.dgea$celltype, "_", "Late.Ide")
dg.ide.vl$res.dgea$cluster = paste0(dg.ide.vl$res.dgea$celltype, "_", "VeryLate.Ide")

dg.cilta.l$res.dgea.sign$cluster = paste0(dg.cilta.l$res.dgea.sign$celltype, "_", "Late.Cilta")
dg.cilta.vl$res.dgea.sign$cluster = paste0(dg.cilta.vl$res.dgea.sign$celltype, "_", "VeryLate.Cilta")
dg.ide.l$res.dgea.sign$cluster = paste0(dg.ide.l$res.dgea.sign$celltype, "_", "Late.Ide")
dg.ide.vl$res.dgea.sign$cluster = paste0(dg.ide.vl$res.dgea.sign$celltype, "_", "VeryLate.Ide")

dg.tp = list(
  res.dgea = rbind(dg.cilta.l$res.dgea, dg.cilta.vl$res.dgea, dg.ide.l$res.dgea, dg.ide.vl$res.dgea),
  res.dgea.sign = rbind(dg.cilta.l$res.dgea.sign, dg.cilta.vl$res.dgea.sign, dg.ide.l$res.dgea.sign, dg.ide.vl$res.dgea.sign)
)

table(dg.tp$res.dgea.sign$cluster)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Cilta vs. Ide | All T-cells | Tile plot | diff LFC
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a = dg.cilta.l$res.dgea.sign
b = dg.ide.l$res.dgea.sign
a$avg_log2FC_2 = b$avg_log2FC[match(paste0(a$celltype, a$feature), paste0(b$celltype, b$feature))]
a = a[!is.na(a$avg_log2FC_2), ]
a = a[!grepl("^RPL|^RPS|^MT-", a$feature), ]
a$avg_log2FC_DIR = sign(a$avg_log2FC)
a$avg_log2FC_2_DIR = sign(a$avg_log2FC_2)
a = a[a$avg_log2FC_DIR != a$avg_log2FC_2_DIR, ]

df = rbind(
  a %>% dplyr::select(avg_log2FC, celltype, feature) %>% mutate(PRD = "Cilta"),
  a %>% dplyr::select(avg_log2FC = avg_log2FC_2, celltype, feature)  %>% mutate(PRD = "Ide")
)

df$ID = factor(
  paste0(df$feature, "_", df$celltype),
  levels = rev(sort(unique(paste0(df$feature, "_", df$celltype))))
)
dflabs1 <- gsub("_.*", "", levels(df$ID))
dflabs2 <- gsub(".*_", "", levels(df$ID))

diff.lfc.pl =
  ggplot(df, aes(PRD, as.numeric(ID), fill = avg_log2FC)) +
  geom_tile(color = "white", lwd = .3) +
  scale_fill_scico(
    palette = "vik", midpoint = 0, begin = .1, end = .9, na.value = "white",
    breaks = breaks_pretty(3), limits = c(-max(abs(df$avg_log2FC)), max(abs(df$avg_log2FC)))
  ) +
  scale_y_continuous(
    breaks = 1:length(dflabs1), labels = dflabs1,
    sec.axis = sec_axis(~., breaks = 1:length(dflabs2), labels = dflabs2),
    expand = c(0, 0)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  guides(
    fill = guide_colorbar(
      title = "Log2FC",
      barwidth = unit(4, 'lines'), barheight = unit(.35, 'lines'),
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    )
  ) + xlab(NULL) + ylab(NULL) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05),
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
    axis.ticks = element_blank(),
    legend.ticks.length = unit(0.05, 'cm'),
    legend.position = "bottom",
    legend.margin = margin(t=-3),
    axis.text.y = element_text(size = 7)
  ) +
  ggtitle("DE genes\nLate vs. LP for all T cells")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Cilta vs. Ide | All T-cells | Barplot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

i = dg.ide.l$res.dgea.sign
i$PRODUCT = "Ide"
c = dg.cilta.l$res.dgea.sign
c$PRODUCT = "Cilta"

res.dgea.sign = rbind(i, c)
res.dgea.pl = res.dgea.sign[res.dgea.sign$feature %in% df.ftrs$feature, ]
res.dgea.pl$Type = df.ftrs$type[match(res.dgea.pl$feature, df.ftrs$feature)]
res.dgea.pl = res.dgea.pl[!grepl("gdT", res.dgea.pl$celltype), ]
res.dgea.pl$LIN = ifelse(grepl("CD4", res.dgea.pl$celltype), "CD4", "CD8")
res.dgea.pl$LIN = paste0(res.dgea.pl$LIN, " | ", res.dgea.pl$PRODUCT)
top.de =
  ggplot(res.dgea.pl, aes(avg_log2FC, feature, fill = celltype)) +
  geom_bar(stat="identity", width = .8) +
  facet_grid2(
    Type ~ LIN,
    scales = "free_y", independent = F, space = "free"
  ) +
  geom_vline(xintercept = 0, linewidth = .2, linetype = "dashed") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(margin = margin(l = 1.5, unit = "pt"), size = 6),
    panel.spacing = unit(.25, "lines"),
    legend.margin = margin(t=-2.5),
    legend.position = "bottom",
    plot.title = element_text(
      hjust = 0.5, margin = margin(3,0,0,0), face = "plain"
    )
  ) +
  labs(x = "Log2 Fold Change") +
  ggtitle("DE genes: Late vs. LP for all T cells") +
  scale_fill_manual(values = til.col) +
  guides(fill = guide_legend(title = NULL, nrow = 3))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top Clonotype composition for (clones wiith CARs)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.target = "PRODUCT"
pd = se.t@meta.data

pd = pd[!is.na(pd$CLONE_PSEUDO_ID), ]
pd$CAR_BY_EXPRS = as.logical(pd$CAR_BY_EXPRS)
pd = pd[pd$CLONE_PSEUDO_ID %in% pd[pd$CAR_BY_EXPRS == TRUE, ]$CLONE_PSEUDO_ID, ]

cl.tbl = table(pd$CLONE_PSEUDO_ID, pd$TIMEPOINT)
cl.tbl = cl.tbl[matrixStats::rowSums2(cl.tbl > 0) > 1, ]
pd = pd[pd$CLONE_PSEUDO_ID %in% rownames(cl.tbl), ]

tmp = pd %>%
  dplyr::group_by(CLONE_PSEUDO_ID) %>%
  dplyr::summarise(Frequency = sum(CAR_BY_EXPRS)) %>%
  dplyr::filter(Frequency > 0) %>%
  dplyr::pull(CLONE_PSEUDO_ID)
pd = pd[pd$CLONE_PSEUDO_ID %in% tmp, ]

set.seed(4321)
pd$barchode = rownames(pd)
pd = pd %>% group_by(orig.ident) %>% slice_sample(n=500) %>% data.frame()
rownames(pd) = pd$barchode

tmp = pd[pd$PRODUCT == "cilta", ]
tmp = droplevels(tmp[tmp$celltype_short_3 == "CD4 T-Cell", ])
table(tmp$PATIENT_ID)

pd.p40 = pd[pd$PATIENT_ID == "Patient040", ]
pd.p40 = pd.p40[pd.p40$celltype_short_3 == "CD4 T-Cell", ]

df = pd[!duplicated(paste0(pd$TIMEPOINT, pd$CLONE_PSEUDO_ID)), ] %>%
  dplyr::group_by(CLONE_PSEUDO_ID, celltype_short_3, .data[[.target]]) %>%
  dplyr::summarise(AVE = mean(clonalProportion)) %>%
  dplyr::group_by(celltype_short_3, .data[[.target]]) %>%
  dplyr::slice_max(order_by = AVE, n = 10) %>%
  data.frame()

df.pl = pd[pd$CLONE_PSEUDO_ID %in% c(unique(df$CLONE_PSEUDO_ID), unique(pd.p40$CLONE_PSEUDO_ID)), ]

df.pl =
  df.pl %>%
  dplyr::group_by(TIMEPOINT, celltype, .data[[.target]], CLONE_PSEUDO_ID) %>%
  dplyr::mutate(value = n()) %>%
  dplyr::distinct(TIMEPOINT, celltype, celltype_short_3, .data[[.target]], value) %>%
  dplyr::group_by(TIMEPOINT, .data[[.target]], CLONE_PSEUDO_ID) %>%
  dplyr::mutate(ct_sum = sum(value))
df.pl$PATIENT_ID = as.character(pd$PATIENT_ID[match(df.pl$CLONE_PSEUDO_ID, pd$CLONE_PSEUDO_ID)])
df.pl$RESPONSE = as.character(pd$BEST_RESPONSE_CONSENSUS[match(df.pl$CLONE_PSEUDO_ID, pd$CLONE_PSEUDO_ID)])

y.lables = df.pl
y.lables = y.lables[!duplicated(y.lables$CLONE_PSEUDO_ID), ]
y.lables$PATIENT_ID = gsub("Patient0", "P", y.lables$PATIENT_ID)
# y.lables$PATIENT_ID = paste0(y.lables$PATIENT_ID, "_", y.lables$RESPONSE)
y.lables = setNames(y.lables$PATIENT_ID, y.lables$CLONE_PSEUDO_ID)

top.cd4.car =
  ggplot(data = df.pl[df.pl$celltype_short_3 == "CD4 T-Cell", ], aes(x = TIMEPOINT, y = CLONE_PSEUDO_ID)) +
  PieGlyph::geom_pie_glyph(slices = 'celltype', values = 'value', aes(radius = log10(ct_sum))) +
  scale_radius(
    range = c(0.05, .2),
    breaks = c(log10(10), log10(25), log10(100), log10(400)),
    labels = c(10^log10(10), 10^log10(25), 10^log10(100), 10^log10(400))
  ) +
  # scale_radius_manual(values  = c(.5, 1 ,2.5))
  ggstats::geom_stripped_rows() +
  facet_grid(.data[[.target]] ~ celltype_short_3, scales = "free", space = "free") +
  scale_fill_manual(
    values = til.col,
    labels = setNames(gsub("CD4.|CD8.", "", df.pl$celltype), df.pl$celltype)
  ) +
  theme(
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(l=-5),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05, size = rel(1)),
    # axis.title = element_blank(),
    axis.ticks  = element_blank()
  ) +
  guides(
    fill = guide_legend(title =  NULL, order = 1),
    radius = guide_legend(title = "Nbr. of cells", order = 2)
  ) +
  scale_y_discrete(labels = y.lables) +
  xlab(NULL) + ylab("Clones")

top.cd8.car =
  ggplot(data = df.pl[df.pl$celltype_short_3 == "CD8 T-Cell", ], aes(x = TIMEPOINT, y = CLONE_PSEUDO_ID)) +
  geom_pie_glyph(slices = 'celltype', values = 'value', aes(radius = log10(ct_sum))) +
  scale_radius(
    range = c(0.05, .2),
    breaks = c(log10(10), log10(50), log10(500), log10(1500)),
    labels = c(10^log10(10), 10^log10(50), 10^log10(500), 10^log10(1500))
  ) +
  ggstats::geom_stripped_rows() +
  facet_grid(.data[[.target]] ~ celltype_short_3, scales = "free", space = "free") +
  scale_fill_manual(
    values = til.col,
    labels = setNames(gsub("CD4.|CD8.", "", df.pl$celltype), df.pl$celltype)
  ) +
  theme(
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(l=-5),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05, size = rel(1)),
    axis.title = element_blank(),
    axis.ticks  = element_blank()
  ) +
  guides(
    fill = guide_legend(title =  NULL, order = 1),
    radius = guide_legend(title = "Nbr. of cells", order = 2)
  ) +
  scale_y_discrete(labels = y.lables)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Barplot: celltype composition for singletons and expanded clones
# grouped by CAR Product
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.t@meta.data

# pd = pd[!is.na(pd$CLONE_PSEUDO_ID), ]
# pd$CAR_BY_EXPRS = as.logical(pd$CAR_BY_EXPRS)
# pd = pd[pd$CLONE_PSEUDO_ID %in% pd[pd$CAR_BY_EXPRS == TRUE, ]$CLONE_PSEUDO_ID, ]
pd = pd[pd$CAR_BY_EXPRS == TRUE, ]

pd$TMP = paste(pd$orig.ident,  pd$celltype_short_3)
f = sort(table(pd$TMP)) > 4
f = f[f]
pd = pd[pd$TMP %in% names(f), ]

set.seed(4321)
pd$barchode = rownames(pd)
pd = pd %>% group_by(orig.ident) %>% slice_sample(n=500) %>% data.frame()
rownames(pd) = pd$barchode

# t = droplevels(pd[pd$PRODUCT == "cilta", ])
# t = t[t$celltype_short_3 == "CD4 T-Cell", ]
# table(t$orig.ident, t$TIMEPOINT)

pd.ct = pd %>%
  group_by(PRODUCT, TIMEPOINT, celltype_short_3, celltype) %>%
  dplyr::summarise(nbr.cells = n()) %>%
  dplyr::mutate(freq = nbr.cells / sum(nbr.cells)) %>%
  data.frame()
pd.ct$freq = pd.ct$freq * 100
pd.ct = droplevels(pd.ct)
# pd.ct[grepl("cilta", pd.ct$PRODUCT), ]

lbl = unique(pd.ct$celltype)
lbl = setNames(gsub("CD4.|CD8.", "", lbl), lbl)
ct.comp.pl.4 =
  ggplot(pd.ct[pd.ct$celltype_short_3 == "CD4 T-Cell", ], aes(x = TIMEPOINT, y = freq, fill = celltype)) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 8, margin = margin(l=2)),
    legend.key.spacing.y = unit(2, "pt"),
    legend.position = "none",
    legend.margin = margin(t=-5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = rel(1)),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
    panel.spacing = unit(.5, "lines")
  ) +
  facet_nested(~ celltype_short_3 +  PRODUCT, scales = "free", space = "free") +
  scale_fill_manual(values = til.col, labels = lbl) +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  geom_text(
    aes(label = ifelse(freq > 10, paste0(round(freq,1)), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 2))

ct.comp.pl.8 =
  ggplot(pd.ct[pd.ct$celltype_short_3 == "CD8 T-Cell", ], aes(x = TIMEPOINT, y = freq, fill = celltype)) +
  geom_bar(stat="identity") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 8, margin = margin(l=2)),
    legend.key.spacing.y = unit(2, "pt"),
    legend.position = "none",
    legend.margin = margin(t=-5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = rel(1)),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
    panel.spacing = unit(.5, "lines")
  ) +
  facet_nested(~ celltype_short_3 +  PRODUCT, scales = "free", space = "free") +
  scale_fill_manual(values = til.col, labels = lbl) +
  ylab("% of cells") + xlab(NULL) + labs(fill = NULL) +
  geom_text(
    aes(label = ifelse(freq > 10, paste0(round(freq,1)), '')),
    position = position_stack(vjust = 0.5), size = 2
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 2))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA  | Cilta vs. Ide | CAR+ cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.w = se.t[, se.t$CAR_BY_EXPRS == "TRUE"]
se.w$TMP = paste0(se.w$orig.ident, se.w$celltype_short_3)
Idents(se.w) = "TMP"
se.w = subset(se.w, downsample = 200)
Idents(se.w) = "orig.ident"

se.w$TIMEPOINT_CAR = paste0(se.w$TIMEPOINT, "_", se.w$CAR_BY_EXPRS)
pd = se.w@meta.data
pd = pd[!duplicated(pd$orig.ident), ]
df = table(se.w$PATIENT_ID, se.w$TIMEPOINT_CAR)
df = data.frame(rbind(df))
df$RESP = pd$BEST_RESPONSE_CONSENSUS[match(rownames(df), pd$PATIENT_ID)]
df = df[order(df$RESP), ]
df$PRODUCT = pd$PRODUCT[match(rownames(df), pd$PATIENT_ID)]
df[df$Late_TRUE >= 10, ]

dg.product = dgea_gex(
  obj = se.w, group = "PRODUCT", ctrs.grp1 = "cilta", ctrs.grp2 = "ide",
  .target = "celltype_short_3",
  latent.vars = c("STUDY", "BEST_RESPONSE_CONSENSUS", "nFeature_RNA"),
  split.by.tp = T, logfc.threshold = log2(1.25), min.pct = .25,
  subsample = F, .min.cells = 50, subsample.n = 200, min.de.genes = 10
)

dg = dg.product
dg$res.dgea.sign = dg$res.dgea.sign[dg$res.dgea.sign$feature %in% c("CAR-BCMA", go.immu.ftrs), ]
table(dg$res.dgea.sign$cluster)

nbr.tops = 12
tops_up = dg$res.dgea.sign %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = nbr.tops, wt = avg_log2FC) %>%
  dplyr::arrange(-avg_log2FC)

tops_down = dg$res.dgea.sign %>%
  dplyr::filter(avg_log2FC < 0) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = -nbr.tops, wt = avg_log2FC) %>%
  dplyr::arrange(-avg_log2FC)

tops = rbind(tops_up, tops_down)
tops$ORDER = 1:nrow(tops)
tops$ORDER = factor(tops$ORDER, levels = levels(as.factor(tops$ORDER)))

dge.car.barpl =
  ggplot(tops, aes(avg_log2FC, ORDER, fill = -log10(p_val_adj))) +
  geom_bar(stat="identity", width = .8) +
  scale_y_discrete(labels = setNames(tops$feature, tops$ORDER)) +
  scale_fill_scico(palette = "navia", direction = -1, end = .9) +
  # facet_wrap(~ cluster, nrow = 1, scales = "free") +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(2, "lines"),
    axis.title.y = element_blank(),
    legend.ticks.length = unit(0.05, 'cm'),
    legend.margin = margin(t=0),
    plot.margin = unit(c(4,4,0,4), "pt"),
    plot.title = element_text(hjust = 0.5, face = "plain")
  ) +
  guides(
    fill = guide_colorbar(
      title =  "-Log10(FDR)", title.vjust = 1,
      barwidth = unit(4, 'lines'),
      barheight = unit(.35, 'lines'),
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    )
  ) +
  xlab("Log2 Fold Change") +
  ggtitle("Late | DE genes\nCilta vs. Ide for CD8+CAR+ cells")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA | Cilta vs. Ide | CAR+ cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dge = dg.product$res.dgea.sign
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
names(ftrs.l.1)

ora.go = parallel::mclapply(ftrs.l.1, function(x){
  run_nmf_ora(
    genes = x, universe = rownames(se.t),
    category = "C5", subcategory = "BP"
  )
}, mc.cores = length(ftrs.l.1))

ora.go.immu = lapply(ora.go, function(x){
  x[x$pathway %in% go.immu.pathways, ]
})

ora.cars =
  ora_barpl(
    gsea.res = ora.go.immu,
    ftrs.list = ftrs.l.1,
    nbr.tops = 14,
    font.size = 8,
    min.genes = 3,
    # max.value = 1.5,
    term.length = 40,
    sort.by.padj = F,
    barwidth = unit(.35, 'lines'),
    barheight = unit(4, 'lines'),
    bar.width = .7,
  ) +
  ggtitle("ORA: Cilta vs. Ide") +
  theme(
    axis.text.x = element_text(size = rel(1)),
    axis.text.y = element_text(size = 7, lineheight = .75),
    legend.position = "bottom",
    legend.margin = margin(t=0, r = 110),
    plot.margin = unit(c(4,4,0,4), "pt"),
    plot.title = element_text(
      hjust = 0.5, margin=margin(0,150,4,0), face = "plain")
  ) +
  guides(
    fill = guide_colorbar(
      title =  "Pathway direction", title.vjust = 1.1,
      barwidth = unit(4, 'lines'),
      barheight = unit(.35, 'lines'),
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    )
  ) +
  ggtitle("Enrichment test for DE genes\ncomparing Cilta with Ide CD8+CAR+ cells")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a = plot_grid(
  NULL,
  plot_grid(
    plot_grid(
      plot_grid(
        pl.exp.product, NULL, da.pl, nrow = 3, rel_heights = c(1, .1, 1.1),
        labels = c("", "", "B"), label_fontface = "bold", label_size = 12
      )
    ),
    NULL,
    plot_grid(NULL, p.val.pl.main, nrow = 2, rel_heights = c(-.01, 1)),
    NULL,
    plot_grid(NULL, top.de, NULL, nrow = 3, rel_heights = c(-.01, 1, 0.01)),
    NULL,
    diff.lfc.pl,
    ncol = 7, rel_widths = c(2, .2, 1, .2,  2.15, .1, 1.16),
    labels = c("A", "", "C", "", "D", "", "E"), label_fontface = "bold",
    label_size = 12, vjust = .6
  ),
  nrow = 2, rel_heights = c(.015, 1)
)

b = plot_grid(
  NULL,
  plot_grid(
    plot_grid(
      NULL,
      ggdraw() +
        draw_label(
          "Clones including CAR+ cells", fontface = 'plain', hjust = 0.8, size = 8
        ),
      plot_grid(
        NULL,
        plot_grid(
          top.cd4.car, NULL, top.cd8.car, ncol = 3, rel_widths = c(1.15, .1, 1)
        ),
        NULL,
        nrow = 3, rel_heights = c(-.01, 1, -.01)
      ),
      nrow = 3, rel_heights = c(.01, .02, 1)
    ),
    NULL,
    dge.car.barpl,
    NULL,
    ora.cars,
    ncol = 5, rel_widths = c(.55, .025, .2, .05, .38),
    labels = c("F", "", "G", "", "H"), label_fontface = "bold",
    label_size = 12, vjust = .5
  ),
  nrow = 2, rel_heights = c(.015, 1)
)

ggsave2(
  filename="publication/figures_main/fig_05.png",
  plot_grid(a, NULL, b, nrow = 3, rel_heights = c(1, .035, 1)),
  width = 165, height = 133, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

ggsave2(
  filename="publication/figures_main/fig_05.pdf",
  plot_grid(a, NULL, b, nrow = 3, rel_heights = c(1, .035, 1)),
  width = 165, height = 133, dpi = 300,
  bg = "white", units = "mm", scale = 1.6,
  device = cairo_pdf
)
