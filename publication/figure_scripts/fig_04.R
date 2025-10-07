print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Figure 4")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

Sys.setenv( TZ="Etc/GMT+1" )

Sys.setenv(RETICULATE_PYTHON = "/homes/olymp/michael.rade/.virtualenvs/r-reticulate/bin/python")
library(reticulate)
# py_config()

.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "HGNChelper", "forcats", "cowplot",
  "rlang", "remotes", "scGate", "patchwork", "openxlsx", "scCustomize", "ggpubr",
  "tidyverse", "scCustomize", "ggh4x", "scico", "ggrepel", "clinfun", "paletteer"
)

.bioc_packages = c(
  "dittoSeq", "glmGamPoi", "SummarizedExperiment", "SingleCellExperiment",
  "org.Hs.eg.db", "clusterProfiler", "speckle", "limma"
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

if (any(!"TCRanker" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('carmonalab/TCRanker')
}
library(TCRanker)

if (any(!"SeuratWrappers" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SeuratData" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

if (any(!"SeuratDisk" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

if (any(!"liana" %in% installed.packages())) {
  Sys.setenv( TZ="Etc/GMT+1" ) # without this, OmnipathR will not be installed
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('saezlab/liana')
}
library(liana)

if (any(!"scRepertoire" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github(repo = "ncborcherding/scRepertoire")
}
library(scRepertoire)

# devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
# devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
# Sys.unsetenv("GITHUB_PAT")
# devtools::install_github("jackbibby1/SCPA")
# library(SCPA)

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")
source("code/helper/adt_rna_gene_mapping.R")
source("code/helper/ora.R")
theme_set(mytheme(base_size = 8))
base.size = 8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

# T-cells
se.t = readRDS(paste0(manifest$meta_pub$work, "integration/06_seurat_harmony_t_all_new.Rds"))
vdj.t = readRDS(paste0(manifest$meta_pub$work, "integration/05_vdj_t_new.Rds"))
se.t = AddMetaData(se.t, vdj.t)

se.t$CRS_GRADE_2 = as.factor(paste0("Gr_", se.t$CRS_GRADE))
se.t@meta.data = se.t@meta.data %>% dplyr::mutate(
  CRS_GRADE_COARSE = dplyr::case_when(
    CRS_GRADE_2 == "Gr_1" | CRS_GRADE_2 == "Gr_2" ~ "Gr_1/2",
    TRUE ~ "Gr_0"
  )
)
se.t@meta.data = droplevels(se.t@meta.data)

###

se.meta = readRDS(
  paste0(manifest$meta_pub, "integration/05_seurat_harmony_all_new.Rds")
)
se.meta$CRS_GRADE_2 = as.factor(paste0("Gr_", se.meta$CRS_GRADE))
se.meta@meta.data = se.meta@meta.data %>% dplyr::mutate(
  CRS_GRADE_COARSE = dplyr::case_when(
    CRS_GRADE_2 == "Gr_1" | CRS_GRADE_2 == "Gr_2" ~ "Gr_1/2",
    TRUE ~ "Gr_0"
  )
)

se.meta@meta.data = se.meta@meta.data %>% dplyr::mutate(
  celltype_short_3 = dplyr::case_when(
    celltype_short_3 == "Progenitor" ~ celltype,
    TRUE ~ celltype_short_3
  )
)
se.meta@meta.data = droplevels(se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# max CRP vs. CRS; max CRP vs. CAR; CAR vs. CRS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata = readRDS("publication/clinicial_data/clinical_table_DF_2024_10_28.Rds")
pdata.clin = pdata$pdata.clin

# Remove samples that received a second CAR therapy
pdata.clin = pdata.clin[!grepl("_2$", pdata.clin$SAMPLE_ID), ]
# Remove batch 5 samples (not in publication)
keep.samples = as.numeric(gsub("Patient", "", pdata.clin$PATIENT_ID)) <= 65
pdata.clin = pdata.clin[keep.samples, ]

pdata.clin = pdata.clin %>% mutate(
  BEST_RESPONSE_CONSENSUS = case_when(
    BEST_RESPONSE == "CR"~ "CR",
    TRUE ~ "non-CR"
  )
)

pdata.crs = pdata$pdata.crs
crp = pdata.crs[, grepl("CRP_", colnames(pdata.crs)), ]
crp = crp %>% dplyr::mutate_if(is.character, as.numeric)
pdata.crs$CRPMAX = apply(crp[, grepl("CRP_", colnames(crp)), ], 1, max, na.rm=TRUE)
pdata.crs = pdata.crs[, !grepl("CRP_|IL6_|PCT_", colnames(pdata.crs)), ]
colnames(pdata.crs) = gsub("CRPMAX", "CRP_MAX", colnames(pdata.crs))

pdata.imm = pdata$pdata.imm.s
pdata.imm = pdata.imm[pdata.imm$DAY == "Day 7", ]

pdata.clin$CRS_GROUP = pdata.crs$CRS_GROUP[match(pdata.clin$PATIENT_ID, pdata.crs$PATIENT_ID)]
pdata.clin$TOCILIZUMAB = pdata.crs$TOCI[match(pdata.clin$PATIENT_ID, pdata.crs$PATIENT_ID)]
pdata.clin$CRS = pdata.crs$CRS[match(pdata.clin$PATIENT_ID, pdata.crs$PATIENT_ID)]
pdata.clin$CRS_GRADE = as.factor(pdata.crs$CRS_GRADE[match(pdata.clin$PATIENT_ID, pdata.crs$PATIENT_ID)])
pdata.clin$CRP_MAX = pdata.crs$CRP_MAX[match(pdata.clin$PATIENT_ID, pdata.crs$PATIENT_ID)]
pdata.clin$CRP_MAX[is.infinite(pdata.clin$CRP_MAX)] = ""
pdata.clin$CRP_MAX = as.numeric(pdata.clin$CRP_MAX)
# pdata.clin = pdata.clin[!is.na(pdata.clin$CRP_MAX) , ]

pdata.clin$CAR_POS_D7 = pdata.imm$CD3_CAR_PERC[match(pdata.clin$PATIENT_ID, pdata.imm$PATIENT_ID)]


table(pdata.clin$CRS)[2] / sum(table(pdata.clin$CRS))
table(pdata.clin$CRS_GROUP)
table(pdata.clin$CRS_GROUP)[3] / sum(table(pdata.clin$CRS_GROUP)[2:3])

table(pdata.clin[pdata.clin$PRODUCT == "ide", ]$CRS)[2] / sum(table(pdata.clin[pdata.clin$PRODUCT == "ide", ]$CRS))
table(pdata.clin[pdata.clin$PRODUCT != "ide", ]$CRS)[2] / sum(table(pdata.clin[pdata.clin$PRODUCT != "ide", ]$CRS))

chisq.test(table(
  pdata.clin[pdata.clin$PRODUCT == "ide", ]$CRS,
  pdata.clin[pdata.clin$PRODUCT == "ide", ]$BEST_RESPONSE_CONSENSUS
), simulate.p.value = TRUE)
chisq.test(table(
  pdata.clin[pdata.clin$PRODUCT != "ide", ]$CRS,
  pdata.clin[pdata.clin$PRODUCT != "ide", ]$BEST_RESPONSE_CONSENSUS
), simulate.p.value = TRUE)

chisq.test(table(pdata.clin$CRS_GROUP, pdata.clin$BEST_RESPONSE_CONSENSUS))
chisq.test(table(pdata.clin$CRS, pdata.clin$BEST_RESPONSE_CONSENSUS))
chisq.test(table(pdata.clin$CRS_GRADE, pdata.clin$PRODUCT), simulate.p.value = TRUE)

set.seed(4321)
jonck = clinfun::jonckheere.test(
  pdata.clin$CRP_MAX, as.numeric(pdata.clin$CRS_GRADE),
  alternative = c("two.sided"), nperm = 1000
)
jonck = data.frame(group1 = 0, group2 = 2, pval = jonck$p.value)
jonck$pval = paste0("p = ", round(jonck$pval, 3))

crp.vs.crs.pl =
  ggplot(pdata.clin, aes(CRS_GRADE, CRP_MAX)) +
  geom_boxplot(outlier.color = NA, linewidth = .3, fatten = 1) +
  geom_dotplot(
    aes(fill = CRS_GROUP), binaxis='y', stackdir='center', stroke = .3, dotsize = 1.5
  ) +
  theme(
    axis.title.y = element_text(vjust = + 2, size = 7),
    axis.title.x = element_text(size = 7),
    legend.box.spacing = unit(0, "pt"),
    legend.text = element_text(margin = margin(l = 1, unit = "pt")),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  ggpubr::stat_pvalue_manual(
    jonck,
    y.position = 360, vjust = .5, x = "group2",
    label = "pval", size = 2.75
  ) +
  ylab("Max. CRP post [ng/l]") +
  xlab("CRS Grade") +
  scale_fill_manual(
    values = c("#6699CC", "#EE99AA", "#994455"),
    labels = c("no" = "no", "yes+Toci" = "yes & Toci", "yes-Toci" = "yes")
  )

###

crp.vs.car.pl =
  ggplot(pdata.clin, aes(CAR_POS_D7, CRP_MAX)) +
  geom_smooth(data = pdata.clin, aes(CAR_POS_D7, CRP_MAX), method=lm, se=T, linewidth = .3, color = "black") +
  geom_point(size = 1, aes(color = CRS_GROUP)) +
  ggpubr::stat_cor(label.sep='\n', method = "spearman", size = 2.75, color = "black", label.y = 320) +
  #  scale_x_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks=c(0, 2, 10, 50)
  ) +
  ylab("Max. CRP post [ng/l]") +
  xlab("Day 7 | % of CD3+CAR+") +
  scale_color_manual(
    values = c("#6699CC", "#EE99AA", "#994455"),
    labels = c("no" = "no", "yes+Toci" = "yes & Tocilizumab", "yes-Toci" = "yes")
  ) +
  theme(
    axis.title.y = element_text(vjust = + 2, size = 7),
    axis.title.x = element_text(size = 7),
    legend.text = element_text(margin = margin(l = 1, unit = "pt")),
    legend.margin = margin(l = 3, t = -2)
  ) +
  guides(
    color = guide_legend(
      title = "CRS", override.aes = list(shape = 16, size = 3)
    )
  )

###

set.seed(4321)
jonck.car = clinfun::jonckheere.test(
  pdata.clin$CAR_POS_D7, as.numeric(pdata.clin$CRS_GRADE),
  alternative = c("two.sided"), nperm = 1000
)
jonck.car = data.frame(group1 = 0, group2 = 2, pval = jonck.car$p.value)
jonck.car$pval = paste0("p = ", round(jonck.car$pval, 3))

car.vs.crs.pl =
  ggplot(pdata.clin, aes(CRS_GRADE, CAR_POS_D7)) +
  geom_boxplot(outlier.color = NA, linewidth = .3, fatten = 1) +
  geom_dotplot(
    aes(fill = CRS_GROUP), binaxis='y', stackdir='center', stroke = .3, dotsize = 1.5
  ) +
  # scale_y_sqrt(breaks = c(0, 5, 25, 50, 100), limits = c(0, 100)) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks=c(0, 2, 10, 50), limits = c(0, 140)
  ) +
  theme(
    axis.title.y = element_text(vjust = + 2, size = 7),
    axis.title.x = element_text(size = 7),
    legend.box.spacing = unit(0, "pt"),
    legend.text = element_text(margin = margin(l = 1, unit = "pt")),
    # legend.title = element_blank(),
    legend.position = "right"
  ) +
  ggpubr::stat_pvalue_manual(
    jonck.car,
    y.position = 130, vjust = .8, x = "group2",
    label = "pval", size = 2.75
  ) +
  ylab("Day 7 | % of CD3+CAR+") +
  xlab("CRS Grade") +
  labs(fill = "CRS") +
  scale_fill_manual(
    values = c("#6699CC", "#EE99AA", "#994455"),
    labels = c("no" = "no", "yes+Toci" = "yes & Tocilizumab", "yes-Toci" = "yes")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Differential abundance analysis (DA)
# Pre-filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
f = sort(rowSums(table(p.data$orig.ident, p.data$celltype) > 5))
p.data = droplevels(p.data[p.data$orig.ident %in% names(f[f > 10]), ])

ct.fltr = p.data %>%
  dplyr::group_by(celltype, TIMEPOINT, CRS) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(thres = n > 200) %>%
  dplyr::summarise(n = sum(thres)) %>%
  dplyr::filter(n >= 1) %>%
  dplyr::mutate(ID = paste0(celltype, TIMEPOINT)) %>%
  data.frame()
p.data = droplevels(
  p.data[paste0(p.data$celltype, p.data$TIMEPOINT) %in% ct.fltr$ID, ]
)

t = p.data[!duplicated(p.data$orig.ident), ]
table(t$TIMEPOINT, t$CRS)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DA
# Base: all cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
celltype = "celltype"
grp1 = "yes"
grp2 = "no"
ctrst = "CRS"
p.data[[ctrst]] = factor(p.data[[ctrst]])

res.speckle = lapply(unique(as.character(p.data$TIMEPOINT)), function(tp){
  pd = droplevels(p.data[p.data$TIMEPOINT == tp, ])
  pd = droplevels(pd[pd[[ctrst]] %in% c(grp1, grp2), ])
  print(tp)
  print(table(pd[!duplicated(pd$orig.ident), ][[ctrst]]))
  pd[[ctrst]] = factor(pd[[ctrst]], levels = c(grp1, grp2))
  res = speckle::propeller(
    clusters = pd[[celltype]], sample = pd$orig.ident, group = pd[[ctrst]], transform = "logit"
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
res.speckle$BaselineProp.clusters = factor(
  res.speckle$BaselineProp.clusters, levels = lvls
)

res.speckle[grepl("CD4 Naive", res.speckle$BaselineProp.clusters), ]

f = as.character(res.speckle[res.speckle$P.Value < .1, ]$BaselineProp.clusters)
res.speckle.sign = res.speckle[res.speckle$BaselineProp.clusters %in% f, ]

# Thresholf min/max value
v.max = 2.5
res.speckle.sign$FC[res.speckle.sign$FC > v.max & !is.infinite(res.speckle.sign$FC)] = v.max
res.speckle.sign$FC[res.speckle.sign$FC < -v.max & !is.infinite(res.speckle.sign$FC)] = -v.max

da.pl =
  da_tile_pl(
    df = res.speckle.sign,
    grp1 = "CRS",
    grp2 = "non-CRS",
    pl.title = "Comparison of cell type composition between patients\nwith CRS Grade 1 or 2 and Grade 0",
    tile.size = 3.6, stroke.size = 1.5
  ) + theme(
    plot.title = element_text(hjust = 0.5, face = "plain"),
    legend.position = "right",
    legend.title = element_text(margin = margin(l = -1, b = 3, unit = "pt"), size = rel(1)),
    legend.margin = margin(b = -28, t = 26, l = 2),
    legend.box.margin = margin(20, 0, 0, 0),
  ) +
  guides(
    shape = guide_legend(order = 3, override.aes = list(size = 3, stroke = 1)),
    color = guide_colorbar(
      title = paste0("Log2 Fold Change\nCRS vs. non-CRS"), order = 1,
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt,
      barwidth = unit(.35, 'lines'), barheight = unit(3, 'lines')
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TCR Ranker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
customSignatures = list(
  "Cytotoxicity" = c("PRF1", "GNLY", "GZMB", "GZMH", "GZMK", "NKG7")
)

clonotyp_einrch = function(se.obj = NULL){
  se.obj = se.obj[, !is.na(se.obj$CTstrict)]
  se.obj@meta.data = droplevels(se.obj@meta.data)
  res = parallel::mclapply(as.character(unique(se.obj$TIMEPOINT)), function(x){
    obj = se.obj[, se.obj$TIMEPOINT == x]
    tcr.ranking <- TCRanker(
      query = obj, tcr = "CTstrict", species = "human", minClonSize = 5,
      filterCell = NULL, group = "CRS", signature = customSignatures
    )
    tcr.ranking$TIMEPOINT = x
    tcr.ranking$singleton = ifelse(tcr.ranking$size == 1, "yes", "no")
    tcr.ranking
  }, mc.cores = 1)
  do.call("rbind", res)
}

se.w = se.t[, se.t$celltype_short_3 == "CD4 T-Cell"]
cl.enrich.cd4 = clonotyp_einrch(se.obj = se.w)
cl.enrich.cd4$LIN = "CD4 T-Cell"
se.w = se.t[, se.t$celltype_short_3 == "CD8 T-Cell"]
cl.enrich.cd8 = clonotyp_einrch(se.obj = se.w)
cl.enrich.cd8$LIN = "CD8 T-Cell"

cl.enrich = rbind(cl.enrich.cd4, cl.enrich.cd8)

df.cl = cl.enrich %>%
  dplyr::group_by(LIN, TIMEPOINT, group) %>% slice_max(freq, n = 50)
df.cl$TIMEPOINT = factor(df.cl$TIMEPOINT, levels = c("LP", "Late", "Very Late"))

target = "Cytotoxicity.score"

cl.pl =
  ggplot(df.cl, aes(group, .data[[target]], color = .data[[target]])) +
  geom_point(
    aes(size = size, group = group),
    position = ggplot2::position_jitter(width = .15, seed = 0),
  ) +
  scale_size(range = c(.1, 2)) +
  facet_nested(~ LIN + TIMEPOINT) +
  stat_summary(
    mapping = aes(x = group, y = .data[[target]]),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, color = "white", size = .2, lwd = .5
  ) +
  stat_summary(
    mapping = aes(x = group, y = .data[[target]]),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median, color = "black", size = .1, lwd = .25
  ) +
  stat_compare_means(label.x = 1.5, hjust = .5, vjust = 1, size = 2.75, label = "p.format") +
  scale_color_scico(
    palette = "bilbao", direction = -1, alpha = .4,
    breaks = scales::pretty_breaks(4)
  ) +
  # guides(color = "none") +
  guides(
    color = guide_colorbar(
      title = "Cytotoxicity score  ", alpha = 1, title.vjust = 0.65,
      barwidth = unit(4.5, 'lines'), barheight = unit(.35, 'lines'),
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    )
  ) +
  labs(size = "Clonotype Size") +
  xlab("CRS") +
  ylab("Cytotoxicity score") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  coord_cartesian(ylim=c(.1, 1)) +
  theme(
    panel.spacing = unit(.5, "lines"),
    axis.title.y = element_text(vjust = + 2),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2),
    legend.position = "bottom",
    legend.margin = margin(l = 5, t = -2),
    legend.ticks.length = unit(0.05, 'cm'),
    legend.title = element_text(margin = margin(l = 3, r = 2, unit = "pt"), size = rel(1)),
    # plot.margin = unit(c(0,4,4,4), "pt")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Polyfunctional classification
# https://github.com/magnessa/EudraCT_2016-004043-36/blob/main/bioinformatics_analysis.ipynb
# https://doi.org/10.1158/1078-0432.CCR-23-0178
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.w = se.t
exprs.data = se.w@assays$RNA@counts

regulatory.genes <- c("IL4", "IL10", "IL13", "IL22", "TGFB1","TNFRSF9", "CD40LG")
regulatory.genes <- intersect(regulatory.genes, rownames(exprs.data))
names(regulatory.genes) = rep("regulatory", length(regulatory.genes))
effector.genes <- c("GZMB", "IFNG", "CCL3", "PRF1", "TNF","LTA")
effector.genes <- intersect(effector.genes, rownames(exprs.data))
names(effector.genes) = rep("effector", length(effector.genes))
stimulatory.genes <- c("CSF2", "IL2", "IL5", "IL7", "IL8", "IL9", "IL12A", "IL15", "IL21")
stimulatory.genes <- intersect(stimulatory.genes, rownames(exprs.data))
names(stimulatory.genes) = rep("stimulatory", length(stimulatory.genes))
chemoattractive.genes <- c("CCL11", "CXCL10", "CCL4", "CCL5")
chemoattractive.genes <- intersect(chemoattractive.genes, rownames(exprs.data))
names(chemoattractive.genes) = rep("chemoattractive", length(chemoattractive.genes))
inflammatory.genes <- c("IL1B", "IL6", "IL17A", "IL17F", "CCL2", "CCL13")
inflammatory.genes <- intersect(inflammatory.genes, rownames(exprs.data))
names(inflammatory.genes) = rep("inflammatory", length(inflammatory.genes))

all.genes <- c(
  regulatory.genes, effector.genes, stimulatory.genes, chemoattractive.genes,
  inflammatory.genes
)

background_calc <- function(gene){
  av <- mean(exprs.data[gene, exprs.data[gene,]>0])
  if(is.na(av)){
    av <- 0
  }
  return(av)
}
background <- function(genelist) (sapply(genelist, background_calc))

n.genes = parallel::mclapply(unique(names(all.genes)), function(x){
  gene.set = all.genes[names(all.genes) == x]
  df = data.frame( colSums(exprs.data[gene.set,] > background(gene.set)) )
  colnames(df) = paste0("n.", x, ".genes")
  df
}, mc.cores = 5)
se.w = AddMetaData(se.w, do.call("cbind", n.genes))

polyf.classes <- c(
  "n.regulatory.genes", "n.effector.genes", "n.stimulatory.genes",
  "n.chemoattractive.genes", "n.inflammatory.genes"
)

thr <- 1 # number of genes needed to consider that cell having that function
se.w@meta.data$polyfunctionality <- rowSums(se.w@meta.data[polyf.classes] >= thr)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Polyfunctional plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.sub = se.w
se.sub = se.sub[, grepl("CD8", se.sub$celltype)]
se.sub = se.sub[, !is.na(se.sub$CTstrict)]
se.sub@meta.data = droplevels(se.sub@meta.data)

orig.idents = unique(as.character(se.sub$orig.ident))

polyf_table_patient = lapply(orig.idents, function(x){
  this_patient <- se.sub@meta.data %>% dplyr::filter(orig.ident == x)
  polyf_profile_this_patient <- this_patient %>% dplyr::count(polyfunctionality)
  n_cells_this_patient <- nrow(this_patient)
  polyf_profile_this_patient$prop <- polyf_profile_this_patient$n / n_cells_this_patient
  polyf_profile_this_patient$orig.ident = x
  polyf_profile_this_patient
})

polyf_table_patient = do.call("rbind", polyf_table_patient)
polyf.samples = polyf_table_patient %>%
  dplyr::filter(polyfunctionality > 1) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarize(number_poly = sum(n, na.rm=TRUE))

polyf_table_patient = polyf_table_patient %>%
  reshape2::dcast(polyfunctionality ~ orig.ident, value.var = "prop")
polyf_table_patient[is.na(polyf_table_patient)] <- 0
polyf_table_patient_long <- reshape2::melt((polyf_table_patient), id.vars = "polyfunctionality")
colnames(polyf_table_patient_long) <- c("polyfunctionality", "orig.ident", "polyf_frac")

df.poly = polyf_table_patient_long %>%
  dplyr::filter(polyfunctionality != 0 & polyfunctionality != 1) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarize(polyf_frac = sum(polyf_frac, na.rm=TRUE)) %>%
  droplevels()

pheno = se.sub@meta.data[!duplicated(se.sub@meta.data$orig.ident), ]
pheno$polyfunctionality = NULL

df.poly = base::merge(df.poly, pheno, by = "orig.ident", all.x = T)
df.poly$polyf_cells = polyf.samples$number_poly[match(df.poly$orig.ident, polyf.samples$orig.ident)]

poly.prop.pl =
  ggplot(df.poly, aes(CRS, polyf_frac)) +
  geom_boxplot(outlier.color = NA, linewidth = .3, fatten = 1.25) +
  geom_dotplot(
    aes(fill = CRS_GRADE), binaxis='y', stackdir='center', stroke = .3, dotsize = 1.5
  ) +
  facet_nested(~ celltype_short_3 + TIMEPOINT) +
  theme(
    panel.spacing = unit(.5, "lines"),
    axis.title.y = element_text(vjust = + 2),
    legend.text = element_text(margin = margin(l = -1, unit = "pt")),
    legend.position = "bottom",
    legend.margin = margin(t = -2),
    plot.title = element_text(hjust = 0.5, face = "plain"),
    ggh4x.facet.nestline = element_line(colour = "black", linewidth = .2)
  ) +
  stat_compare_means(label.x = 1.5, hjust = .5, vjust = 1, size = 2.75, label = "p.format") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  ylab("Fraction of\npolyfunctional cells") +
  # facet_grid(~ TIMEPOINT) +
  xlab("CRS") + labs(fill = "CRS Grade") +
  scale_fill_manual(
    values = c("#6699CC", "#EE99AA", "#994455"),
    labels = c("no" = "no", "yes+Toci" = "yes &\nToci", "yes-Toci" = "yes")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: RNA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dg.crs = dgea_gex(
  obj = se.meta, group = "CRS", ctrs.grp1 = "yes", ctrs.grp2 = "no",
  .target = "celltype",
  latent.vars = c("STUDY", "nFeature_RNA"),
  split.by.tp = T, logfc.threshold = log2(1.1), min.pct = .2,
  subsample = T, threads = 25, min.de.genes = 10
)

dg = dg.crs
dg$res.dgea.sign = dg$res.dgea.sign[
  abs(dg$res.dgea.sign$avg_log2FC) > log2(1.25),
]
table(dg$res.dgea.sign$celltype, dg$res.dgea.sign$timepoint)

# Supp Table
supps.de.sign = dg$res.dgea.sign %>%
  dplyr::select(
    Gene_name = feature, log2FC = avg_log2FC, Pvalue = p_val, FDR = p_val_adj,
    Time_point = timepoint, Cell_identity = celltype
  ) %>%
  dplyr::arrange(Time_point, Cell_identity, FDR) %>%
  dplyr::group_by(Time_point, Cell_identity) %>%
  dplyr::slice_head(n = 100)

sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info/table_de_crs_vs_noncrs.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  supps.de.sign,
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)

ftrs = c("CAR-BCMA", bm.ftrs)
ftrs = unique(ftrs[ftrs != ""])

k = sort(table(dg$res.dgea.sign$cluster)) >= 10
k = names(k[k])
dg$res.dgea.sign = dg$res.dgea.sign[
  dg$res.dgea.sign$cluster %in% k,
]
table(dg$res.dgea.sign$cluster)


cts = as.character(unique(se.meta$celltype))
ct.b = cts[grepl("^B|Plasma|CMP", cts)]
ct.t = cts[grepl("^CD|gd", cts)]
ct.moma = cts[grepl("Mono|DC", cts)]
ct.nk = cts[grepl("NK", cts)]

# T cell
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea.pl.t.lp =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "LP", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "LP", ], nbr.tops = 5,
    subset.celltypes = ct.t,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("T cell subtypes | LP")

dgea.pl.t.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ], nbr.tops = 5,
    subset.celltypes = ct.t,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -10,
    axis.min = -3 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ) +
  ggtitle("T cell subtypes | Late")

dgea.pl.t.very.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Very Late", ], nbr.tops = 5,
    subset.celltypes = ct.t,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("T cell subtypes | Very Late")

# NK
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea.pl.nk.lp =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "LP", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "LP", ], nbr.tops = 5,
    subset.celltypes = ct.nk,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("NK cell subtypes | LP")

dgea.pl.nk.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ], nbr.tops = 5,
    subset.celltypes = ct.nk,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("NK cell subtypes | Late")

dgea.pl.nk.very.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Very Late", ], nbr.tops = 5,
    subset.celltypes = ct.nk,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("NK cell subtypes | Very Late")

# MoMA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea.pl.moma.lp =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "LP", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "LP", ], nbr.tops = 5,
    subset.celltypes = ct.moma,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("Monocyte and DC subtypes | LP")

dgea.pl.moma.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ], nbr.tops = 5,
    subset.celltypes = ct.moma,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("Monocyte and DC subtypes | Late")

dgea.pl.moma.very.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Very Late", ], nbr.tops = 5,
    subset.celltypes = ct.moma,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("Monocyte and DC subtypes | Very Late")


# B Plasma
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea.pl.b.lp =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "LP", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "LP", ], nbr.tops = 5,
    subset.celltypes = ct.b,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("B cell and Plasma cell subtypes | LP")

# dgea.pl.b.late =
#   dgea_plot(
#     dg$res.dgea[dg$res.dgea$timepoint == "Late", ],
#     dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ], nbr.tops = 5,
#     subset.celltypes = ct.b,
#     box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
#     ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
#   ) +
#   ggtitle("B cell and Plasma cell subtypes | Late")

dgea.pl.b.very.late =
  dgea_plot(
    dg$res.dgea[dg$res.dgea$timepoint == "Very Late", ],
    dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Very Late", ], nbr.tops = 5,
    subset.celltypes = ct.b,
    box.padding = .1, text.repel.size = 2, text.de.nbr.size = 2.2,
    ylim.extend.up = 0, subset.tops = ftrs, legend.margin.t = -20
  ) +
  ggtitle("B cell and Plasma cell subtypes | Very Late")

#################################

ggsave2(
  filename="publication/extended_data_files/Fig_4_crs_dgea.png",
  plot_grid(
    plot_grid(
      plot_grid(
        dgea.pl.t.lp, NULL, dgea.pl.t.late, NULL, dgea.pl.t.very.late,
        nrow = 5, rel_heights = c(1, .01, 1, .01, 1),
        labels = c("a", "", "", "", ""), label_fontface = "bold",
        label_size = 11, vjust = 1.1
      ),
      NULL,
      plot_grid(
        dgea.pl.nk.lp, NULL, dgea.pl.nk.late, NULL, dgea.pl.nk.very.late,
        nrow = 5, rel_heights = c(1, .1, 1, .1, 1),
        labels = c("b", "", "", "", ""), label_fontface = "bold",
        label_size = 11, vjust = 1.1
      ),
      ncol = 3, rel_widths = c(1.75, .1, 1)
    ),
    NULL,
    plot_grid(
      plot_grid(
        dgea.pl.moma.lp, NULL, dgea.pl.moma.late, NULL, dgea.pl.moma.very.late,
        nrow = 5, rel_heights = c(1, .01, 1, .01, 1),
        labels = c("c", "", "", "", ""), label_fontface = "bold",
        label_size = 11, vjust = 1.1
      ),
      NULL,
      plot_grid(
        dgea.pl.b.lp, NULL, NULL, NULL, dgea.pl.b.very.late,
        nrow = 5, rel_heights = c(1, .1, 1, .1, 1),
        labels = c("d", "", "", "", ""), label_fontface = "bold",
        label_size = 11, vjust = 1.1
      ),
      ncol = 3, rel_widths = c(1.75, .1, 1)
    ),
    nrow = 3, rel_heights = c(1, .05, 1)
  ),
  width = 190, height = 240, dpi = 500, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)

df.t = dg$res.dgea[dg$res.dgea$timepoint == "Late", ]
df.t = df.t[grepl("CD4|CD8", df.t$celltype), ]
df.sign.t = dg$res.dgea.sign[dg$res.dgea.sign$timepoint == "Late", ]
df.sign.t = df.sign.t[grepl("CD4|CD8", df.sign.t$celltype), ]
max(df.sign.t$avg_log2FC)
min(df.sign.t$avg_log2FC)
dgea.pl.l.t =
  dgea_plot(
    dgea.res = df.t,
    dgea.res.sign = df.sign.t, nbr.tops  = 5,
    box.padding = .25,
    text.repel.size = 1.9, text.de.nbr.size = 2,
    cluster = "celltype",
    axis.max = 10, axis.min = -2.7,
    legend.margin.t = -15, cluster.label.y = .2, subset.tops = ftrs
  )  +
  ggtitle("Late | DE genes comparing CRS Grade 1 or 2 with Grade 0\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dg = dg.crs$res.dgea.sign
# dg = dg[dg$timepoint == "Late", ]
dg = dg[abs(dg$avg_log2FC) > log2(1.25), ]
dg$cluster = paste0(dg$celltype, "_", dg$timepoint)
dg$cluster = gsub("CD4 ", "CD4\n", dg$cluster)
dg$cluster = gsub("CD8 ", "CD8\n", dg$cluster)
dg$cluster = gsub("CD4\nCTL.", "CD4 CTL\n", dg$cluster)

ftrs.l = dg
ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
ftrs.l = ftrs.l[grepl("CD4|CD8", names(ftrs.l))]

ora.go.cd8 = parallel::mclapply(ftrs.l[grepl("CD8", names(ftrs.l))], function(x){
  run_nmf_ora(
    genes = x, universe = rownames(se.meta),
    category = "CD8"
  )
}, mc.cores = 1)

ora.go.cd4 = parallel::mclapply(ftrs.l[grepl("CD4", names(ftrs.l))], function(x){
  run_nmf_ora(
    genes = x, universe = rownames(se.meta),
    category = "CD4"
  )
}, mc.cores = 1)

ora.go.cd4.8 = lapply(c(ora.go.cd4, ora.go.cd8), function(x){
  x = x[x$overlap > 2, ]
  if(nrow(x) > 1){return(x)}
})
ora.go.cd4.8 = ora.go.cd4.8[lengths(ora.go.cd4.8) > 0]

ora.go.t.pl =
  ora_bubble(
    gsea.res = ora.go.cd4.8,
    ftrs.list = ftrs.l,
    nbr.tops = 15,
    min.genes = 3,
    dot.range = c(1, 4),
    # max.value = 1,
    term.length = 50,
    sort.by.padj = F,
    facet.spit = T
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.margin = margin(t=-10, b = 10),
    strip.text.x = element_text(size = rel(.7)),
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1))
  ) +
  guides(
    fill = guide_colorbar(
      title =  "Pathway\ndirection", title.vjust = 1,
      barwidth = unit(.35, 'lines'),
      barheight = unit(3, 'lines'),
      order = 1,
      ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  ) +
  ggtitle(
    "Enrichment test for DE genes comparing CRS Grade 1 or 2 with Grade 0"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GSEA: non-CRS vs. CRS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs.l = dg.crs$res.dgea

ftrs.l = split(ftrs.l, ftrs.l$cluster)
ftrs.l = lapply(ftrs.l, function(x){
  ftrs = x$feature
  names(ftrs) = x$avg_log2FC
  ftrs
})
lengths(ftrs.l)

dg.gsea = parallel::mclapply(ftrs.l, function(x){
  run_nmf_ora(genes = x, universe = rownames(se.meta), category = "H", ora = F)
}, mc.cores = 25)

gsea.df = data.table::rbindlist(dg.gsea, idcol = "CELLTYPE")
gsea.df$TIMEPOINT = gsub(".*_", "", gsea.df$CELLTYPE)
gsea.df$CELLTYPE = gsub("_.*", "", gsea.df$CELLTYPE)

alpha = .1

fltr = gsea.df %>%
  dplyr::filter(padj < alpha) %>%
  dplyr::group_by(CELLTYPE, pathway) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(pathway)

gsea.df = droplevels(gsea.df[gsea.df$pathway %in% fltr, ])

fltr = gsea.df %>%
  dplyr::filter(padj < alpha) %>%
  dplyr::group_by(CELLTYPE, TIMEPOINT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1)

gsea.df = droplevels(gsea.df[paste0(gsea.df$CELLTYPE, gsea.df$TIMEPOINT) %in% paste0(fltr$CELLTYPE, fltr$TIMEPOINT), ])

lvls = levels(factor(gsea.df$TIMEPOINT))
gsea.df$TIMEPOINT = factor(
  gsea.df$TIMEPOINT,
  levels = c(lvls[grepl("LP", lvls)], lvls[grepl("^Late", lvls)], lvls[grepl("Very Late", lvls)])
)

gsea.df$pathway = factor(gsea.df$pathway, levels = rev(levels(factor(gsea.df$pathway))))

thres = quantile(gsea.df$NES, .999, na.rm = T)
gsea.df$NES[gsea.df$NES > thres & !is.infinite(gsea.df$NES)] = thres
gsea.df$NES[gsea.df$NES < -thres & !is.infinite(gsea.df$NES)] = -thres
max.value = thres

gsea.df$LIN = ifelse(grepl("CD4|CD8|NK|Plasma", gsea.df$CELLTYPE), "Lymphoid", "Myeloid")
gsea.sub = gsea.df[!grepl("CD4|CD8|NK|Plasma", gsea.df$CELLTYPE), ]
gsea.hm =
  ggplot() +
  geom_tile(data = subset(gsea.df, padj < alpha), aes(TIMEPOINT, pathway, fill = NES), color = "white", lwd = 1.5) +
  geom_tile(data = subset(gsea.df, padj >= alpha), aes(TIMEPOINT, pathway), fill = "#EBEBEB", color = "white", lwd = 1.5) +
  # geom_point(data = subset(gsea.df, padj < .05), aes(x = TIMEPOINT, y = pathway), size = 5.5, shape = "*") +
  scale_color_manual(name="FDR<0.05", labels=NULL, values="#EBEBEB") +
  scale_fill_scico(
    palette = "vik", midpoint = 0, na.value = "#BBBBBB",
    begin = .1, end = .9, limits = c(-max.value, max.value)
  ) +
  facet_grid2( ~ CELLTYPE, space = "free", scales = "free_x") +
  guides(
    fill = guide_colorbar(
      title = "NES", order = 1,
      title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(5, 'lines'),
      ticks.linewidth = 1.5/.pt
    )
  ) + xlab(NULL) + ylab(NULL) +
  theme(
    aspect.ratio = 1,
    strip.text.x = element_text(angle = 90, hjust = 0),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05),
    plot.title = element_text(hjust = 0.5, face = "plain", size = rel(1)),
    axis.ticks.y = element_blank(),
    legend.key.size = unit(9,"pt")
  )

ggsave2(
  filename="publication/extended_data_files/Fig_4_crs_gsea.png",
  gsea.hm,
  width = 180, height = 80, dpi = 500, bg = "white", units = "mm", scale = 1.6,
  device = png, type = "cairo"
)








