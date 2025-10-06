# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "data.table", "ggplot2",
  "scales", "ggridges", "scCustomize", "cowplot", "xtable"
)
.bioc_packages = c("dittoSeq", "scDblFinder")

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

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)


source("code/helper/functions.R")
source("code/helper/styles.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

seurat.path = paste0(manifest$ukl_b3$seurat)

out.path = paste0(seurat.path, "seurat_post_souporcell_cilta_ida.Rds")

se.ide = readRDS(paste0(seurat.path, "seurat_ori_ide.Rds"))
se.cilta = readRDS(paste0(seurat.path, "seurat_ori_cilta.Rds"))

se.ide = se.ide[, intersect(colnames(se.cilta), colnames(se.ide))]
se.cilta = se.cilta[, intersect(colnames(se.cilta), colnames(se.ide))]

se.meta = se.cilta
se.meta.ori = se.meta
se.meta = NormalizeData(se.meta)

pdata.assoc = read.csv("data/metadata_basic_run_3_pub.csv", na.strings = "")
pdata.assoc$LANE_NAME = gsub("_", "-", pdata.assoc$LANE_NAME)
pdata.assoc$GENDER = ifelse(pdata.assoc$GENDER == "m", "MALE", "FEMALE")

souporcell.path = paste0(manifest$ukl_b3$work, "souporcell/output")
souporcell.clusters = list.files(
  path = souporcell.path, pattern = "clusters.tsv", full.names = T, recursive = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^MT-")
se.meta@meta.data$log10GenesPerUMI = log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)

nFeature_low_cutoff = 250
nFeature_high_cutoff = 8000
nCount_low_cutoff = 1000
nCount_high_cutoff = 100000
mt_cutoff = 15
complx_cutoff = .8

label_cells_rm = function(obj) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
        (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
        (Perc_of_mito_genes > mt_cutoff) | (log10GenesPerUMI < complx_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  tbl = as.data.frame.matrix(table(obj$orig.ident, obj$KEEP_CELL))
  tbl$PERC = round(tbl[, 2] / rowSums(tbl)* 100, 2)
  print(tbl)
  obj
}

se.meta = label_cells_rm(se.meta)
se.meta = subset(se.meta, subset = KEEP_CELL == TRUE)
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta$KEEP_CELL = NULL

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parse Souporcell output
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sex.exprs = FetchData(se.meta, vars = c("XIST", "RPS4Y1"), slot = "data")
colnames(sex.exprs) = c("FEMALE", "MALE")
# gene_list_plot = c("XIST", "RPS4Y1")
# Stacked_VlnPlot(seurat_object = se.meta, features = gene_list_plot, x_lab_rotate = TRUE, add.noise = F)

l = list()
for (i in souporcell.clusters) {
  lane = basename(dirname(i))
  lane = gsub("multi_", "", lane)
  df = read.table(i, header = T)
  df$barcode = paste0(lane, "_", df$barcode)
  df$lane = lane
  # print(table(df$status))
  l[[i]] = df
}
souporcell = do.call("rbind", l); rownames(souporcell) = NULL
souporcell.ori = souporcell

souporcell$MALE = sex.exprs$MALE[match(souporcell$barcode, rownames(sex.exprs))]
souporcell$FEMALE = sex.exprs$FEMALE[match(souporcell$barcode, rownames(sex.exprs))]
souporcell = souporcell[souporcell$barcode %in% rownames(sex.exprs), ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign clusters to samples based on average expression (XIST, RPS4Y1)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
smpl = rowSums(table(pdata.assoc$LANE_NAME, pdata.assoc$SAMPLE_NAME))
sc.multi = souporcell[!souporcell$lane %in% names(smpl[smpl == 1]), ]
sc.multi = sc.multi[sc.multi$barcode %in% colnames(se.meta), ]

sc.multi.melt = sc.multi %>%
  dplyr::filter(status == "singlet") %>%
  dplyr::select(barcode, assignment, lane, MALE, FEMALE) %>%
  reshape2::melt(id = c("barcode", "assignment", "lane"))

sc.multi.ave = sc.multi.melt %>%
  dplyr::group_by(lane, assignment, variable) %>%
  dplyr::summarise(ave = mean(value)) %>%
  data.frame()

sc.multi.max = sc.multi.ave %>%
  dplyr::group_by(lane, assignment) %>%
  slice_max(ave, n = 1) %>%
  data.frame()

sc.multi.max$sample_name = pdata.assoc$SAMPLE_NAME[match(
  paste0(sc.multi.max$lane, "_", sc.multi.max$variable),
  paste0(pdata.assoc$LANE_NAME, "_", pdata.assoc$GENDER)
)]

sc.multi.melt$sample_name = sc.multi.max$sample_name[match(
  paste0(sc.multi.melt$lane, "_", sc.multi.melt$assignment),
  paste0(sc.multi.max$lane, "_", sc.multi.max$assignment)
)]

sc.multi$sample_name = sc.multi.melt$sample_name[match(sc.multi$barcode, sc.multi.melt$barcode)]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign clusters to samples based on number of cells expressing XIST or RPS4Y1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sc.multi.max.abs = sc.multi %>%
  dplyr::filter(status == "singlet") %>%
  dplyr::select(lane, assignment, lane, MALE, FEMALE)
sc.multi.max.abs = reshape2::melt(sc.multi.max.abs ,id = c("lane", "assignment"))
sc.multi.max.abs$value = sc.multi.max.abs$value > 0
sc.multi.max.abs =
  sc.multi.max.abs %>% dplyr::group_by(lane, assignment, variable) %>%
  dplyr::summarise(sum = sum(value)) %>%
  dplyr::group_by(lane, assignment) %>%
  dplyr::filter(sum == max(sum)) %>%
  data.frame()
sc.multi.max.abs$sample_name = pdata.assoc$SAMPLE_NAME[match(
  paste0(sc.multi.max.abs$lane, "_", sc.multi.max$variable),
  paste0(pdata.assoc$LANE_NAME, "_", pdata.assoc$GENDER)
)]

sc.multi.max.abs$sample_name_2 = sc.multi.max$sample_name[match(
  paste0(sc.multi.max.abs$lane, "_", sc.multi.max.abs$assignment),
  paste0(sc.multi.max$lane, "_", sc.multi.max$assignment)
)]

stopifnot(identical(sc.multi.max.abs$sample_name_2, sc.multi.max.abs$sample_name))

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# df = sc.multi
# df = subset(df, status == "singlet")
# df = df %>% dplyr::select(lane, cluster = assignment, RPS4Y1 = MALE, XIST = FEMALE)
# df = reshape2::melt(df ,id = c("lane", "cluster"))
# df$value = df$value > 0
# df = df %>% dplyr::group_by(lane, cluster, variable) %>%
#   dplyr::summarise(sum = sum(value))
#
# ggplot(data=df, aes(x = cluster, y=sum, fill = variable)) +
#   geom_bar(stat="identity", position=position_dodge()) +
#   mytheme(base_size = 10) +
#   theme(
#     aspect.ratio = 1,
#     panel.spacing = unit(1, "lines")
#   ) +
#   facet_wrap(~ lane, scales = "free_y") +
#   ylab("Nbr. of cells expressing RPS4Y1 or XIST (>0 counts)") +
#   xlab("Cluster assignment by Souporcell") +
#   labs(fill = NULL)

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# df = sc.multi
# df = subset(df, status == "singlet")
#
# df = df %>% dplyr::select(MALE, FEMALE, sample_name, lane, assignment)
# df = reshape2::melt(df ,id = c("sample_name", "lane", "assignment"))
# axis.cut = quantile(c(df$value), .99)
# df$value[df$value >= axis.cut] = axis.cut
# df$sample_name = paste0(df$assignment, " | ", df$sample_name)
# df$lane = gsub("Lane", "Well", df$lane)
# df = df[df$value != 0, ]
#
# df$variable = ifelse(df$variable == "MALE", "RPS4Y1", "XIST")
#
# # pl.mf.dens =
#   ggplot(df, aes(x = value, y = variable, fill = variable)) +
#   ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5, size = .2, alpha = .5) +
#   theme(
#     # aspect.ratio = 1,
#     panel.spacing = unit(1, "lines"),
#     legend.position = "none"
#   ) +
#   xlab("Gene expression") +
#   facet_wrap( ~ lane + sample_name, ncol = 6)

# ggsave2(
#   filename="analysis/souporcell/ukl_b3_dens_m_vs_r_samples.pdf",
#   pl.mf.dens,
#   width = 160, height = 200, dpi = 100, units = "mm", scale = 1
# )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sc = souporcell.ori
sc$ID = paste0(sc$lane, "_", sc$assignment)

sc$sample_name = sc.multi.max$sample_name[match(
  sc$ID,
  paste0(sc.multi.max$lane, "_", sc.multi.max$assignment)
)]

nbr.smpls.well = rowSums(table(pdata.assoc$LANE_NAME, pdata.assoc$SAMPLE_NAME))
sc$WELL_SPLIT = ifelse(sc$lane %in% names(nbr.smpls.well[nbr.smpls.well != 1]), "multi", "single")
sc$sample_name = ifelse(sc$WELL_SPLIT == "single", sc$lane, sc$sample_name)
sc$sample_name = ifelse(
  sc$WELL_SPLIT == "single",
  pdata.assoc$SAMPLE_NAME[match(sc$sample_name, pdata.assoc$LANE_NAME)],
  sc$sample_name
)
sc = subset(sc, WELL_SPLIT == "single" | sc$status == "singlet")
sc$ID = NULL
rownames(sc) = sc$barcode; sc$barcode = NULL
table(sc$sample_name, sc$assignment, useNA = "always")

se.cilta.anno = se.cilta
se.cilta.anno@meta.data$WELL = se.cilta.anno@meta.data$orig.ident
se.cilta.anno = SeuratObject::AddMetaData(se.cilta.anno, sc)
se.cilta.anno@meta.data$orig.ident = se.cilta.anno@meta.data$sample_name
se.cilta.anno@meta.data$orig.ident = factor(se.cilta.anno@meta.data$orig.ident)
se.cilta.anno@meta.data$sample_name = NULL
Idents(se.cilta.anno) = "orig.ident"

se.ide.anno = se.ide
se.ide.anno@meta.data$WELL = se.ide.anno@meta.data$orig.ident
se.ide.anno = SeuratObject::AddMetaData(se.ide.anno, sc)
se.ide.anno@meta.data$orig.ident = se.ide.anno@meta.data$sample_name
se.ide.anno@meta.data$orig.ident = factor(se.ide.anno@meta.data$orig.ident)
se.ide.anno@meta.data$sample_name = NULL
Idents(se.ide.anno) = "orig.ident"

se.ide.anno = subset(se.ide.anno, ident = subset(pdata.assoc, PRODUCT == "abecma")$SAMPLE_NAME)
se.cilta.anno = subset(se.cilta.anno, ident = subset(pdata.assoc, PRODUCT == "abecma")$SAMPLE_NAME, invert = T)

se.meta.fin = merge(se.ide.anno, se.cilta.anno)
se.meta.fin@meta.data = droplevels(se.meta.fin@meta.data)
se.meta.fin$orig.ident = as.factor(se.meta.fin$orig.ident)
DefaultAssay(se.meta.fin) = "RNA"
Idents(se.meta.fin) = "orig.ident"
se.meta.fin[["RNA"]] <- JoinLayers(se.meta.fin[["RNA"]])

table(is.na(se.meta.fin$orig.ident))
se.meta.fin = se.meta.fin[, !is.na(se.meta.fin$orig.ident)] # Cells that could not be assigned by souporcell are removed

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scDblFinder
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj.l = Split_Object(se.meta.fin, split.by = "orig.ident", threads = 40)
scDb.res = parallel::mclapply(obj.l, function(se){
  suppressWarnings({
    suppressMessages({
  set.seed(1234)
  sce = scDblFinder(GetAssayData(se, slot="counts"))
    })
  })
  df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  df$barcode = rownames(df)
  df
}, mc.cores = 50)
scDb.res = do.call("rbind", scDb.res)
rownames(scDb.res) = scDb.res$barcode
scDb.res$barcode = NULL
se.meta.fin = AddMetaData(se.meta.fin, scDb.res)
DefaultAssay(se.meta.fin) = "RNA"
Idents(se.meta.fin) = "orig.ident"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.fin
saveRDS(se.meta.fin, file = out.path)
