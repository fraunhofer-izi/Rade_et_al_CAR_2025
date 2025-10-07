# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2", "reshape2", "dplyr", "foreach", "naturalsort", "ggthemes",
  "cowplot", "clustree", "devtools", "scales", "stringr", "harmony", "MetBrewer",
  "Seurat", "future", "scCustomize", "scGate"
)
.bioc_packages = c("dittoSeq", "scds", "SingleCellExperiment", "UCell")

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
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SignatuR" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/SignatuR")
}
library(SignatuR)

if (any(!"STACAS" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/STACAS")
}
library(STACAS)

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("LOAD DATA")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

output.file = paste0(manifest$meta_pub, "integration/05_seurat_harmony_all_new.Rds")
output.file.t = paste0(manifest$meta_pub, "integration/06_seurat_harmony_t_all_new.Rds")
output.file.t.vdj = paste0(manifest$meta_pub, "integration/05_vdj_t_new.Rds")
output.file.b.vdj = paste0(manifest$meta_pub, "integration/05_vdj_b_new.Rds")

se.meta = readRDS(paste0(manifest$meta_pub, "04_seurat_anno_2_new.Rds"))

# source("code/01_cohorts/prep_pheno.R")
pdata = readRDS("data/clinical_data_for_seurat.Rds")

pd = se.meta@meta.data
pd = pd[, !colnames(pd) %in% colnames(pdata)[colnames(pdata) != "orig.ident"]]
pd$barcode = rownames(pd)

pd = merge(pd, pdata, by.x = "orig.ident", by.y = "orig.ident", all.x = T)
rownames(pd) = pd$barcode
pd$barcode = NULL
pd = pd[rownames(se.meta@meta.data), ]
stopifnot(identical(rownames(pd), rownames(se.meta@meta.data)))
se.meta@meta.data = pd

se.meta = subset(se.meta, TIMEPOINT != "ASCT")
se.meta = subset(se.meta, SOURCE != "BM")
se.meta = subset(se.meta, TIMEPOINT != "pre_Tec")
se.meta = se.meta[, !grepl("_2", se.meta$SAMPLE_ID)]
se.meta@meta.data = droplevels(se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: all celltypes")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = integration(
  obj = se.meta,
  no.ftrs = 2000,
  threads = 30,
  .nbr.dims = 20,
  run.integration = T,
  harmony.group.vars = c("PATIENT_ID")
)

DimPlot_scCustom(
  se.meta, reduction = "umap", group.by = "celltype_short_3",
  pt.size = 1, colors_use = ct.col
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# se.meta = readRDS(output.file)
se.meta.t = se.meta[, grepl("^CD4|^CD8|^gd|^dp", se.meta$celltype_short_3)]
se.meta.t@meta.data = droplevels(se.meta.t@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: CD4 -> re-anno Tregs")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.cd4 = se.meta.t[, grepl("CD4", se.meta.t$celltype)]
se.meta.cd4 = integration(
  obj = se.meta.cd4,
  no.ftrs = 1000,
  threads = 15,
  .nbr.dims = 25, run.integration = T,
  harmony.group.vars = c("PATIENT_ID")
)

DimPlot_scCustom(
  se.meta.cd4, reduction = "umap", group.by = "celltype", pt.size = .1,
  colors_use = til.col
) |
DimPlot_scCustom(
  se.meta.cd4, reduction = "umap", group.by = "RNA_snn_res.0.1", pt.size = .1
)

pd.cd4 = se.meta.cd4@meta.data
treg.cluster = pd.cd4 %>% dplyr::count(celltype, RNA_snn_res.0.1) %>%
  dplyr::group_by(RNA_snn_res.0.1) %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::filter(celltype == "CD4 Treg") %>%
  pull(RNA_snn_res.0.1) %>%
  as.character()
pd.cd4$celltype = ifelse(
  pd.cd4$RNA_snn_res.0.1 == treg.cluster,
  "CD4 Treg",
  as.character(pd.cd4$celltype)
)
treg = pd.cd4 %>% dplyr::select(celltype)
colnames(treg) = "celltype_2"
se.meta.cd4 = AddMetaData(se.meta.cd4, treg)

print(table(se.meta.cd4$celltype, se.meta.cd4$celltype_2))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: CD8 -> re-anno TEMRA")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.cd8 = se.meta.t[, grepl("CD8", se.meta.t$celltype)]
se.meta.cd8 = integration(
  obj = se.meta.cd8[, se.meta.cd8$CellCycle == "FALSE"],
  no.ftrs = 1000,
  threads = 15,
  .nbr.dims = 25, run.integration = T,
  harmony.group.vars = c("PATIENT_ID")
)

DimPlot_scCustom(
  se.meta.cd8, reduction = "umap", group.by = "celltype", pt.size = .1,
  colors_use = til.col
) |
DimPlot_scCustom(
  se.meta.cd8, reduction = "umap", group.by = "RNA_snn_res.0.1", pt.size = .1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Annotate TEMRA subcluster
# TEMRA cluster2: KLRC2 and IKZF2 high
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta.cd8) = "RNA"
exprs = FetchData(se.meta.cd8, vars = c("RNA_snn_res.0.1", "KLRC2", "IKZF2"))

tbl = table(exprs$RNA_snn_res.0.1, (exprs$KLRC2 > 1 | exprs$IKZF2 > 1))
df = as.data.frame.matrix(tbl)
temra.subcluster = rownames(df[df$`TRUE` > df$`FALSE`, ])
pd.cd8 = se.meta.cd8@meta.data
pd.cd8$celltype = ifelse(
  pd.cd8$RNA_snn_res.0.1 == temra.subcluster,
  "CD8.EMRA.2", as.character(pd.cd8$celltype)
)
pd.cd8$celltype = gsub("CD8.TEMRA", "CD8.EMRA.1", pd.cd8$celltype)
temra = pd.cd8 %>% dplyr::select(celltype)
colnames(temra) = "celltype_2"
se.meta.cd8 = AddMetaData(se.meta.cd8, temra)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tmp = se.meta
# se.meta = tmp
tmp.t = se.meta.t
# se.meta.t = tmp.t

pd.t = rbind(
  pd.cd8 %>% dplyr::select(celltype),
  pd.cd4 %>% dplyr::select(celltype)
)

se.meta.t$celltype = as.character(se.meta.t$celltype)
se.meta.t$celltype_ori = se.meta.t$celltype
se.meta.t$celltype = pd.t$celltype[match(rownames(se.meta.t@meta.data), rownames(pd.t))]
se.meta.t$celltype = ifelse(is.na(se.meta.t$celltype), se.meta.t$celltype_ori, se.meta.t$celltype)
se.meta.t$celltype = factor(se.meta.t$celltype)

se.meta$celltype_t = as.character(se.meta.t$celltype)[
  match(rownames(se.meta@meta.data), rownames(se.meta.t@meta.data))
]
se.meta$celltype = ifelse(
  is.na(se.meta$celltype_t),
  as.character(se.meta$celltype),
  as.character(se.meta$celltype_t)
)
se.meta$celltype = factor(se.meta$celltype)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: All celltypes -> Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
saveRDS(se.meta, output.file)
print("done")
# se.meta = readRDS(output.file)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Table with samples analyzed in the paper -> Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = se.meta@meta.data
pd = pd[!duplicated(pd$orig.ident), ]
saveRDS(
  pd %>% dplyr::select(PATIENT_ID, SAMPLE_ID, orig.ident),
  file = "data/single_cell_samples.Rds"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("VDJ objects for T and B-cells -> Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
vdj.df = assign_vdj(
  obj = se.meta,
  vdj = "vdj_t",
  batch = c(
    paste0(manifest$ukl_b1$data_dl, "cellranger/"),
    paste0(manifest$ukl_b2$data_dl),
    paste0(manifest$ukl_b3$data_dl, "cilta/"),
    paste0(manifest$ukl_b4$data_dl)
  ),
  present.bool = F, export.table = T
)
vdj.df$barcode = rownames(vdj.df)
vdj.df$orig.ident = se.meta$orig.ident[match(rownames(vdj.df), rownames(se.meta@meta.data))]
vdj.df = vdj.df %>% dplyr::group_by(orig.ident) %>%
  dplyr::mutate(orig.ident.nbr = sum(!is.na(CTstrict))) %>%
  data.frame()
vdj.df$clonalProportion = vdj.df$clonalFrequency / vdj.df$orig.ident.nbr

vdj.df$orig.ident.nbr = NULL
vdj.df$orig.ident = NULL
rownames(vdj.df) = vdj.df$barcode
vdj.df$barcode = NULL

saveRDS(vdj.df[!is.na(vdj.df$CTnt), ], file = output.file.t.vdj)

vdj.df = assign_vdj(
  obj = se.meta,
  vdj = "vdj_b",
  batch = c(
    paste0(manifest$ukl_b1$data_dl, "cellranger/"),
    paste0(manifest$ukl_b2$data_dl),
    paste0(manifest$ukl_b3$data_dl, "cilta/"),
    paste0(manifest$ukl_b4$data_dl)
  ),
  present.bool = F, export.table = T
)
vdj.df$barcode = rownames(vdj.df)
vdj.df$orig.ident = se.meta$orig.ident[match(rownames(vdj.df), rownames(se.meta@meta.data))]
vdj.df = vdj.df %>% dplyr::group_by(orig.ident) %>%
  dplyr::mutate(orig.ident.nbr = sum(!is.na(CTstrict))) %>%
  data.frame()
vdj.df$clonalProportion = vdj.df$clonalFrequency / vdj.df$orig.ident.nbr

vdj.df$orig.ident.nbr = NULL
vdj.df$orig.ident = NULL
rownames(vdj.df) = vdj.df$barcode
vdj.df$barcode = NULL

saveRDS(vdj.df[!is.na(vdj.df$CTnt), ], file = output.file.b.vdj)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: all T-cells")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.t = integration(
  obj = se.meta.t,
  no.ftrs = 1000,
  threads = 15,
  .nbr.dims = 20,
  run.integration = T,
  harmony.group.vars = c("PATIENT_ID")
)

# DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = "STUDY", pt.size = .1)
# DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = "celltype", pt.size = .1, colors_use = til.col) & coord_flip()
# DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = "RNA_snn_res.0.1", pt.size = .1)
# FeaturePlot_scCustom(se.meta.t, features = "nFeature_RNA")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Integration: all T-cells -> Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
saveRDS(se.meta.t, output.file.t)
# se.meta.t = readRDS(output.file.t)


