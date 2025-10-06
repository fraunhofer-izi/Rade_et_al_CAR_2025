# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "cowplot", "data.table",
  "ggplot2", "ggthemes", "scGate", "patchwork", "Signac", "devtools", "scico"
)
.bioc_packages = c("UCell")

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

if (any(!"ProjecTILs" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("carmonalab/STACAS")
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

if (any(!"Azimuth" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github('satijalab/azimuth', ref = 'master')
}

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (Sys.info()["nodename"] == "ribnode020") {
  ncores = 40
} else {
  ncores = 20
}

manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta_pub, "03_seurat_anno_1.Rds"))
se.meta[["RNA"]] <- as(se.meta[["RNA"]], Class = "Assay")

output.file = paste0(manifest$meta_pub, "04_seurat_anno_2.Rds")

scGate_models_DB = get_scGateDB("data/scGateDB", force_update = F)
Idents(se.meta) = "orig.ident"

# source("code/01_cohorts/prep_pheno.R")
pdata = readRDS("data/clinical_data_for_seurat.Rds")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Coarse celltype consistency with scGATE models")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate(
    CT_L1_COARSE = dplyr::case_when(
      grepl("CD4|CD8|other T|T cell", CT_L1) ~ "T-Cell",
      grepl("Mono|DC", CT_L1) ~ "MoDC",
      TRUE ~ CT_L1
    )
  )

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    CT_L2_COARSE  = case_when(
      grepl("Eryth|Platelet|Prog_RBC", CT_L2) ~ "Erythrocyte",
      grepl("HSPC|CLP|EMP|GMP|LMPP|HSC|Prog_Mk|Prog_DC|^Prog_B|Prog Mk|pro B|pre B|BaEoMa", CT_L2) ~ "Progenitor",
      grepl("ILC|Stromal", CT_L2) ~ "Other",
      grepl("B cell|B-Cell Memory|B-Cell Naive|transitional B|B intermediate", CT_L2) ~ "B-Cell",
      grepl("Plasma|Plasmablast", CT_L2) ~ "Plasma cell",
      grepl("^cDC", CT_L2) ~ "MoDCMac",
      grepl("^pDC", CT_L2) ~ "MoDCMac",
      grepl("^ASDC|^mDC|pre-pDC|pre-mDC", CT_L2) ~ "MoDCMac",
      grepl("CD14", CT_L2) ~ "MoDCMac",
      grepl("CD16", CT_L2) ~ "MoDCMac",
      grepl("Macrophage", CT_L2) ~ "MoDCMac",
      grepl("^NK|CD56 bright NK", CT_L2) ~ "NK",
      grepl("^CD4", CT_L2) ~ "T-Cell",
      grepl("^CD8", CT_L2) ~ "T-Cell",
      grepl("gdT|dpT|dnT|T Proliferating|Treg|MAIT", CT_L2) ~ "T-Cell",
      TRUE ~ CT_L2
    )
  )

se.meta.w = se.meta

# table(se.meta.w$CT_L2, se.meta.w$CT_L2_COARSE)

se.meta.w = assign_vdj(
  obj = se.meta.w,
  vdj = "vdj_t",
  batch = c(
    paste0(manifest$ukl_b1$data_dl, "cellranger/"),
    paste0(manifest$ukl_b2$data_dl),
    paste0(manifest$ukl_b3$data_dl, "cilta/"),
    paste0(manifest$ukl_b4$data_dl)
  ),
  present.bool = T
)

se.meta.w = assign_vdj(
  obj = se.meta.w,
  vdj = "vdj_b",
  batch = c(
    paste0(manifest$ukl_b1$data_dl, "cellranger/"),
    paste0(manifest$ukl_b2$data_dl),
    paste0(manifest$ukl_b3$data_dl, "cilta/"),
    paste0(manifest$ukl_b4$data_dl)
  ),
  present.bool = T
)

hbb.exprs = FetchData(se.meta.w, vars = c("HBB"), layer = "counts")
se.meta.w = AddMetaData(se.meta.w, hbb.exprs)

pd = se.meta.w@meta.data
pd = pd[pd$SCGATE_IMMUNE == "Pure", ]; print(nrow(pd))
pd = pd[pd$DOUBLETS_CONSENSUS == FALSE | pd$CAR_BY_EXPRS == "TRUE", ]; print(nrow(pd))
pd = pd[pd$SCGATE_PLATELET == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_ERYTHROCYTE == "Impure", ]; print(nrow(pd))
pd = pd[pd$SCGATE_MEGAKARYOCYTE == "Impure", ]; print(nrow(pd))
pd = pd[!pd$CT_L2_COARSE == "Other", ]; print(nrow(pd))
pd = pd[!pd$CT_L2_COARSE == "Erythrocyte", ]; print(nrow(pd))
pd = pd[pd$HBB < 2, ]; print(nrow(pd))
pd = pd[pd$PLATELET_UCell < .35, ]; print(nrow(pd))

pd = pd[pd$CT_L2_SCORE > .5 | pd$CAR_BY_EXPRS == "TRUE", ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE != "NK" & pd$SCGATE_NK == "Pure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE == "NK" & pd$SCGATE_NK == "Impure"), ]; print(nrow(pd))
pd = pd[!(!grepl("MoDCMac|Progenitor", pd$CT_L2_COARSE) & pd$SCGATE_MYELOID == "Pure"), ]; print(nrow(pd))
pd = pd[!(grepl("Mono|ASDC|cDC|Macro|pre-mDC|pre-pDC", pd$CT_L2) & pd$SCGATE_MYELOID == "Impure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE != "T-Cell" & pd$SCGATE_TCELL == "Pure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE == "T-Cell" & pd$SCGATE_TCELL == "Impure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE != "T-Cell" & pd$CAR_BY_EXPRS == "TRUE"), ]; print(nrow(pd))
pd = pd[pd$CD4CD8_BY_EXPRS != "CD4+CD8+", ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE != "T-Cell" & pd$VDJ_T_AVAIL == TRUE), ]; print(nrow(pd))
pd = pd[!(!grepl("B-Cell|Plasma|Progenitor", pd$CT_L2_COARSE) & pd$VDJ_B_AVAIL == TRUE), ]; print(nrow(pd))
pd = pd[!(!grepl("B-Cell|Progenitor", pd$CT_L2_COARSE) & pd$SCGATE_BCELL == "Pure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE == "B-Cell" & pd$SCGATE_BCELL == "Impure"), ]; print(nrow(pd))
pd = pd[!(pd$CT_L2_COARSE != "Plasma cell" & pd$SCGATE_PLASMACELL == "Pure"), ]; print(nrow(pd))

# table(pd$CT_L2_COARSE, pd$SCGATE_MYELOID)
# table(pd$CT_L2_COARSE, pd$SCGATE_NK)
# table(pd$CT_L2_COARSE, pd$SCGATE_TCELL)
# table(pd$CT_L2_COARSE, pd$SCGATE_BCELL)
# table(pd$CT_L2_COARSE, pd$SCGATE_PLASMACELL)
# table(pd$CT_L2_COARSE, pd$CAR_BY_EXPRS)
# table(pd$CT_L2_COARSE, pd$VDJ_T_AVAIL)
# table(pd$CT_L2_COARSE, pd$VDJ_B_AVAIL)

se.meta = se.meta[ , rownames(se.meta@meta.data) %in% rownames(pd)]
se.meta@meta.data = droplevels(se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("T cells of a clone should consist of CD4 and CD8 cells.")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.t = subset(se.meta, CT_L1_COARSE == "T-Cell")
se.t = assign_vdj(
  obj = se.t,
  vdj = "vdj_t",
  batch = c(
    paste0(manifest$ukl_b1$data_dl, "cellranger/"),
    paste0(manifest$ukl_b2$data_dl),
    paste0(manifest$ukl_b3$data_dl, "cilta/"),
    paste0(manifest$ukl_b4$data_dl)
  ),
  present.bool = F
)
print(paste0("Nbr of T-Cells: ", ncol(se.t)))

# Remove inconsistency in clones between CD4 and CD8 cells
se.obj = se.t
se.obj = se.obj[, se.obj$CD4CD8_BY_EXPRS == "CD4+CD8-" | se.obj$CD4CD8_BY_EXPRS == "CD4-CD8+"]
se.obj@meta.data = droplevels(se.obj@meta.data)
se.obj$CD4CD8 = ifelse(se.obj$CD4CD8_BY_EXPRS == "CD4-CD8+", "CD8", "CD4")
rm.clones = clean_clonotypes(obj = se.obj, celltype = "CD4CD8")
length(rm.clones)
se.t = se.t[, !rownames(se.t@meta.data) %in% rm.clones]

# A clonotype must not be present in more than one patient.
se.t$PATIENT_ID = pdata$PATIENT_ID[match(se.t$orig.ident, pdata$orig.ident)]
rm.clones = clean_clonotypes_inter(obj = se.t)
length(rm.clones)
se.t = se.t[, !rownames(se.t@meta.data) %in% rm.clones]
se.t$PATIENT_ID = NULL

# Assignment of CD4- and CD8-negative cells to CD4 or CD8 clones
pd = se.t@meta.data
pd = pd[!is.na(pd$CTstrict), ]
pd = pd %>%
  dplyr::group_by(CTstrict) %>%
  dplyr::count(name = "cloneFreqAll") %>%
  dplyr::arrange(desc(cloneFreqAll))
pd$pseudo_id = paste0("Clone_", seq(1:nrow(pd)))
se.t$clonePseudoID = pd$pseudo_id[match(se.t$CTstrict, pd$CTstrict)]

df = as.data.frame.matrix(table(se.t$clonePseudoID, se.t$CD4CD8_BY_EXPRS))
stopifnot(!any(df$`CD4-CD8+` > 0 & df$`CD4+CD8-` > 0 ))
df = df %>%
  dplyr::mutate(LIN = dplyr::case_when(
    `CD4-CD8+` > 0 ~ "CD8",
    `CD4+CD8-` > 0 ~ "CD4"
  ))

lin.df = data.frame(
  row.names = rownames(se.t@meta.data),
  CD4CD8_BY_EXPRS = se.t@meta.data$CD4CD8_BY_EXPRS,
  clonePseudoID = se.t@meta.data$clonePseudoID
)
lin.df$T_LIN = df$LIN[match(lin.df$clonePseudoID, rownames(df))]
lin.df = lin.df %>% dplyr::mutate(
  T_LIN = dplyr::case_when(
    CD4CD8_BY_EXPRS == "CD4-CD8+" ~ "CD8",
    CD4CD8_BY_EXPRS == "CD4+CD8-" ~ "CD4",
    TRUE ~ T_LIN
  )
)
lin.df$T_LIN[is.na(lin.df$T_LIN)] = "T"
# table(lin.df$T_LIN, useNA = "always")
# table(se.t$CD4CD8_BY_EXPRS, useNA = "always")
se.t = AddMetaData(se.t, lin.df %>% dplyr::select(T_LIN))

# Remove CAR+ in LP samples
pd = se.t@meta.data
pd$TIMEPOINT = pdata$TIMEPOINT[match(pd$orig.ident, pdata$orig.ident)]
rm.car.fp = rownames(pd[pd$TIMEPOINT == "LP" & pd$CAR_BY_EXPRS == "TRUE", ])
se.t = se.t[, !colnames(se.t) %in% rm.car.fp]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("CD4/CD8 Imputation")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(paste0("Nbr of T-Cells: ", ncol(se.t)))
cell.idents = sort(table(se.t$orig.ident))
# print(cell.idents)
se.t = se.t[, se.t$orig.ident %in% names(cell.idents[cell.idents > 100])]
se.t@meta.data = droplevels(se.t@meta.data)

exprs.smoothed.l = parallel::mclapply(as.character(unique(se.t$orig.ident)), function(x){

  se = subset(se.t, orig.ident == x)
  # print(as.character(unique(se$orig.ident))[1])
  se = se %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F)

  exprs.smoothed = UCell::SmoothKNN(
    obj = se,
    signature.names = c("CD4", "CD8B"),
    assay="RNA", reduction="pca",
    k = 10, suffix = "_CD4CD8_smooth"
  )
  DefaultAssay(exprs.smoothed) = "RNA_CD4CD8_smooth"
  exprs.smoothed = FetchData(exprs.smoothed, vars = c("CD4", "CD8B"), layer = "data")
  colnames(exprs.smoothed) = c("CD4", "CD8")
  exprs.smoothed$GRP = se$CD4CD8_BY_EXPRS[match(rownames(exprs.smoothed), rownames(se@meta.data))]
  exprs.smoothed$SAMPLE = as.character(se$orig.ident[1])
  # noise <- rnorm(n = nrow(exprs.smoothed)) / 1e10
  # exprs.smoothed$CD4 = exprs.smoothed$CD4 + noise
  # exprs.smoothed$CD8 = exprs.smoothed$CD8 + noise
  exprs.smoothed

}, mc.cores = 40)

exprs.smoothed = do.call("rbind", exprs.smoothed.l)
exprs.smoothed$CT_L1 = se.t$CT_L1[match(rownames(exprs.smoothed), colnames(se.t))]
exprs.smoothed$T_LIN = se.t$T_LIN[match(rownames(exprs.smoothed), colnames(se.t))]
DefaultAssay(se.t) = "ADT"
adt.exprs = FetchData(se.t, vars = c("CD4", "CD8A"), layer = "data")
DefaultAssay(se.t) = "RNA"
exprs.smoothed$ADT_CD4 = adt.exprs$CD4[match(rownames(exprs.smoothed), rownames(adt.exprs))]
# exprs.smoothed$ADT_CD4[exprs.smoothed$ADT_CD4 > quantile(exprs.smoothed$ADT_CD4, .99)] = quantile(exprs.smoothed$ADT_CD4, .99)
# exprs.smoothed$ADT_CD4[exprs.smoothed$ADT_CD4 < 1] = NA
exprs.smoothed$ADT_CD8 = adt.exprs$CD8A[match(rownames(exprs.smoothed), rownames(adt.exprs))]
# exprs.smoothed$ADT_CD8[exprs.smoothed$ADT_CD8 > quantile(exprs.smoothed$ADT_CD8, .99)] = quantile(exprs.smoothed$ADT_CD8, .99)
# exprs.smoothed$ADT_CD8[exprs.smoothed$ADT_CD8 < 1] = NA

thres.x = .2
thres.y = .3

# ggplot(exprs.smoothed[1:100000, ], aes(CD4, CD8)) +
#   geom_point(aes(color = ADT_CD4), size = .5, alpha = .7) +
#   guides(alpha = 'none') +
#   scico::scale_color_scico(palette = "acton", direction = -1) +
#   # guides(colour = guide_legend(ncol = 1, override.aes = list(size=6, shape = 16, alpha = 1))) +
#   # scale_color_manual(values = ggthemes::tableau_color_pal()(10)) +
#   mytheme() +
#   theme(
#     aspect.ratio = 1,
#     panel.spacing = unit(1, "lines"),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#   ) +
#   geom_hline(yintercept = thres.y, lwd = .2, linetype = "dashed") +
#   geom_vline(xintercept = thres.x, lwd = .2, linetype = "dashed") +
#   facet_wrap( ~ T_LIN)

# Annotate T cell as CD4/8 cell based on knn smoothing
df = exprs.smoothed %>% mutate(
  T_LIN_WORK = case_when(
    (CD4 > thres.x & CD8 < thres.y) & T_LIN == "T" ~ "CD4",
    (CD4 < thres.x & CD8 > thres.y)  & T_LIN == "T" ~ "CD8",
    TRUE ~ T_LIN
  )
)
table(df$T_LIN_WORK, df$T_LIN)
sum(table(df$T_LIN_WORK, df$T_LIN)[1:2,3])  / sum(table(df$T_LIN_WORK, df$T_LIN)[, 3])

# Annotate T cell as CD4/8 cell based on ADT
df = df %>% mutate(
  T_LIN_WORK = case_when(
    (ADT_CD4 > 1 & ADT_CD8 < 1) & T_LIN_WORK == "T" ~ "CD4",
    (ADT_CD4 < 1 & ADT_CD8 > 1)  & T_LIN_WORK == "T" ~ "CD8",
    TRUE ~ T_LIN_WORK
  )
)
table(df$T_LIN_WORK, df$T_LIN)

df = df %>% mutate(
  T_LIN_WORK = case_when(
    (ADT_CD4 > 1) & T_LIN_WORK == "CD8" ~ "T",
    (ADT_CD8 > 1)  & T_LIN_WORK == "CD4" ~ "T",
    TRUE ~ T_LIN_WORK
  )
)
table(df$T_LIN_WORK, df$T_LIN)
sum(table(df$T_LIN_WORK, df$T_LIN)[1:2,3])  / sum(table(df$T_LIN_WORK, df$T_LIN)[, 3])

add.meta = df %>% dplyr::select(
  T_LIN_WORK, CD4_SMOOTH = CD4, CD8_SMOOTH = CD8
)
se.t = AddMetaData(se.t, add.meta)

# Annotate gamma delta T cells
pd = se.t@meta.data
pd = pd[pd$SCGATE_GDT == "Pure", ]
pd = pd[is.na(pd$CTstrict), ]
pd = pd[pd$T_LIN_WORK != "CD4", ]
se.t@meta.data$T_LIN_WORK = ifelse(
  rownames(se.t@meta.data) %in% rownames(pd),
  "gdT", se.t@meta.data$T_LIN_WORK
)
se.t = subset(se.t, T_LIN_WORK != "T")

# tbl = FetchData(
#   se.t,
#   vars = c("CD4_SMOOTH", "CD8_SMOOTH", "T_LIN_WORK", "T_LIN", "orig.ident", "CD4CD8_BY_EXPRS")
# )
# ggplot(tbl, aes(CD4_SMOOTH, CD8_SMOOTH)) +
#   geom_point(aes(color = T_LIN), shape = ".", alpha = .5) +
#   guides(alpha = 'none') +
#   guides(colour = guide_legend(ncol = 1, override.aes = list(size=6, shape = 16, alpha = 1))) +
#   scale_color_manual(values = ggthemes::tableau_color_pal()(10)) +
#   mytheme() +
#   theme(
#     aspect.ratio = 1,
#     panel.spacing = unit(1, "lines"),
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#   ) +
#   geom_hline(yintercept = thres.y, lwd = .2, linetype = "dashed") +
#   geom_vline(xintercept = thres.x, lwd = .2, linetype = "dashed") +
#   facet_grid( ~ T_LIN_WORK)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("ProjecTILs (T-cells)")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(paste0("Nbr of T-Cells: ", ncol(se.t)))
cd8.path = paste0(manifest$base$workdata, "references/atlases/CD8T_human_ref_v1.rds")
cd4.path = paste0(manifest$base$workdata, "references/atlases/CD4T_human_ref_v1.rds")
if(!file.exists(cd8.path)){
  options(timeout = max(900, getOption("timeout")))
  download.file("https://figshare.com/ndownloader/files/41414556", destfile = cd8.path)
}
if(!file.exists(cd4.path)){
  options(timeout = max(900, getOption("timeout")))
  download.file("https://figshare.com/ndownloader/files/39012395", destfile = cd4.path)
}
ref.cd8 <- load.reference.map(cd8.path)
ref.cd4 <- load.reference.map(cd4.path)
# DimPlot(ref.cd8, cols = til.col, label = T) + theme(aspect.ratio = 1) + ggtitle("CD8 T reference") |
# DimPlot(ref.cd4, cols = til.col, label = T) + theme(aspect.ratio = 1) + ggtitle("CD4 T reference")

# Classify CD8 T subtypes
se.cd8 <- ProjecTILs.classifier(
  subset(se.t, T_LIN_WORK == "CD8"), ref.cd8, ncores = ncores,
  split.by = "orig.ident", filter.cells = FALSE, min.confidence = 0
)

# Classify CD4 T subtypes
se.cd4 <- ProjecTILs.classifier(
  subset(se.t, T_LIN_WORK == "CD4"), ref.cd4, ncores = ncores,
  split.by = "orig.ident", filter.cells = FALSE, min.confidence = 0
)

se.tmp = merge(se.cd4, se.cd8)
# if min.confidence < .2 (ProjecTILs.classifier)
se.tmp = se.tmp[, !is.na(se.tmp$functional.cluster)]
se.t = merge(se.tmp, se.t[, se.t$T_LIN_WORK == "gdT"])
colnames(se.t@meta.data)[colnames(se.t@meta.data) == "functional.cluster"] = "SPICA_TCELL"
colnames(se.t@meta.data)[colnames(se.t@meta.data) == "functional.cluster.conf"] = "SPICA_TCELL_CONF"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge T object with se.meta
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(paste0("Nbr of T-Cells: ", ncol(se.t)))
add.meta = se.t@meta.data %>% dplyr::select(
  T_LIN = T_LIN_WORK, SPICA_TCELL, SPICA_TCELL_CONF
)
stopifnot(identical(
  unname(table(is.na(add.meta$SPICA_TCELL))[2]), ncol(se.t[, se.t$T_LIN_WORK == "gdT"])
))
add.meta = add.meta %>% dplyr::mutate(
  SPICA_TCELL = dplyr::case_when(
    is.na(SPICA_TCELL)  ~ T_LIN,
    TRUE ~ SPICA_TCELL
  )
)

se.meta.tmp = se.meta
# se.meta = se.meta.tmp

se.meta = AddMetaData(se.meta, add.meta)

pd = se.meta@meta.data
rm.barcodes = rownames(pd[is.na(pd$T_LIN) & pd$CT_L1_COARSE == "T-Cell", ])
se.meta = se.meta[, !colnames(se.meta) %in% rm.barcodes]

table(se.meta$T_LIN, se.meta$CT_L1_COARSE, useNA = "always")

# cycling cell will be annotate separately
se.meta$celltype = trimws(gsub(" Proliferating", "", se.meta$CT_L2))
se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate(
    celltype = dplyr::case_when(
      is.na(SPICA_TCELL) ~ celltype,
      TRUE ~ SPICA_TCELL
    )
  )

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    celltype_short_2 = case_when(
      grepl("^CD4", celltype) ~ "CD4 T-Cell",
      grepl("^CD8", celltype) ~ "CD8 T-Cell",
      grepl("gdT", celltype) ~ "gd T-Cell",
      grepl("dpT", celltype) ~ "dp T-Cell",
      TRUE ~ celltype
    )
  )

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    celltype_short_3 = case_when(
      grepl("Eryth", celltype) ~ "Erythrocyte",
      grepl("HSPC|CLP|EMP|GMP|LMPP|HSC|Prog_Mk|Prog Mk|Prog_DC|^Prog_B|pro B|pre B|BaEoMa", celltype) ~ "Progenitor",
      grepl("ILC|BaEoMa|Stromal", celltype) ~ "Other",
      grepl("B cell|Memory B|B memory|B-Cell Memory|Naive B|B-Cell Naive|B naive|transitional B|B intermediate", celltype) ~ "B-Cell",
      grepl("Plasma|Plasmablast", celltype) ~ "Plasma cell",
      grepl("^cDC", celltype) ~ "cDC",
      grepl("^pDC", celltype) ~ "pDC",
      grepl("^ASDC|^mDC|pre-pDC|pre-mDC", celltype) ~ "other DC",
      grepl("^CD14", celltype) ~ "Mono CD14",
      grepl("^CD16", celltype) ~ "Mono CD16",
      grepl("Macrophage", celltype) ~ "Macrophage",
      grepl("^NK|CD56 bright NK", celltype) ~ "NK",
      grepl("^CD4", celltype) ~ "CD4 T-Cell",
      grepl("^CD8", celltype) ~ "CD8 T-Cell",
      grepl("gdT", celltype) ~ "gd T-Cell",
      grepl("dpT", celltype) ~ "dp T-Cell",
      TRUE ~ celltype
    )
  )
# table(se.meta$celltype, se.meta$celltype_short_3)

pd = se.meta@meta.data
rm.ct = table(pd$celltype, pd$CellCycle)[, 1] < 100
se.meta = se.meta[, !se.meta$celltype %in% names(rm.ct[rm.ct])]
se.meta@meta.data = droplevels(se.meta@meta.data)

se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta@meta.data$celltype = factor(se.meta@meta.data$celltype)
se.meta@meta.data$celltype_short_2 = factor(se.meta@meta.data$celltype_short_2)
se.meta@meta.data$celltype_short_3 = factor(se.meta@meta.data$celltype_short_3)

se.meta$VDJ_T_AVAIL = se.meta.w$VDJ_T_AVAIL[match(rownames(se.meta@meta.data), rownames(se.meta.w@meta.data))]
se.meta$VDJ_B_AVAIL = se.meta.w$VDJ_B_AVAIL[match(rownames(se.meta@meta.data), rownames(se.meta.w@meta.data))]

se.meta

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate_if(is.character, as.factor)

saveRDS(se.meta, output.file)
