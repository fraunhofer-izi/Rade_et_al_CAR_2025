# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "cowplot", "data.table",
  "ggplot2", "ggthemes", "scGate", "patchwork", "Signac", "devtools"
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
library(Azimuth)

if (any(!"SeuratData" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"SeuratDisk" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (Sys.info()["nodename"] == "ribnode020") {
  ncores = 35
} else {
  ncores = 5
}

manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta_pub, "02_seurat_post_low_qual_cell_filter.Rds"))
se.meta[["RNA"]] <- as(se.meta[["RNA"]], Class = "Assay")

output.file = paste0(manifest$meta_pub, "03_seurat_anno_1.Rds")

scGate_models_DB = get_scGateDB("data/scGateDB", force_update = F)
Idents(se.meta) = "orig.ident"

# pd = se.meta@meta.data
# set.seed(123)
# s.1 = sample(unique(as.character(pd[pd$SOURCE == "BM", ]$orig.ident)), size = 2)
# s.2 = sample(unique(as.character(pd[pd$SOURCE == "PB", ]$orig.ident)), size = 2)
# se.meta = se.meta[, se.meta$orig.ident %in% c(s.1, s.2)]
# se.meta = subset(se.meta, downsample = 1000)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Azimuth")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
anno_azimuth = function(obj,  ref = "pbmcref") {

  run = function(query){

    options(future.globals.maxSize = 8000 * 1024^2)

    suppressWarnings({
      suppressMessages({
        res = RunAzimuth(
          query, reference = ref,
          annotation.levels = list("celltype.l1", "celltype.l2"),
          verbose = F
        )
      })
    })

    if (class(res) == "try-error") {
      query@meta.data$predicted.celltype.l1 = "error"
      query@meta.data$predicted.celltype.l1.score = NA
      query@meta.data$predicted.celltype.l2 = "error"
      query@meta.data$predicted.celltype.l2.score = NA
      query@meta.data$mapping.score = NA
      d = query@meta.data %>%
        dplyr::select(
          predicted.celltype.l1, predicted.celltype.l1.score,
          predicted.celltype.l2, predicted.celltype.l2.score,
          mapping.score
        )
      d$barcode = rownames(d)
      d
    } else {
      d = res@meta.data %>% dplyr::select(
        predicted.celltype.l1, predicted.celltype.l1.score,
        predicted.celltype.l2, predicted.celltype.l2.score,
        mapping.score
      )
      d$barcode = rownames(d)
      d
    }
  }

  obj.l = Split_Object(obj, split.by = "orig.ident", threads = 30)
  df = parallel::mclapply(obj.l, function(se){run(se)}, mc.cores = 40)

  d = df[names(obj.l)]
  d = do.call("rbind", d)
  rownames(d) = d$barcode
  d= d[rownames(obj@meta.data), ]
  stopifnot(identical(d$barcode, rownames(obj@meta.data)))

  obj@meta.data$CT_L1 = d$predicted.celltype.l1
  obj@meta.data$CT_L1_SCORE = d$predicted.celltype.l1.score
  obj@meta.data$CT_L2 = d$predicted.celltype.l2
  obj@meta.data$CT_L2_SCORE = d$predicted.celltype.l2.score
  obj@meta.data$CT_MAPPING_SCORE = d$mapping.score
  obj@meta.data = obj@meta.data %>% dplyr::mutate_if(is.character, as.factor)

  return(obj)
}

se.meta.pb = subset(se.meta, SOURCE == "PB")
se.meta.pb@meta.data = droplevels(se.meta.pb@meta.data)
se.meta.pb = anno_azimuth(obj = se.meta.pb)

se.meta.bm = subset(se.meta, SOURCE == "BM")
se.meta.bm@meta.data = droplevels(se.meta.bm@meta.data)
se.meta.bm = anno_azimuth(obj = se.meta.bm, ref = "bonemarrowref")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
azi.anno.l1 = rbind(
  se.meta.pb@meta.data %>% dplyr::select(CT_L1, CT_L1_SCORE),
  se.meta.bm@meta.data %>% dplyr::select(CT_L1, CT_L1_SCORE)
)
azi.anno.l2 = rbind(
  se.meta.pb@meta.data %>% dplyr::select(CT_L2, CT_L2_SCORE, CT_MAPPING_SCORE),
  se.meta.bm@meta.data %>% dplyr::select(CT_L2, CT_L2_SCORE, CT_MAPPING_SCORE)
)

se.meta = AddMetaData(se.meta, azi.anno.l1)
se.meta = AddMetaData(se.meta, azi.anno.l2)
se.meta$CT_L2_ORI = se.meta$CT_L2

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    CT_L2 = case_when(
      grepl("Eryth", CT_L2) ~ "Erythrocyte",
      grepl("Naive B|B naive", CT_L2) ~ "B-Cell Naive",
      grepl("Memory B|B memory", CT_L2) ~ "B-Cell Memory",
      grepl("Plasma|Plasmablast", CT_L2) ~ "Plasma cell",
      grepl("NK CD56|NK_CD56bright", CT_L2) ~ "NK CD56bright",
      grepl("^cDC", CT_L2) ~ "cDC",
      grepl("^pDC", CT_L2) ~ "pDC",
      grepl("^CD14", CT_L2) ~ "Mono CD14",
      grepl("^CD16", CT_L2) ~ "Mono CD16",
      grepl("^CD8 Effector", CT_L2) ~ "CD8 Effector",
      TRUE ~ CT_L2
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta.l = Split_Object(se.meta, split.by = "orig.ident", threads = 30)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MISC
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
models.list.1 <- scGate_models_DB$human$PBMC[c("Platelet", "Erythrocyte")]
models.list.2 <- scGate_models_DB$human$generic[c(
  "Immune", "Myeloid", "MoMacDC", "Megakaryocyte"
)]
models.list.3 <- scGate_models_DB$human$TME_HiRes[c("Macrophage")]
models.list.4 <- scGate_models_DB$human$PBMC[c("NK")]
models.list.5 <- scGate_models_DB$human$PBMC[c("PlasmaCell")]
models.list.6 <- scGate_models_DB$human$PBMC[c("Bcell")]
models.list.7 <- scGate_models_DB$human$PBMC[c("gdT")]
models.list.8 <- scGate_models_DB$human$generic[c("Tcell")]

sc_gating = function(obj, obj.l, model) {
  suppressWarnings({
    suppressMessages({
      # x = obj.l[[1]]
      obj.l = parallel::mclapply(obj.l, function(x) {
        x = scGate::scGate(
          x, model = model, assay = "RNA", slot = "data",
          output.col.name = "SCGATE", ncores = 1
        )
        df = x@meta.data[, grepl("^SCGATE", colnames(x@meta.data)), drop = F]
        colnames(df) = toupper(colnames(df))
        df$barcode = rownames(df)
        df
      }, mc.cores = 40)
      df = do.call("rbind", obj.l)
      rownames(df) = df$barcode
      df$barcode = NULL
      obj = AddMetaData(obj, df)
    })
  })
  obj
}
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.1)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = c(models.list.2, models.list.3))
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.4)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.5)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.6)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.7)
se.meta = sc_gating(obj = se.meta, obj.l = se.meta.l, model = models.list.8)

se.meta = UCell::AddModuleScore_UCell(
  se.meta,
  features = list("PLATELET" = c("NRGN", "PF4", "PPBP", "ITGA2B")),
  assay = "RNA", slot = "data",
  ncores = 50, force.gc = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate_if(is.character, as.factor)

se.meta
saveRDS(se.meta, output.file)

