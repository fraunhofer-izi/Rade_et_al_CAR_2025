# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "data.table", "ggplot2",
  "scales", "openxlsx", "cowplot", "scCustomize"
)
.bioc_packages = c(
  "dittoSeq", "clustifyr", "scds", "scDblFinder", "ComplexHeatmap", "UCell"
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

if (any(!"Azimuth" %in% installed.packages())) {
  remotes::install_github('satijalab/azimuth', ref = 'master')
}
library(Azimuth)

if (any(!"Azimuth" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"SeuratDisk" %in% installed.packages())) {
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

# if (any(!"ProjecTILs" %in% installed.packages())) {
#   Sys.unsetenv("GITHUB_PAT")
#   remotes::install_github("carmonalab/ProjecTILs")
# }
library(ProjecTILs)

source("code/helper/functions.R")
source("code/helper/styles.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.b1 = readRDS(paste0(manifest$ukl_b1$seurat, "seurat_ori_pub.Rds"))
se.b2 = readRDS(paste0(manifest$ukl_b2$seurat, "seurat_ori.Rds"))
se.b3 = readRDS(paste0(manifest$ukl_b3$seurat, "seurat_post_souporcell_cilta_ida.Rds"))
se.b4 = readRDS(paste0(manifest$ukl_b4$seurat, "seurat_post_souporcell.Rds"))

m = list(
  a = rownames(se.b1[["ADT"]]),
  b = rownames(se.b2[["ADT"]]),
  c = rownames(se.b3[["ADT"]]),
  d = rownames(se.b3[["ADT"]])
)
m = ComplexHeatmap::list_to_matrix(m)

# Intersection set of ADT features
# Exception: hashtags (only run2, maybe later for QC)
keep.adt.ftrs = c(rownames(m[rowSums(m) == ncol(m), ]), "HASHTAG-1", "HASHTAG-2")

output.file.1 = paste0(manifest$meta_pub, "01_seurat_pre_low_qual_cell_filter.Rds")
output.file.2 = paste0(manifest$meta_pub, "02_seurat_post_low_qual_cell_filter.Rds")

# Batch 1
pdata.batch1 = read.xlsx(
  "data/clinicopathological_discovery.xlsx",
  sheet = 1, rowNames = F, startRow = 1
)

se.b1@meta.data$STUDY = "UKL_Batch_1"
se.b1@meta.data$SOURCE = pdata.batch1$source[match(se.b1@meta.data$orig.ident, pdata.batch1$SAMPLE_NAME)]
se.b1@assays$ADT@counts = se.b1@assays$ADT@counts[rownames(se.b1@assays$ADT) %in% keep.adt.ftrs, ]
se.b1@assays$ADT@data = se.b1@assays$ADT@counts

# Batch 2
se.b2@meta.data$STUDY = "UKL_Batch_2"
se.b2@meta.data$SOURCE = "PB"
se.b2@assays$ADT@counts = se.b2@assays$ADT@counts[rownames(se.b2@assays$ADT) %in% keep.adt.ftrs, ]
se.b2@assays$ADT@data = se.b2@assays$ADT@counts

# Batch 3
se.b3@meta.data$STUDY = "UKL_Batch_3"
se.b3@meta.data$SOURCE = "PB"
se.b3@assays$ADT@counts = se.b3@assays$ADT@counts[rownames(se.b3@assays$ADT) %in% keep.adt.ftrs, ]
se.b3@assays$ADT@data = se.b3@assays$ADT@counts

# Batch 4
se.b4@meta.data$STUDY = "UKL_Batch_4"
se.b4@meta.data$SOURCE = "PB"
se.b4@assays$ADT@counts = se.b4@assays$ADT@counts[rownames(se.b4@assays$ADT) %in% keep.adt.ftrs, ]
se.b4@assays$ADT@data = se.b4@assays$ADT@counts

table(se.b2$WELL_SPLIT, se.b2$status, useNA = "always")
table(se.b3$WELL_SPLIT, se.b3$status, useNA = "always")
table(se.b4$WELL_SPLIT, se.b4$status, useNA = "always")

se.meta = merge(se.b1, se.b2)
se.meta = merge(se.meta, se.b3)
se.meta = merge(se.meta, se.b4)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

# tmp = se.meta
# se.meta = tmp

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Harmonize CAR construct gene name
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mat.counts = GetAssayData(object = se.meta, assay = "RNA", slot = "counts")
car.exprs = mat.counts[c("ciltacel", "idecel"), ]
car.exprs = car.exprs[1,  , drop = F] + car.exprs[2,  , drop = F]
rownames(car.exprs) = "CAR-BCMA"

new.counts = rbind(
  car.exprs,
  mat.counts[!rownames(se.meta) %in% c("ciltacel", "idecel"), ]
)

adt.counts = GetAssayData(object = se.meta, assay = "ADT", slot = "counts")
se.meta <- CreateSeuratObject(counts = new.counts, meta.data = se.meta@meta.data)
se.meta[["ADT"]] = CreateAssayObject(counts = adt.counts)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = cd4cd8_car_present(se.meta)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = NormalizeData(se.meta, assay = 'RNA', normalization.method = "LogNormalize")
se.meta = NormalizeData(se.meta, assay = "ADT", normalization.method = 'CLR', margin = 2)
slot(object = se.meta[["ADT"]], name = 'data') = as(se.meta@assays$ADT@data, "dgCMatrix")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add %MT, %Ribosomal and complexity values
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^MT-")
se.meta[["Perc_of_ribosomal_genes"]] = Seurat::PercentageFeatureSet(se.meta, pattern = "^RPL|^RPS")
se.meta@meta.data$log10GenesPerUMI = log10(se.meta$nFeature_RNA) / log10(se.meta$nCount_RNA)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
Idents(se.meta) = "orig.ident"

saveRDS(se.meta, output.file.1)

# se.meta = readRDS(output.file.1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nFeature_low_cutoff = 250
nFeature_high_cutoff = 8000
nCount_low_cutoff = 1000
nCount_high_cutoff = 100000
mt_cutoff = 15
complx_cutoff = 0.8

label_cells_rm = function(obj) {
  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
      (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
      (Perc_of_mito_genes > mt_cutoff) | (log10GenesPerUMI < complx_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  # tbl = as.data.frame.matrix(table(obj$orig.ident, obj$KEEP_CELL))
  # tbl$PERC = round(tbl[, 2] / rowSums(tbl)* 100, 2)
  # print(tbl)
  obj
}

cell.track = count_cells_per_sample(c(se.meta))
se.meta = label_cells_rm(obj = se.meta)
se.meta = subset(se.meta, subset = KEEP_CELL == TRUE)
se.meta@meta.data = droplevels(se.meta@meta.data)
se.meta$KEEP_CELL = NULL
cell.track = count_cells_per_sample(c(se.meta), cell.track, "n1")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CellCycleScoring (Part 1)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("data/signatures/cellCycleMarkers.R")
se.meta = CellCycleScoring(
  se.meta, s.features = s.genes, g2m.features = g2m.genes,
  assay = 'RNA', search = TRUE
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj.l = Split_Object(se.meta, split.by = "orig.ident", threads = 20)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("scds")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scds_doublets = function(se){
  suppressWarnings({
    suppressMessages({
      sce = as.SingleCellExperiment(se)
      sce = scds::cxds(sce)
      sce = scds::bcds(sce)
      sce = scds::cxds_bcds_hybrid(sce, estNdbl = T)
      CD  = data.frame(sce@colData)
      CD = CD %>% dplyr::select(cxds_score, bcds_score, hybrid_score, cxds_call, bcds_call, hybrid_call)
      CD$barcode = rownames(CD)
      CD
    })
  })
}

scds_res = parallel::mclapply(obj.l, function(se){
  scds_doublets(se)
}, mc.cores = 20)

scds_res = do.call("rbind", scds_res)
rownames(scds_res) = scds_res$barcode
scds_res$barcode = NULL
se.meta = AddMetaData(se.meta, scds_res)
se.meta$DOUBLETS_CONSENSUS = se.meta$scDblFinder_class == "doublet" & se.meta$hybrid_call == TRUE

cell.track = count_cells_per_sample(c(subset(se.meta, DOUBLETS_CONSENSUS == FALSE)), cell.track, "n2")

saveRDS(
  cell.track, file = "data/qc/stats_celltrack.Rds"
)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CellCycleScoring (Part 2)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# se = obj.l$MXMERZ002A_03
estimate_cc = function(se){

  suppressWarnings({
    suppressMessages({

      data(cell.cycle.obj) # ProjecTILs package
      DefaultAssay(se) = "RNA"

      se = se %>%
        FindVariableFeatures(verbose = F) %>%
        ScaleData(verbose = F) %>%
        RunPCA(npcs = 20, verbose = F) %>%
        FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

      tmp =  tryCatch(
        FindClusters(se, resolution = 2, verbose = F),
        error=function(e) "error"
      )
      if(class(tmp) != "Seurat") {
        se = FindClusters(se, resolution = 1, verbose = F)
      } else {
        se = tmp
      }

      # if(ncol(se) < 31) {
      #   se = RunUMAP(se, reduction = "pca", dims = 1:20, seed.use = 1234, verbose = F)
      # } else {
      #   se = RunUMAP(se, reduction = "pca", dims = 1:20, seed.use = 1234, verbose = F)
      # }

      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }

      cc.phase = quiet(
        clustifyr::run_gsea(
          GetAssayData(object = se, assay = "RNA", slot = "data"),
          query_genes = cell.cycle.obj$human$cycling,
          cluster_ids =  se@meta.data[["seurat_clusters"]], n_perm = 1000
        )
      )

      cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
      cc.cl = rownames(cc.phase[cc.phase$pval < 0.05, ])

      cc.cells = (se@meta.data[["seurat_clusters"]] %in% cc.cl)
      se@meta.data$CellCycle = factor(cc.cells)

      pd.cc = se@meta.data[cc.cells, ]
      pd.cc = pd.cc %>% dplyr::mutate(
        CellCycle_Phase = dplyr::case_when(
          G2M.Score > S.Score ~ "G2M",
          TRUE ~ "S"
        )
      )
      se$CellCycle_Phase = pd.cc$CellCycle_Phase[match(rownames(se@meta.data), rownames(pd.cc))]
      se$CellCycle_Phase[is.na(se$CellCycle_Phase)] = "G1M"
      se$CellCycle_Phase = factor(se$CellCycle_Phase, levels = c("G1M", "S", "G2M"))

      # (DimPlot_scCustom(se, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) & mytheme() & theme(legend.position = "none") |
      #   DimPlot_scCustom(se, reduction = "umap", group.by = "CellCycle", pt.size = 1) & mytheme() |
      #   DimPlot_scCustom(se, reduction = "umap", group.by = "CellCycle_Phase", pt.size = 1) & mytheme() ) /
      #   (
      #   FeaturePlot_scCustom(
      #     se,
      #     reduction = "umap",
      #     features = c("S.Score", "G2M.Score"),
      #     pt.size = .5, na_cutoff = 0.1,
      #     colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100))
      #   ) & mytheme())

      res = se@meta.data %>% dplyr::select(CellCycle, CellCycle_Phase)
      res$barcode = rownames(res)
      res
    })
  })
}

cc.res = parallel::mclapply(obj.l, function(se){
  estimate_cc(se)
}, mc.cores = 7)
cc.res = do.call("rbind", cc.res)
rownames(cc.res) = cc.res$barcode
se.meta = AddMetaData(se.meta, cc.res)
se.meta$barcode = NULL

se.meta@meta.data = se.meta@meta.data %>%
  dplyr::mutate_if(is.character, as.factor)

# Annotate cells that were missed by the cluster approach
cutoff = .25
se.meta$CellCycle_Phase = dplyr::case_when(
  se.meta$S.Score > cutoff & se.meta$G2M.Score <= cutoff ~ "S",
  se.meta$G2M.Score > cutoff & se.meta$S.Score <= cutoff ~ "G2M",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$S.Score > se.meta$G2M.Score) ~ "S",
  (se.meta$S.Score > cutoff & se.meta$G2M.Score > cutoff) & (se.meta$G2M.Score > se.meta$S.Score) ~ "G2M",
  TRUE ~ se.meta$CellCycle_Phase
)
se.meta$CellCycle_Phase = factor(se.meta$CellCycle_Phase, levels = c("G1M", "S", "G2M"))
se.meta@meta.data$CellCycle = factor(ifelse(se.meta$CellCycle_Phase == "G1M", "FALSE", "TRUE"))

data(cell.cycle.obj) # ProjecTILs package
se.meta = AddModuleScore_UCell(
  se.meta,
  features = list("CellCycle_SCORE" = cell.cycle.obj$human$cycling),
  assay = "RNA", ncores = 25, force.gc = T
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
Idents(se.meta) = "orig.ident"

saveRDS(se.meta, output.file.2)
