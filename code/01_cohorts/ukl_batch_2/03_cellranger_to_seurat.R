print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Batch 2")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("Seurat", "yaml", "dplyr", "doParallel", "parallel", "data.table", "Matrix")
.bioc_packages = c("biomaRt", "scDblFinder", "SingleCellExperiment")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = manifest$ukl_b2$data_dl
seurat.path = paste0(manifest$base$workdata, "cohorts/ukl_batch_2/seurat/")
dir.create(seurat.path, recursive = T)

obj.path = paste0(seurat.path, "seurat_ori.Rds")

cellranger.dirs = list.dirs(path = work.path, full.names = T, recursive = T)

fltrd.dirs = cellranger.dirs[grepl("sample_filtered_feature_bc_matrix", cellranger.dirs)]
tmp = gsub(work.path, "", fltrd.dirs)
tmp = gsub("/", "", gsub("out.*", "", tmp))
names(fltrd.dirs) = tmp
print(length(fltrd.dirs))

# i = names(fltrd.dirs)[1]
bpparam = BiocParallel::MulticoreParam(workers = 5)
seurat.l = BiocParallel::bplapply(names(fltrd.dirs), function(i) {

  id = i
  fltrd.counts = Read10X(data.dir = fltrd.dirs[names(fltrd.dirs) == id], gene.column = 2)

  seu.obj = CreateSeuratObject(counts = fltrd.counts[[1]], project = id)
  seu.obj[["ADT"]] = CreateAssayObject(counts = fltrd.counts[[2]])

  seu.obj = RenameCells(seu.obj, new.names = gsub("multi_", "", colnames(seu.obj)))
  seu.obj@meta.data$orig.ident = gsub("multi_", "", id)

  # # scDblFinder
  # sce = scDblFinder(GetAssayData(seu.obj, slot="counts"))
  # stopifnot(identical(colnames(seu.obj), colnames(sce)))
  # seu.obj@meta.data$scDblFinder_score = sce$scDblFinder.score
  # seu.obj@meta.data$scDblFinder_class = sce$scDblFinder.class

  seu.obj

}, BPPARAM = bpparam)

se.meta = merge(seurat.l[[1]], y = seurat.l[2:length(seurat.l)], add.cell.ids = names(seurat.l), project = "ukl_batch2")
se.meta@meta.data$orig.ident = factor(se.meta@meta.data$orig.ident)
se.meta[["RNA"]] <- JoinLayers(se.meta[["RNA"]])

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assay ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "ADT"
adt.ftrs = gsub("totalseqC-", "", rownames(se.meta))

rownames(se.meta@assays$ADT@counts) = toupper(adt.ftrs)
se.meta@assays$ADT@data = se.meta@assays$ADT@counts
rownames(se.meta@assays$ADT@meta.features) = rownames(se.meta@assays$ADT@counts)
DefaultAssay(se.meta) = "RNA"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Souporcell
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
souporcell.path = paste0(manifest$ukl_b2$work, "souporcell/output")
souporcell.clusters = list.files(
  path = souporcell.path, pattern = "clusters.tsv", full.names = T, recursive = T
)

pdata.sing = read.csv("data/metadata_basic_run_2_pub.csv")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parse souporcell outpult
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "ADT"
se.meta = NormalizeData(se.meta, normalization.method = 'CLR', margin = 2)
slot(object = se.meta[["ADT"]], name = 'data') = as(se.meta@assays$ADT@data, "dgCMatrix")
hashtags.exprs = FetchData(se.meta, vars = c("HASHTAG-1", "HASHTAG-2"), slot = "data")
DefaultAssay(se.meta) = "RNA"

l = list()
for (i in souporcell.clusters) {
  lane = basename(dirname(i))
  lane = gsub("multi_", "", lane)
  df = read.table(i, header = T)
  df$barcode = paste0(lane, "_", df$barcode)
  df$lane = lane
  l[[i]] = df
}
souporcell = do.call("rbind", l); rownames(souporcell) = NULL

souporcell$hashtag_1 = hashtags.exprs$`HASHTAG-1`[match(souporcell$barcode, rownames(hashtags.exprs))]
souporcell$hashtag_2 = hashtags.exprs$`HASHTAG-2`[match(souporcell$barcode, rownames(hashtags.exprs))]

sc.multi = souporcell[!souporcell$lane %in% pdata.sing[pdata.sing$HASHTAG_NAME == "",]$LANE_NAME, ]

# Based on the singleron QC report: UC3CGX (T3_Lane_4) had insufficient cells post-processing and was aborted.
# Therefore, this well contains only one sample
sc.multi = sc.multi[!sc.multi$lane %in% "T3_Lane_4", ]

sc.multi = sc.multi[sc.multi$barcode %in% colnames(se.meta), ]

sc.multi.cp = sc.multi %>%
  dplyr::filter(assignment == "0" | assignment == "1") %>%
  dplyr::select(barcode, assignment, lane, hashtag_1, hashtag_2)
sc.multi.cp = reshape2::melt(sc.multi.cp, id = c("barcode", "assignment", "lane"))

sc.multi.ave = sc.multi.cp %>%
  dplyr::group_by(lane, assignment, variable) %>%
  dplyr::summarise(ave = median(value))

sc.multi.max = sc.multi.ave %>%
  dplyr::group_by(lane, assignment) %>%
  slice_max(ave, n = 1)

sc.multi.max$sample_name = pdata.sing$SAMPLE_NAME[match(
  paste0(sc.multi.max$lane, "_", sc.multi.max$variable),
  paste0(pdata.sing$LANE_NAME, "_", tolower(pdata.sing$HASHTAG_NAME))
)]

sc.multi.cp$sample_name = sc.multi.max$sample_name[match(
  paste0(sc.multi.cp$lane, "_", sc.multi.cp$assignment),
  paste0(sc.multi.max$lane, "_", sc.multi.max$assignment)
)]

sc.multi$sample_name = sc.multi.cp$sample_name[match(sc.multi$barcode, sc.multi.cp$barcode)]

df = sc.multi
df = subset(df, assignment == "1" | assignment == "0")

df = df %>% dplyr::select(hashtag_1, hashtag_2, sample_name, lane, assignment)
df = reshape2::melt(df ,id = c("sample_name", "lane", "assignment"))
axis.cut = quantile(c(df$value), .999)
df$value[df$value >= axis.cut] = axis.cut
df$variable = gsub("hashtag", "ht", df$variable)
df$sample_name = paste0(df$assignment, " | ", df$sample_name)
df$lane = gsub("Lane", "Well", df$lane)

# ggplot(df, aes(x = value, y = variable, fill = variable)) +
#   ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 0.5, size = .2, alpha = .5) +
#   theme(
#     # aspect.ratio = 1,
#     panel.spacing = unit(1, "lines"),
#     legend.position = "none"
#   ) +
#   xlab("Hashtag expression") +
#   facet_wrap( ~ lane + sample_name, ncol = 6)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TCR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fltrd.vdj = list(
  T5_Lane_1 = paste0(manifest$ukl_b2$data_dl, "multi_T5_Lane_1/outs/per_sample_outs/multi_T5_Lane_1/vdj_t/filtered_contig_annotations.csv"),
  T5_Lane_2 = paste0(manifest$ukl_b2$data_dl, "multi_T5_Lane_2/outs/per_sample_outs/multi_T5_Lane_2/vdj_t/filtered_contig_annotations.csv")

)

contig_list <- lapply(fltrd.vdj, function(x) {
  tryCatch(read.csv(x), error=function(e) NULL)
})

combined <- combineTCR(
  contig_list,
  samples = paste0(names(contig_list))
)

souporcell.sub = subset(souporcell, lane == "T5_Lane_1")
combined.sub = combined$T5_Lane_1[combined$T5_Lane_1$CTstrict %in% combined$T5_Lane_2$CTstrict, ]
souporcell.sub = souporcell.sub[souporcell.sub$barcode %in% combined.sub$barcode, ]
table(souporcell.sub$assignment)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Merge and save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pheno = se.meta@meta.data
pheno = pheno %>% dplyr::select(orig.ident)
pheno$barcode = rownames(pheno)
pheno$barcode_short = gsub(".+_", "", rownames(pheno))

sc.other = souporcell[!souporcell$barcode %in% sc.multi$barcode, ]
sc.t5.lane.1 = subset(sc.other, lane == "T5_Lane_1")
sc.other = subset(sc.other, lane != "T5_Lane_1")
sc.other$sample_name = pdata.sing$SAMPLE_NAME[match(sc.other$lane, pdata.sing$LANE_NAME)]

sc.t5.lane.1 = sc.t5.lane.1 %>% dplyr::mutate(
  sample_name = dplyr::case_when(
    assignment == "1" ~ "P23_1_PB_VL",
    assignment == "0" ~ "P15_1_PB_LP"
  )
)

sc.fin = rbind(
  sc.multi %>% dplyr::mutate(WELL_SPLIT = "multi"),
  sc.t5.lane.1 %>% dplyr::mutate(WELL_SPLIT = "multi"),
  sc.other %>% dplyr::mutate(WELL_SPLIT = "single")
)
sc.fin = sc.fin %>% dplyr::select(-lane, -hashtag_1, -hashtag_2)

se.meta.clean = se.meta
se.meta.clean@meta.data$WELL = se.meta.clean@meta.data$orig.ident
se.meta.clean@meta.data$barcode = rownames(se.meta.clean@meta.data)
se.meta.clean@meta.data = merge(
  se.meta.clean@meta.data,
  sc.fin, by.x = "barcode", by.y = "barcode"
)
se.meta.clean@meta.data$orig.ident = se.meta.clean@meta.data$sample_name
se.meta.clean@meta.data$orig.ident = factor(se.meta.clean@meta.data$orig.ident)
se.meta.clean@meta.data$sample_name = NULL
stopifnot(identical(se.meta.clean@meta.data$barcode, colnames(se.meta)))
rownames(se.meta.clean@meta.data) = se.meta.clean@meta.data$barcode
se.meta.clean@meta.data$barcode = NULL

stopifnot(identical(colnames(se.meta.clean), colnames(se.meta)))
stopifnot(identical(rownames(se.meta.clean@meta.data), rownames(se.meta@meta.data)))

DefaultAssay(se.meta.clean) = "RNA"
Idents(se.meta.clean) = "orig.ident"

se.meta.clean = se.meta.clean[, !is.na(se.meta.clean$orig.ident)]
se.meta.clean$RM = paste0(se.meta.clean$WELL_SPLIT, "_", se.meta.clean$status)
se.meta.clean = subset(se.meta.clean, RM != "multi_unassigned")
se.meta.clean$RM = NULL
se.meta.clean@meta.data = droplevels(se.meta.clean@meta.data)

# table(se.meta.clean$WELL_SPLIT, se.meta.clean$status)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scDblFinder
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj.l = Split_Object(se.meta.clean, split.by = "orig.ident", threads = 40)
scDb.res = parallel::mclapply(obj.l, function(se){
  set.seed(1234)
  sce = scDblFinder(GetAssayData(se, slot="counts"))
  df = data.frame(sce@colData) %>% dplyr::select(scDblFinder.score, scDblFinder.class)
  colnames(df) = c("scDblFinder_score", "scDblFinder_class")
  df$barcode = rownames(df)
  df
}, mc.cores = 50)
scDb.res = do.call("rbind", scDb.res)
rownames(scDb.res) = scDb.res$barcode
scDb.res$barcode = NULL
se.meta.clean = AddMetaData(se.meta.clean, scDb.res)
DefaultAssay(se.meta.clean) = "RNA"
Idents(se.meta.clean) = "orig.ident"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
saveRDS(se.meta.clean, file = obj.path)
