# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig1: Survival analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
surv_plot = function(
    obj,
    group,
    legend.labs,
    leg.title,
    col.pal = c("#7B9AB6", "#9B740A"),
    cutoff = 15,
    x.breaks = 3,
    x_anno = -1.5,
    y_anno = .12,
    conf.int = T,
    font.size = 9,
    anno.size = 2.3,
    nr.size = 2.5,
    ...
) {

  surv.obj <<- obj

  formula = as.formula(paste0("Surv(PFS_M, PROGRESSION)~", group))
  survModel = survfit(formula, data=surv.obj)
  survModel$call$formula = formula

  surv.pl = survminer::ggsurvplot(
    survModel,
    conf.int = conf.int,
    conf.int.alpha = .4,
    legend.labs = legend.labs,
    legend.title = " ",
    legend = "right",
    xlab = "Time (months)",
    ylab = "PFS probability",
    title =  " ",
    censor.size = 3,
    censor.shape = 124,
    size = .7,
    ggtheme = mytheme(base_size = font.size),
    tables.theme = theme_cleantable(),
    risk.table.fontsize = nr.size,
    risk.table.title = " ",
    risk.table.pos = "in",
    risk.table = "absolute",
    palette = col.pal,
    break.x.by = x.breaks,
    xlim = c(0,cutoff)
  )

  surv.pl$plot = surv.pl$plot +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.justification = c(.5, 0),
      legend.direction = "horizontal",
      legend.margin = margin(5, 10, 10, 5),
      legend.box.margin=margin(-15,-15,-15,-15),
      legend.key.spacing.x = unit(10, "pt"),
      plot.margin =  margin(4, 4, -2, 4, unit = "pt"),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(vjust = + 2),
    ) +
    guides(colour = guide_legend(
      title.position="top",
      title.hjust = 0,
      # keywidth = 1.5,
      # keyheight = .7,
      override.aes = list(size = .6))
    ) +
    coord_cartesian(xlim=c(0, cutoff)) +
    scale_x_continuous(expand = c(0, 2), breaks = seq(0, cutoff, by = x.breaks))

  # Annotate KM plot
  log.rank.pval = survival::survdiff(formula, data = surv.obj)

  anno = paste0(
    "**", leg.title, "**",
    "<br>Log-rank p-value = ", signif(log.rank.pval$pvalue, 2),
    "<br>n=", sum(survModel$n), "; number of events: ", sum(survModel$n.event)
  )

  surv.pl$plot = surv.pl$plot + annotate(
    "richtext", x = x_anno, y = y_anno, hjust = 0, vjust=0, label = anno,
    size = anno.size, label.color=NA, fill=NA
  )

  surv.pl$table = surv.pl$table +
    theme(
      panel.border = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank()
    ) +
    coord_cartesian(xlim=c(0, cutoff)) +
    scale_x_continuous(expand = c(0, 2), breaks = seq(0, cutoff, by = x.breaks))

  surv.pl$plot + inset_element(
    surv.pl$table, left = -.01, bottom = -0.03, right = 1.0075, top =.3
  )

  #  return(surv.gg)
}


cox_table = function(
    obj,
    group,
    leg.txt = "cilta-cel vs. ide-cel") {

  tab_obj <<- obj

  survModel = as.formula(paste0("Surv(PFS_M, PROGRESSION) ~", group))
  nfit = coxph(survModel, data=tab_obj)

  eventText = capture.output(summary(nfit))[4]
  eventText = gsub("\\s+", " ", eventText)
  eventText = gsub("= ", " = ", eventText)

  tabtext = summary(nfit)$coefficients[, - (2:4)]
  tabtext = matrix(tabtext, ncol=2)

  target = tabtext[, 2, drop = F]
  l = list()
  for (i in 1:nrow(target)) {
    if (target[i, 1] <= 0.0001) {
      l[[i]] = formatC(target[i, 1], format = "e", digits= 2)
    } else if (target[i, 1] > 0.0001 & target[i, 1] <= 0.001) {
      l[[i]] = format(round(target[i, 1], digits=5), nsmall = 5)
    } else {
      l[[i]] = format(round(target[i, 1], digits=3), nsmall = 3)
    }
  }
  target = matrix(unlist(l), ncol = 1)

  tabtext = cbind(
    format(round(exp(tabtext[, 1, drop = F]), digits=2)),
    target
  )

  coi = (summary(nfit)$conf.int[, 3:4])
  coi = matrix(coi, ncol=2)
  coi = format(round(coi, digits=2))

  lhr = paste0(tabtext[, 1], " [", coi[, 1], ", ", coi[, 2], "]")
  tabtext[, 1] = lhr
  rownames(tabtext) = leg.txt
  colnames(tabtext) = c("  HR [95% CI]  ", "  p-value  ")

  ggtexttable(tabtext, theme = ttheme("light", base_size = 8))

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize, Harmony, Clustering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integration = function(
    obj,
    obj.l = NULL,
    no.ftrs = 500,
    .assay = "RNA",
    max.cl = 1,
    min.cells.per.sample = 25,
    threads = 5,
    .nbr.dims = 15,
    do.cluster = T,
    do.dimreduc = T,
    cc.regr = F,
    hvg.union = T,
    custom.features = NULL,
    run.integration = T,
    harmony.group.vars = NULL,
    obj.split.by = "orig.ident",
    .perp = 50,
    dmreduc.dims = NULL,
    n.neighbors = 30,
    min.dist = 0.3) {

  library(SignatuR)
  library(parallel)
  library(BiocParallel)
  library(harmony)

  if (any(!"readgmt" %in% installed.packages())) {
    Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("jhrcook/readgmt")
  }
  library(readgmt)

  start.time <- Sys.time()

  if(is.null(dmreduc.dims)) {
    dmreduc.dims = .nbr.dims
  }

  DefaultAssay(obj) = "RNA"
  obj@meta.data = droplevels(obj@meta.data)
  obj = DietSeurat(obj, counts = TRUE, data = TRUE)
  obj = NormalizeData(obj)

  # Gene categories to exclude from variable genes
  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito)[[1]],
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)[[1]],
    SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR)[[1]]
  )
  bl = c(bl, c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "XIST")) # gender genes
  bl <- unique(bl)

  if (hvg.union == T) {

    if(is.null(obj.l)) {
      print("Split object")
      obj.l = Split_Object(obj, split.by = obj.split.by, threads = threads)
    }

    select.bool = unlist(lapply(obj.l, function(x){ncol(x) >= min.cells.per.sample}))
    print(table(select.bool))
    obj.l = obj.l[select.bool]
    length(obj.l)

    print("HVG")
    obj.l = parallel::mclapply(obj.l, function(x) {
      x = x[!rownames(x) %in% bl, ]
      x = FindVariableFeatures(x, selection.method = "vst",  assay = .assay, verbose = FALSE)
      x
    }, mc.cores = threads)

    features = SelectIntegrationFeatures(object.list = obj.l, nfeatures = no.ftrs)

    VariableFeatures(obj) = features

    rm(obj.l); gc()

  } else if (!is.null(custom.features)) {
    VariableFeatures(obj) = custom.features
  } else {
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = no.ftrs, assay = .assay)
  }

  if (cc.regr == T) {
    obj = ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), assay = .assay)
  } else {
    obj = ScaleData(obj, assay = .assay)
  }

  obj = RunPCA(obj, assay = .assay)
  # plot(ElbowPlot(obj, ndims = 50))

  print(paste0("### PCs used for harmony clustering and dim reduc: ", .nbr.dims, " ###"))

  if(run.integration == T){
    obj = RunHarmony(
      obj, group.by.vars = harmony.group.vars, # theta = c(2,3),
      reduction.use ='pca', dims.use = 1:.nbr.dims, max_iter = 15, ncores = threads,
      verbose = T
    )
    comp.wrk = 'harmony'
  } else {
    comp.wrk = 'pca'
  }

  if (do.cluster == T) {
    print("Find clusters")
    obj = FindNeighbors(obj, reduction = comp.wrk, dims = 1:.nbr.dims, verbose = F)

    reso = seq(0,max.cl,.1)
    names(reso) = reso
    suppressWarnings({
      suppressMessages({
        findclusters.res = parallel::mclapply(reso, function(x) {
          FindClusters(obj, resolution = x, verbose = F)@meta.data[, "seurat_clusters", drop = F]
        }, mc.cores = length(reso))
      })
    })
    res.names = names(findclusters.res)
    findclusters.res = do.call("cbind", findclusters.res)
    colnames(findclusters.res) = paste0("RNA_snn_res.", res.names)
    stopifnot(identical(rownames(obj@meta.data), rownames(findclusters.res)))
    obj = AddMetaData(obj, findclusters.res)
  }
  if (do.dimreduc == T) {
    # print("tSNE")
    # obj = RunTSNE(
    #   obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
    #   nthreads = threads # tsne.method = "FIt-SNE"
    # )
    print("UMAP")
    obj = RunUMAP(
      obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
      min.dist = min.dist, n.neighbors = n.neighbors, verbose = F
    )
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  return(obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: GEX
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_gex = function(
    obj,
    .target = "celltype",
    group = "BEST_RESPONSE_CONSENSUS",
    ctrs.grp1 = "non-CR",
    ctrs.grp2 = "CR",
    split.by.tp = T,
    min.pct = .25,
    logfc.threshold = log2(1.25),
    .min.cells = 200,
    subsample = F,
    subsample.n = 250,
    test.method = "MAST",
    latent.vars = c("STUDY", "nFeature_RNA"),
    threads = 20,
    min.de.genes = NULL
){

  se.w = obj
  se.w@meta.data = droplevels(se.w@meta.data)

  se.w@meta.data[[.target]] = gsub("_", "\\.", se.w@meta.data[[.target]])
  if(split.by.tp == T) {
    se.w$TMP  = paste0(se.w@meta.data[[.target]], "_", se.w@meta.data$TIMEPOINT)
    tbl = as.data.frame.matrix(table(se.w$TMP, se.w@meta.data[[group]]))
    tbl = tbl[, colnames(tbl) %in% c(ctrs.grp1, ctrs.grp2)]
    tbl = tbl[rowSums2(tbl > .min.cells) == 2, ]
    print(tbl)

    # pd = se.w@meta.data
    # pd = pd[pd$TMP %in% rownames(tbl), ]
    # pd = pd[pd[[group]] %in% c(ctrs.grp1, ctrs.grp2), ]
    # for (i in naturalsort(unique(pd$TMP))) {
    #   d = droplevels(pd[pd$TMP %in% i, ])
    #   p = lapply(split(d, d[[group]]), function(x){
    #     x = droplevels(x)
    #     unname(table(x$orig.ident))
    #   })
    #   print("--------------------")
    #   print(i)
    #   print(p)
    # }

  } else {
    se.w$TMP = se.w@meta.data[[.target]]
    print(table(se.w@meta.data[[.target]], se.w@meta.data[[group]]))
  }

  res.dgea = run_wilx(
    obj = se.w, target = "TMP", min.cells = .min.cells,
    lfc.thresh = logfc.threshold, min.pct.thres = min.pct,
    contrast.group = group, contrast = c(ctrs.grp1, ctrs.grp2),
    rm.var.chains = T, subsample = subsample, subsample.n = subsample.n,
    test.method = test.method, latent.vars = latent.vars,
    threads = threads
  )

  res.dgea$celltype = gsub("_.+", "", res.dgea$cluster)
  res.dgea$timepoint = gsub(".+_", "", res.dgea$cluster)
  res.dgea.sign = subset(res.dgea, significant == T)
  res.dgea.sign = res.dgea.sign[abs(res.dgea.sign$avg_log2FC) > logfc.threshold, ]
  res.dgea.sign = res.dgea.sign[res.dgea.sign$pct.1 > min.pct | res.dgea.sign$pct.2 > min.pct, ]

  if(!is.null(min.de.genes)){
    ct.keep = names(table(res.dgea.sign$cluster)[table(res.dgea.sign$cluster) > min.de.genes])
    res.dgea.sign = res.dgea.sign[res.dgea.sign$cluster %in% ct.keep, ]
  }

  if(nrow(res.dgea.sign) == 0){
    return(print("No DE genes found"))
  }

  res.dgea$ID = paste0(res.dgea$cluster, "_", res.dgea$feature)
  res.dgea.sign$ID = paste0(res.dgea.sign$cluster, "_", res.dgea.sign$feature)
  res.dgea$significant = ifelse(
    res.dgea$ID %in% res.dgea.sign$ID, TRUE, FALSE
  )
  gc()
  print(table(res.dgea.sign$cluster, res.dgea.sign$group))
  list(res.dgea = res.dgea, res.dgea.sign = res.dgea.sign)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_adt = function(
    obj,
    .target = "celltype_short_3",
    group = "BEST_RESPONSE_CONSENSUS",
    ctrs.grp1 = "non-CR",
    ctrs.grp2 = "CR",
    split.by.tp = T,
    min.pct = .5,
    logfc.threshold = log2(1.25),
    .min.cells = 500,
    subsample = F,
    subsample.n = 250,
    .assay = "ADT",
    threads = 20
){

  se.w = obj
  DefaultAssay(se.w) = "ADT"
  se.w@meta.data = droplevels(se.w@meta.data)

  se.w = se.w[, se.w@meta.data[[group]] %in% c(ctrs.grp1, ctrs.grp2)]
  se.w@meta.data = droplevels(se.w@meta.data)

  se.w@meta.data[[.target]] = gsub("_", "\\.", se.w@meta.data[[.target]])
  if(split.by.tp == T) {
    se.w$TMP  = paste0(se.w@meta.data[[.target]], "_", se.w@meta.data$TIMEPOINT)
    tbl = as.data.frame.matrix(table(se.w$TMP, se.w@meta.data[[group]]))
    tbl = tbl[, colnames(tbl) %in% c(ctrs.grp1, ctrs.grp2)]
    tbl = tbl[rowSums2(tbl > .min.cells) == 2, ]
    print(tbl)

  } else {
    se.w$TMP = se.w@meta.data[[.target]]
    print(table(se.w@meta.data[[.target]], se.w@meta.data[[group]]))
  }

  fltr = rowSums2(table(se.w@meta.data$TMP, se.w@meta.data[[group]]) > .min.cells)
  se.w$TMP = ifelse(se.w$TMP %in% names(fltr[fltr == 2]), se.w$TMP, NA)
  se.w = se.w[, !is.na(se.w$TMP)]

  res.wlx = run_wilx(
    obj = se.w, target = "TMP", min.cells = .min.cells, assay = .assay,
    lfc.thresh = 0, min.pct.thres = 0, rm.var.chains = T,
    contrast.group = group, contrast = c(ctrs.grp1, ctrs.grp2),
    subsample = subsample, subsample.n = subsample.n,
    threads = threads
  )

  res.wlx$celltype = gsub("_.+", "", res.wlx$cluster)
  res.wlx$timepoint = gsub(".+_", "", res.wlx$cluster)
  if(nrow(res.wlx[res.wlx$pct.1 == 0 & res.wlx$pct.2 == 0, ]) != 0){
    res.wlx[res.wlx$pct.1 == 0 & res.wlx$pct.2 == 0, ]$avg_log2FC = 0
  }

  # split by celltype
  obj.l = Split_Object(se.w, split.by = .target, threads = threads)

  # x = obj.l$CD8.EM
  ave.stats = parallel::mclapply(obj.l, function(x){

    ave.stats = GetMatrixFromSeuratByGroupMulti(
      obj = x, features = rownames(x), .assay = .assay,
      group1 = group, group2 = "TIMEPOINT", perc_expr_thres = 1
    )

    ave.stats.grp1 = ave.stats$exp_mat[, grepl(paste0("^", ctrs.grp1), colnames(ave.stats$exp_mat)), drop = F] %>%
      reshape2::melt()
    colnames(ave.stats.grp1) = c("FTR", "GROUP", "ave.1")
    pct.stats.grp1 = ave.stats$percent_mat[, grepl(paste0("^", ctrs.grp1), colnames(ave.stats$percent_mat)), drop = F] %>%
      reshape2::melt()
    colnames(pct.stats.grp1) = c("FTR", "GROUP", "pct.1")
    ave.stats.grp1$pct.1 = pct.stats.grp1$pct.1[match(rownames(ave.stats.grp1), rownames(pct.stats.grp1))]
    ave.stats.grp2 = ave.stats$exp_mat[, grepl(paste0("^", ctrs.grp2), colnames(ave.stats$exp_mat)), drop = F] %>%
      reshape2::melt()
    pct.stats.grp2 = ave.stats$percent_mat[, grepl(paste0("^", ctrs.grp2), colnames(ave.stats$percent_mat)), drop = F] %>%
      reshape2::melt()
    ave.stats.grp1$ave.2 = ave.stats.grp2$value[match(rownames(ave.stats.grp1), rownames(ave.stats.grp2))]
    ave.stats.grp1$pct.2 = pct.stats.grp2$value[match(rownames(ave.stats.grp1), rownames(pct.stats.grp2))]
    ave.stats.grp1$TIMEPOINT = gsub(".+\\|", "", ave.stats.grp1$GROUP)
    ave.stats.grp1$group.1 = gsub("\\|.+", "", ave.stats.grp1$GROUP)

    if(split.by.tp == T) {
      ave.stats.grp1$FTR = paste0(
        ave.stats.grp1$FTR, "_", x@meta.data[[.target]][1], "_",
        ave.stats.grp1$TIMEPOINT
      )
    } else {
      ave.stats.grp1$FTR = paste0(
        ave.stats.grp1$FTR, "_", x@meta.data[[.target]][1]
      )
    }
    ave.stats.grp1
  }, mc.cores = 30)
  ave.stats = data.table::rbindlist(ave.stats)

  res.wlx$pct.1 = ave.stats$pct.1[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$pct.2 = ave.stats$pct.2[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$ave.1 = ave.stats$ave.1[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$ave.2 = ave.stats$ave.2[match(rownames(res.wlx), ave.stats$FTR)]

  res.wlx.sign = subset(res.wlx, significant == T)
  res.wlx.sign = res.wlx.sign[abs(res.wlx.sign$avg_log2FC) > logfc.threshold, ]
  res.wlx.sign = res.wlx.sign[res.wlx.sign$pct.1 > min.pct | res.wlx.sign$pct.2 > min.pct, ]
  res.wlx$ID = paste0(res.wlx$cluster, "_", res.wlx$feature)
  res.wlx.sign$ID = paste0(res.wlx.sign$cluster, "_", res.wlx.sign$feature)
  res.wlx$significant = ifelse(
    res.wlx$ID %in% res.wlx.sign$ID, TRUE, FALSE
  )
  print(table(res.wlx.sign$cluster, res.wlx.sign$group))
  list(res.wlx = res.wlx, res.wlx.sign = res.wlx.sign)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FindMarkers: GEX/ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_wilx = function(
    obj = NULL,
    target = NULL,
    contrast.group = NULL,
    contrast = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    rm.var.chains = F,
    padj.thresh = 0.05,
    lfc.thresh = .1,
    min.pct.thres = 0,
    test.method = "wilcox",
    latent.vars = NULL,
    subsample = F,
    subsample.n = 200,
    threads = 20
){

  DefaultAssay(obj) = assay

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)

  obj = obj[, obj@meta.data[[contrast.group]] %in% contrast]
  obj@meta.data = droplevels(obj@meta.data)

  print("Split object")
  obj.l = Split_Object(obj, split.by = target, threads = threads)

  tbl = table(obj@meta.data[[contrast.group]], obj@meta.data[[target]])

  # Filter out celltypes with less than x cells in one group
  obj.l = obj.l[colnames(tbl)[colSums(tbl >= min.cells) == 2]]

  options(warn = 1)
  bpparam = BiocParallel::MulticoreParam(workers = threads)

  print(paste0("Run DGEA: ", test.method))
  wil.res = BiocParallel::bplapply(names(obj.l), function(x) {

    o = obj.l[[x]]
    if(subsample == T){
      Idents(o) = "orig.ident"
      o = subset(o, downsample = subsample.n)
    }

    o@meta.data[[contrast.group]]
    markers = Seurat::FindMarkers(
      object = o,
      ident.1 = contrast[1], ident.2 = contrast[2],
      group.by = contrast.group, assay = assay,
      logfc.threshold = lfc.thresh, min.pct = min.pct.thres,
      test.use = test.method, latent.vars = latent.vars
    )

    markers$feature = rownames(markers)
    markers$cluster = x
    markers$group.1 = contrast[1]
    rownames(markers) = paste0(markers$feature, "_", markers$cluster)
    markers

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res)
  if (rm.var.chains == T) {
    wil.res = wil.res[!grepl('^IGHV|^IGK|^IGL|^IGL|^TCRAB|TRBC|^TR(G|B|A)V|TRGC|TRGV|JCHAIN|HASHTAG', wil.res$feature), ]
  }
  wil.res = wil.res[order(wil.res$p_val_adj, decreasing = F), ]
  wil.res$significant = (wil.res$p_val_adj < padj.thresh) & (abs(wil.res$avg_log2FC) > lfc.thresh)
  wil.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FindAllMarker with DOR
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
wlx_test_dor = function(
    obj = NULL,
    se.grp = NULL,
    downsample = TRUE,
    downsample.nbr = 1000,
    only.pos = T
) {

  obj@meta.data = droplevels(obj@meta.data)
  Idents(obj) = se.grp
  # res.wlx = FindAllMarkers(
  #   object = obj, only.pos = only.pos, min.diff.pct = 0, logfc.threshold = 0,
  #   return.thresh = 1, min.pct = 0
  # )
  # # table(res.wlx$cluster)
  # grps = as.character(unique(res.wlx$cluster))

  if (downsample == TRUE) {
    print(paste0("Downsample identity to ", downsample.nbr))
    obj.sub = subset(x = obj, downsample = downsample.nbr)
  } else {
    print("No Downsampling")
    obj.sub = obj
  }

  obj.sub@meta.data = droplevels(obj.sub@meta.data)
  Idents(obj.sub) = se.grp
  res.wlx = FindAllMarkers(
    object = obj.sub, only.pos = only.pos, min.diff.pct = 0, logfc.threshold = 0,
    return.thresh = 1, min.pct = 0
  )
  # table(res.wlx$cluster)
  grps = as.character(unique(res.wlx$cluster))

  # i = "0"
  l = list()
  for (i in grps) {

    grps.s = subset(res.wlx, cluster == i)

    grp.cells = rownames(obj.sub@meta.data %>% filter(.data[[se.grp]] == i))

    num.cellsInGroup = length(grp.cells)
    num.cellsOutGroup = nrow(obj.sub@meta.data) - length(grp.cells)

    # get number of cells within & outside the group with these genes
    num.TruePos = rowSums(obj.sub@assays$RNA@data[grps.s$gene , grp.cells] > 0)
    num.FalsePos = rowSums(obj.sub@assays$RNA@data[grps.s$gene, !colnames(obj.sub@assays$RNA@data) %in% grp.cells] > 0)

    # get number of cells in and outside group without genes
    num.FalseNeg = num.cellsInGroup - num.TruePos
    num.TrueNeg <- num.cellsOutGroup - num.FalsePos

    # use these values to calculate log(DOR) w/ a pseudocount of 0.5 to avoid +/- infinity values
    logDOR = log((num.TruePos+0.5)/(num.FalsePos+0.5)/((num.FalseNeg+0.5)/(num.TrueNeg+0.5)))

    stopifnot(identical(names(logDOR), grps.s$gene))

    grps.s$logDOR = logDOR
    l[[i]] = grps.s
  }

  if(nrow(res.wlx) > 0) {
    res = do.call("rbind", l) %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
    return(res)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA | Spectra
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_spectra = function(
    obj,
    .target = "celltype_short_3",
    group = "BEST_RESPONSE_CONSENSUS",
    ctrs.grp1 = "non-CR",
    ctrs.grp2 = "CR",
    split.by.tp = T,
    min.pct = .25,
    min.ave = 0.001,
    logfc.threshold = log2(1.25),
    .min.cells = 500,
    subsample = F,
    subsample.n = 200
){

  se.w = obj
  se.w@meta.data = droplevels(se.w@meta.data)

  se.w@meta.data[[.target]] = gsub("_", "\\.", se.w@meta.data[[.target]])
  if(split.by.tp == T) {
    se.w$TMP  = paste0(se.w@meta.data[[.target]], "_", se.w@meta.data$TIMEPOINT)
    tbl = as.data.frame.matrix(table(se.w$TMP, se.w@meta.data[[group]]))
    tbl = tbl[, colnames(tbl) %in% c(ctrs.grp1, ctrs.grp2)]
    tbl = tbl[rowSums2(tbl > .min.cells) == 2, ]
    print(tbl)

  } else {
    se.w$TMP = se.w@meta.data[[.target]]
    print(table(se.w@meta.data[[.target]], se.w@meta.data[[group]]))
  }

  res.wlx = run_wilx(
    obj = se.w, target = "TMP", min.cells = .min.cells,
    lfc.thresh = 0, min.pct.thres = 0,
    contrast.group = group, contrast = c(ctrs.grp1, ctrs.grp2),
    subsample = subsample, subsample.n = subsample.n
  )

  res.wlx$celltype = gsub("_.+", "", res.wlx$cluster)
  res.wlx$timepoint = gsub(".+_", "", res.wlx$cluster)
  res.wlx[res.wlx$pct.1 == 0 & res.wlx$pct.2 == 0, ]$avg_log2FC = 0

  if(split.by.tp == T){
    obj.l = Split_Object(se.w, split.by = .target, threads = 3)
  } else {
    obj.l = list(se.w)
  }
  ave.stats = parallel::mclapply(obj.l, function(x){
    ave.stats = GetMatrixFromSeuratByGroupMulti(
      obj = x, features = rownames(x),
      group1 = group, group2 = "TIMEPOINT", perc_expr_thres = min.ave
    )
    ave.stats.grp1 = ave.stats$exp_mat[, grepl(paste0("^", ctrs.grp1), colnames(ave.stats$exp_mat))] %>%
      reshape2::melt()
    colnames(ave.stats.grp1) = c("FTR", "GROUP", "ave.1")
    pct.stats.grp1 = ave.stats$percent_mat[, grepl(paste0("^", ctrs.grp1), colnames(ave.stats$percent_mat))] %>%
      reshape2::melt()
    colnames(pct.stats.grp1) = c("FTR", "GROUP", "pct.1")
    ave.stats.grp1$pct.1 = pct.stats.grp1$pct.1[match(rownames(ave.stats.grp1), rownames(pct.stats.grp1))]
    ave.stats.grp2 = ave.stats$exp_mat[, grepl(paste0("^", ctrs.grp2), colnames(ave.stats$exp_mat))] %>%
      reshape2::melt()
    pct.stats.grp2 = ave.stats$percent_mat[, grepl(paste0("^", ctrs.grp2), colnames(ave.stats$percent_mat))] %>%
      reshape2::melt()
    ave.stats.grp1$ave.2 = ave.stats.grp2$value[match(rownames(ave.stats.grp1), rownames(ave.stats.grp2))]
    ave.stats.grp1$pct.2 = pct.stats.grp2$value[match(rownames(ave.stats.grp1), rownames(pct.stats.grp2))]
    ave.stats.grp1$TIMEPOINT = gsub(".+\\|", "", ave.stats.grp1$GROUP)
    ave.stats.grp1$group.1 = gsub("\\|.+", "", ave.stats.grp1$GROUP)

    ave.stats.grp1$FTR = paste0(
      ave.stats.grp1$FTR, "_", x@meta.data[[.target]][1], "_",
      ave.stats.grp1$TIMEPOINT
    )
    ave.stats.grp1
  }, mc.cores = length(obj.l))
  ave.stats = data.table::rbindlist(ave.stats)

  res.wlx$pct.1 = ave.stats$pct.1[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$pct.2 = ave.stats$pct.2[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$ave.1 = ave.stats$ave.1[match(rownames(res.wlx), ave.stats$FTR)]
  res.wlx$ave.2 = ave.stats$ave.2[match(rownames(res.wlx), ave.stats$FTR)]

  res.wlx.sign = subset(res.wlx, significant == T)
  res.wlx.sign = res.wlx.sign[abs(res.wlx.sign$avg_log2FC) > logfc.threshold, ]
  res.wlx.sign = res.wlx.sign[res.wlx.sign$pct.1 > min.pct | res.wlx.sign$pct.2 > min.pct, ]
  res.wlx$ID = paste0(res.wlx$cluster, "_", res.wlx$feature)
  res.wlx.sign$ID = paste0(res.wlx.sign$cluster, "_", res.wlx.sign$feature)
  res.wlx$significant = ifelse(
    res.wlx$ID %in% res.wlx.sign$ID, TRUE, FALSE
  )
  print(table(res.wlx.sign$cluster, res.wlx.sign$group))
  list(res.wlx = res.wlx, res.wlx.sign = res.wlx.sign)
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4cd8_car_present = function(
    obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4", "CD274"), slot = "counts")
  ct.cd8cd4 = cd8cd4 %>% mutate(
    CD4CD8_BY_EXPRS = case_when(
      CD4 > 0 & CD8A == 0 & CD8B == 0 ~ "CD4+CD8-",
      CD4 == 0 & (CD8A > 0 | CD8B > 0) ~ "CD4-CD8+",
      CD4 == 0 & CD8A == 0 & CD8B == 0 ~ "CD4-CD8-",
      CD4 > 0 & (CD8A > 0 | CD8B > 0) ~ "CD4+CD8+",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD4CD8_BY_EXPRS)
  rownames(ct.cd8cd4) = rownames(cd8cd4)
  obj = AddMetaData(obj, ct.cd8cd4)

  cd3 = FetchData(obj, c("CD3D", "CD3E", "CD3G"), slot = "counts")
  cd3 = cd3 %>% mutate(
    CD3_BY_EXPRS = case_when(
      CD3D > 0 | CD3E > 0 | CD3G > 0 ~ "CD3",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD3_BY_EXPRS)
  obj = AddMetaData(obj, cd3)

  if("CAR-BCMA" %in% rownames(GetAssayData(obj, slot = c("counts"), assay = "RNA"))){
    car.ftr = FetchData(obj, c("CAR-BCMA"), slot = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr$`CAR-BCMA` > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }
  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Track number of cell (pre-processing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(STUDY, orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(STUDY, orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get metadata from Seurat
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_metadata <- function(obj, ..., embedding = names(obj@reductions), nbr.dim = 2) {

  res = as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (any(!embedding %in% names(obj@reductions))) {
      stop(paste0(embedding, " not found in seurat object\n"), call. = FALSE)
    }
    embed_dat = purrr::map(names(obj@reductions), ~obj@reductions[[.x]]@cell.embeddings[, 1:nbr.dim]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    res = dplyr::left_join(res, embed_dat, by = "cell")
  }

  if (length(list(...)) > 0) {
    cols_to_get <- setdiff(..., colnames(obj@meta.data))
    if (length(cols_to_get) > 0) {
      res = Seurat::FetchData(obj, vars = cols_to_get) %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::left_join(res, ., by = "cell")
    }
  }
  res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split Seurat object (BiocParallel)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  library(parallel)
  library(BiocParallel)

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# flag levels of significance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_signif <- function(
    data, p.col = NULL, output.col = NULL,
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,  1),
    symbols = c("****", "***", "**", "*",  "ns"),
    pval.relax = F
){

  if(pval.relax == T) {
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
    symbols = c("*****", "****", "***", "**", "*",  "")
  }

  if(is.null(output.col)) {
    output.col <- paste0(p.col, ".signif")
  }
  .p.values <- data %>% pull(!!p.col)
  if(all(is.na(.p.values))) {
    .p.signif <- rep("", length(.p.values))
  }
  else{
    .p.signif <- .p.values %>%
      stats::symnum(cutpoints = cutpoints, symbols = symbols, na = "") %>%
      as.character()
  }
  data %>%
    dplyr::mutate(!!output.col := .p.signif)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Write DGEA results to xlsx
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_to_xlsx = function(
    df = NULL,
    sheet = "sheet",
    filename = "Supplementaltable3.xlsx",
    wb = NULL
) {
  library("openxlsx")

  if(is.null(wb)){
    wb <- createWorkbook()
  }

  addWorksheet(wb, sheet)
  writeData(
    wb, sheet,
    df %>% dplyr::select(
      "Gene_symbol" = feature, "LogFC" = logFC,
      "Pvalue" = pval, "FDR" = padj, "Cell_identity" = cluster
    ) %>%
      dplyr::arrange(Cell_identity, desc(LogFC)),
    startRow = 1, startCol = 1)
  saveWorkbook(wb, filename, overwrite = T)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# standardize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
standardize = function(z) {
  rowmean = apply(z, 1, mean, na.rm=TRUE)
  rowsd = apply(z, 1, sd, na.rm=TRUE)
  rv = sweep(z, 1, rowmean,"-")
  rv = sweep(rv, 1, rowsd, "/")
  return(rv)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ave Expression and % of present cells
# https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
GetMatrixFromSeuratByGroupMulti<- function(
    obj,
    features,
    group1,
    group2,
    .assay = "RNA",
    perc_expr_thres = 0
){
  exp_mat<- obj@assays[[.assay]]@data[features, ,drop=FALSE]
  if(.assay == "ADT") {
    count_mat<- obj@assays[[.assay]]@data[features,,drop=FALSE ]
  } else {
    count_mat<- obj@assays[[.assay]]@counts[features,,drop=FALSE ]
  }

  if(group1 == group2){

    colnames(obj@meta.data)[colnames(obj@meta.data) == group2] = paste0(group2, "_2")
    group2 = paste0(group2, "_2")
    obj@meta.data[[group1]] = obj@meta.data[[group2]]
  }

  meta<- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")

  # get the average expression matrix
  exp_df <- as.matrix(exp_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    dplyr::left_join(meta) %>%
    dplyr::group_by(gene, .data[[group1]], .data[[group2]]) %>%
    dplyr::summarise(average_expression = mean(expression)) %>%
    # the trick is to make the data wider in columns: cell_type|group
    tidyr::pivot_wider(names_from = c(.data[[group1]], .data[[group2]]),
                       values_from= average_expression,
                       names_sep="|")

  # convert to a matrix
  exp_mat<- exp_df[, -1] %>% as.matrix()
  rownames(exp_mat)<- exp_df$gene

  # get percentage of positive cells matrix
  count_df <- as.matrix(count_mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    dplyr::left_join(meta) %>%
    dplyr::group_by(gene, .data[[group1]], .data[[group2]]) %>%
    dplyr::summarise(percentage = mean(count > perc_expr_thres)) %>%
    tidyr::pivot_wider(names_from = c(.data[[group1]], .data[[group2]]),
                       values_from= percentage,
                       names_sep="|")

  percent_mat<- count_df[, -1] %>% as.matrix()
  rownames(percent_mat)<- count_df$gene

  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }

  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }

  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# clean Clonotypes: Clonotypes must not have CD4 and CD8 cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clean_clonotypes = function(
    obj = NULL,
    celltype = "celltype_short_3",
    clone_id = "CTstrict"
){

  obj$celltype = obj[[celltype]]
  obj$clone_id = obj[[clone_id]]

  pd = obj@meta.data %>%
    dplyr::select(celltype, clone_id)
  pd$barcode = rownames(pd)
  cl.lin.diff =
    pd %>%
    dplyr::filter(!is.na(clone_id)) %>%
    dplyr::filter(grepl("CD4|CD8", celltype)) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(clone_size = n()) %>%
    dplyr::mutate(cd4 = grepl("CD4", celltype)) %>%
    dplyr::mutate(cd8 = grepl("CD8", celltype)) %>%
    dplyr::mutate(nbr_cd4 = sum(cd4)) %>%
    dplyr::mutate(nbr_cd8 = sum(cd8)) %>%
    dplyr::filter(nbr_cd4 > 0 & nbr_cd8 > 0)

  rm.cl.1 = cl.lin.diff %>%
    dplyr::filter(clone_size <= 5) %>%
    dplyr::pull(barcode) %>%
    unique()

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.1, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1))
  }

  cl.lin.diff$ratio = NA
  for (i in unique(cl.lin.diff$clone_id)) {
    cl.sub = cl.lin.diff[cl.lin.diff$clone_id == i, ]
    cl.sub = cl.sub[!duplicated(cl.sub$clone_id), ]
    .max = names(which.max(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    .min = names(which.min(cl.sub[, c("nbr_cd4", "nbr_cd8")]))
    cl.lin.diff[cl.lin.diff$clone_id == i, ]$ratio = cl.sub[[.max]] / cl.sub[[.min]]
  }

  rm.cl.2 = cl.lin.diff %>%
    dplyr::filter(ratio <= 3) %>%
    data.frame() %>%
    dplyr::pull(barcode) %>%
    unique()

  # tmp = cl.lin.diff[cl.lin.diff$barcode %in% rm.cl.2, ]
  # tmp[!duplicated(tmp$clone_id), ] %>% data.frame() %>% select(nbr_cd4, nbr_cd8, ratio)

  cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.2, ]

  if(nrow(cl.lin.diff) == 0){
    return(c(rm.cl.1, rm.cl.2))
  }

  # cl.lin.diff[!duplicated(cl.lin.diff$clone_id), ] %>% data.frame() %>% select(nbr_cd4, nbr_cd8, ratio)

  rm.cl.3 = cl.lin.diff %>%
    dplyr::group_by(clone_id, celltype) %>%
    dplyr::mutate(lin = n()) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::filter(lin == min(lin)) %>%
    dplyr::pull(barcode) %>%
    unique()

  return(c(rm.cl.1, rm.cl.2, rm.cl.3))
  # cl.lin.diff = cl.lin.diff[!cl.lin.diff$barcode %in% rm.cl.3, ]
  # cl.lin.diff %>%
  #   dplyr::group_by(clone_id) %>%
  #   dplyr::mutate(clone_size = n()) %>%
  #   dplyr::mutate(cd4 = grepl("CD4", celltype)) %>%
  #   dplyr::mutate(cd8 = grepl("CD8", celltype)) %>%
  #   dplyr::mutate(nbr_cd4 = sum(cd4)) %>%
  #   dplyr::mutate(nbr_cd8 = sum(cd8)) %>%
  #   dplyr::filter(nbr_cd4 > 0 & nbr_cd8 > 0) %>%
  #   data.frame()
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# A clonotype must not be present in more than one patient.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clean_clonotypes_inter = function(
    obj = NULL,
    clone_id = "CTstrict"
){

  obj$clone_id = obj[[clone_id]]
  pd = obj@meta.data

  clone.fltr = pd %>%
    dplyr::filter(!is.na(clone_id)) %>%
    dplyr::select(PATIENT_ID, clone_id) %>%
    dplyr::group_by(PATIENT_ID, clone_id) %>%
    dplyr::summarise(clone_size = n()) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(clone_present = n()) %>%
    dplyr::filter(clone_present > 1) %>%
    dplyr::arrange(clone_id) %>%
    data.frame()

  rm.cl.1 = clone.fltr %>%
    dplyr::filter(clone_present > 2) %>%
    dplyr::pull(clone_id) %>% unique()

  rm.cl.2 = clone.fltr %>%
    dplyr::group_by(clone_id) %>%
    dplyr::slice_max(clone_size, n = 1) %>%
    dplyr::filter(clone_size < 25) %>%
    dplyr::pull(clone_id) %>% unique()
  # rownames(pd[pd[[clone_id]] %in% unique(c(rm.cl.1, rm.cl.2)), ])

  rm.cl.3 = clone.fltr %>%
    dplyr::filter(clone_present == 2) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::slice_min(clone_size ) %>%
    data.frame()

  a = paste0(pd$PATIENT_ID, pd$clone_id )
  b = paste0(rm.cl.3$PATIENT_ID, rm.cl.3$clone_id)
  # rownames(pd[a %in% b, ])

  return(
    unique(
      c(
        rownames(pd[pd$clone_id %in% unique(c(rm.cl.1, rm.cl.2)), ]),
        rownames(pd[a %in% b, ])
      )
    )
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assign VDJ
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
assign_vdj = function(
    obj,
    vdj = "vdj_t",
    batch = NULL,
    present.bool = TRUE,
    export.table = FALSE,
    return.combined = F
){

  stopifnot(vdj %in% c("vdj_t", "vdj_b"))

  library(scRepertoire)

  l = list()
  for (i in batch) {
    print(i)
    cellranger.dirs = list.dirs(
      path = i, full.names = T, recursive = F
    )
    cellranger.samples = basename(cellranger.dirs)
    fltrd.vdj = paste0(
      cellranger.dirs, "/outs/per_sample_outs/", cellranger.samples, "/", vdj,
      "/filtered_contig_annotations.csv"
    )
    names(fltrd.vdj) = gsub("multi_", "", cellranger.samples)
    print(length(fltrd.vdj))
    l[[i]] = fltrd.vdj
  }
  names(l) = NULL
  fltrd.vdj = unlist(l, use.names = T)

  ###

  contig_list <- lapply(fltrd.vdj, function(x) {
    tryCatch(read.csv(x), error=function(e) NULL)
  })
  length(contig_list)
  print(paste0("Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])))
  contig_list = contig_list[lengths(contig_list) != 0]

  if(vdj == "vdj_t"){
    combined <- scRepertoire::combineTCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  } else {
    combined <- scRepertoire::combineBCR(
      contig_list,
      samples = paste0(names(contig_list))
    )
  }

  if(return.combined == T){
    return(combined)
  }


  combined = data.table::rbindlist(combined) %>% data.frame()

  if(present.bool == TRUE){
    if(vdj == "vdj_t"){
      obj$VDJ_T_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_T_AVAIL = ifelse(is.na(obj$VDJ_T_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_T_AVAIL))
    } else {
      obj$VDJ_B_AVAIL = combined$CTnt[match(rownames(obj@meta.data), combined$barcode)]
      obj$VDJ_B_AVAIL = ifelse(is.na(obj$VDJ_B_AVAIL), FALSE, TRUE)
      print(table(obj$VDJ_B_AVAIL))
    }
    obj
  } else {

    combined$sample = obj$orig.ident[match(combined$barcode, rownames(obj@meta.data))]
    combined = combined[!is.na(combined$sample), ]
    combined = split(combined, combined$sample)

    combined.keep = lapply(combined, function(x){nrow(x)}) > 0
    combined = combined[combined.keep]

    max.clonotypes = max(unlist(lapply(combined, function(x){max(unname(table(x$CTstrict)))})))
    obj <- combineExpression(
      combined, obj,
      cloneCall = "strict",
      group.by = "sample",
      proportion = F,
      cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=max.clonotypes)
    )
    obj$cloneSize = droplevels(obj$cloneSize)
    # obj
    if(export.table == T){
      obj@meta.data[, (length(colnames(obj@meta.data))-6):length(colnames(obj@meta.data))]
    } else{
      obj
    }
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# UCell: Gene set enrichment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ucell_enrich = function(se.w = NULL, cust.gene.sets = NULL, ncores = 20){

  library(UCell)

  se.w = AddModuleScore_UCell(
    se.w,
    features = cust.gene.sets,
    assay = "RNA", slot = "data",
    ncores = ncores, force.gc = T
  )
  ucell.res = se.w@meta.data[, grepl("UCell$", colnames(se.w@meta.data), ignore.case = T), drop = F]

  # ucell.res.smooth = UCell::SmoothKNN(
  #   obj = se.w,
  #   signature.names = colnames(ucell.res),
  #   reduction="pca", k=10
  # )
  # ucell.res.smooth = ucell.res.smooth@meta.data[, grepl("_UCell_kNN", colnames(ucell.res.smooth@meta.data))]

  colnames(ucell.res) = gsub("_UCell", "", colnames(ucell.res))
  ucell.res = ucell.res[, names(cust.gene.sets), drop = F]

  # colnames(ucell.res.smooth) = gsub("_UCell_kNN", "", colnames(ucell.res.smooth))
  # ucell.res.smooth = ucell.res.smooth[, colnames(ucell.res), drop = F]

  stopifnot(identical(rownames(ucell.res), colnames(se.w)))
  se.w[['custom_UCell_score']] = Seurat::CreateAssayObject(t(ucell.res))
  # slot(object = se.w[["custom_UCell_score"]], name = 'counts') <- new(Class = 'matrix')
  # se.w[['custom_UCell_smooth_score']] = Seurat::CreateAssayObject(t(ucell.res.smooth))
  # # slot(object = se.w[["custom_UCell_smooth_score"]], name = 'counts') <- new(Class = 'matrix')

  se.w@meta.data = se.w@meta.data[ , !grepl("UCell|nCount_UCell|nFeature_UCell", colnames(se.w@meta.data), ignore.case = T)]
  se.w

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Wlx for spectra results
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
spectra_wlx = function(
    df,
    .target,
    min.smpls = 4
){

  if(.target == "TIMEPOINT"){
    m = 3
  } else{
    m = 2
  }

  df$FACTOR_TP = paste0(df$FACTOR, "_", df$TIMEPOINT)
  rm.f =
    df %>%
    dplyr::group_by(FACTOR_TP, .data[[.target]]) %>%
    dplyr::count() %>%
    dplyr::mutate(min_smpls = n > min.smpls) %>%
    dplyr::group_by(FACTOR_TP) %>%
    dplyr::summarise(min_smpls = sum(min_smpls)) %>%
    dplyr::filter(min_smpls != m) %>%
    dplyr::pull(FACTOR_TP)
  df = df[!df$FACTOR_TP %in% rm.f, ]

  bpparam = BiocParallel::MulticoreParam(workers = 20)
  # i = "83-X-CD4_T-X-TNK_IL2-STAT5-signaling"
  wil.res = BiocParallel::bplapply(as.character(unique(df$FACTOR)), function(x) {

    # print(x)
    df.sub = df[df$FACTOR == x, ]

    if(.target == "TIMEPOINT") {
      fc = aggregate(df.sub$SCORE, list(df.sub$TIMEPOINT), median)
      colnames(fc) = c("TIMEPOINT", "mean_exp")
      fc = fc %>% dplyr::group_by(TIMEPOINT) %>%
        tidyr::pivot_wider(names_from = all_of(.target), values_from = mean_exp) %>%
        data.frame()
      fc =  cbind(fc[, 1, drop = F ], log2(fc[, 2:ncol(fc)] + 0.000001))
      ctrst.l = combn(colnames(fc)[1:length(colnames(fc))], 2, simplify = F)
    } else {
      fc = aggregate(df.sub$SCORE, list(df.sub$TIMEPOINT, df.sub[[.target]]), median)
      colnames(fc) = c("TIMEPOINT", .target, "mean_exp")
      fc = fc %>% dplyr::group_by(TIMEPOINT) %>%
        tidyr::pivot_wider(names_from = all_of(.target), values_from = mean_exp) %>%
        data.frame()
      fc =  cbind(fc[, 1, drop = F ], log2(fc[, 2:ncol(fc)] + 0.000001))
      if(colnames(fc)[2] == "CR" & colnames(fc)[3] == "non.CR"){
        fc = fc[, c(1, 3, 2)]
      }
      ctrst.l = combn(colnames(fc)[2:length(colnames(fc))], 2, simplify = F)

    }

    l = list()
    # i = 1
    for (i in 1:nrow(fc)) {
      max.l = list()
      for (j in 1:length(ctrst.l)) {
        max.l[[j]] = fc[i, ctrst.l[[j]][1]] - fc[i, ctrst.l[[j]][2]]
      }
      l[[i]] = unlist(max.l)[which.max(abs(unlist(max.l)))]
    }

    if(.target == "TIMEPOINT") {
      res = rstatix::kruskal_test(
        data = df.sub,
        as.formula(paste0("SCORE ~ ", "TIMEPOINT"))
      )
      res$FACTOR = x
    } else if(length(unique(as.character(df.sub[[.target]]))) > 2) {
      res = df.sub %>%
        dplyr::group_by(TIMEPOINT) %>%
        rstatix::kruskal_test(
          data =.,
          as.formula(paste0("SCORE ~ ", .target))
        )
      res$FACTOR = x
    } else if(length(unique(as.character(df.sub[[.target]]))) == 2) {
      res = df.sub %>%
        dplyr::group_by(TIMEPOINT) %>%
        rstatix::wilcox_test(
          data = .,
          as.formula(paste0("SCORE ~ ", .target))
        )
      res$FACTOR = x
    }

    res$LFC_MAX = unlist(l)

    res

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res) %>% data.frame()
  wil.res$FACTOR_TP = paste0(wil.res$FACTOR, "_", wil.res$TIMEPOINT)
  wil.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DA | late vs. lp, very late vs. lp | paired desing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
da_paired = function(
    obj = NULL,
    transform = "asin",
    only.paired = F,
    block = "PATIENT_ID"
){

  l = list()
  # i = "CR"
  for (i in names(obj)) {
    print(i)
    x = obj[[i]]
    x = droplevels(x[x$TIMEPOINT %in% c(grp1, grp2), ])
    paired = table(x[!duplicated(x$orig.ident), ]$PATIENT_ID)
    print(table(paired))
    if(only.paired == T){
      x = x[x$PATIENT_ID %in% names(paired[paired == 2]), ]
    }
    x = droplevels(x)

    clusters <- as.character(x$celltype)
    sample <- as.character(x$orig.ident)

    # Get transformed proportions
    prop.list <- getTransformedProps(clusters, sample, transform)

    pheno = x[!duplicated(x$orig.ident), ]
    rownames(pheno) = pheno$orig.ident
    pheno = pheno[colnames(prop.list$TransformedProps), ]
    stopifnot(identical(colnames(prop.list$Counts), rownames(pheno)))

    # tmp = prop.list$TransformedProps
    # tmp = reshape2::melt(tmp)
    # tmp$TP = pheno$TIMEPOINT[match(tmp$sample, pheno$orig.ident)]
    # tmp$PATIENT_ID = pheno$PATIENT_ID[match(tmp$sample, pheno$orig.ident)]
    # ggplot(tmp[tmp$clusters == "CD4.CTL_GNLY", ], aes(TP, value)) +
    #   geom_boxplot() +
    #   geom_point() +
    #   geom_line(aes(group=PATIENT_ID))

    if(!is.null(block)) {
      covars = pheno %>% dplyr::select(ctrst = TIMEPOINT, block = all_of(block))
      design = model.matrix(~ 0 + ctrst + block, data = covars)
      colnames(design) = gsub("ctrst", "", colnames(design))
      mycontr <- makeContrasts(paste0(grp1, " - ", grp2), levels=design)
    } else {
      groups = pheno %>% dplyr::select(ctrst = all_of(ctrst)) %>% data.frame()
      design = model.matrix(~ 0 + ctrst, data = covars)
      colnames(design) = gsub("ctrst", "", colnames(design))
      mycontr <- makeContrasts(paste0(grp1, " - ", grp2), levels=design)
    }

    # out <- speckle::propeller.ttest(
    #   prop.list, design, contrasts = mycontr,
    #   robust=T, trend=F, sort=T
    # )
    # out

    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    # Add check for fewer than 3 cell types
    # Robust eBayes doesn't work with fewer than 3 cell types
    if(nrow(prop.trans)<=2){
      message("Setting robust to FALSE for eBayes for less than 3 cell types")
      robust <- FALSE
    }

    fit <- lmFit(prop.trans, design)
    # gene_coef = data.frame(design %*% fit$coefficients["CD8.NaiveLike", ])
    # colnames(gene_coef) = "value"
    # gene_coef$PATIENT = x$PATIENT_ID[match(rownames(gene_coef), x$orig.ident)]
    # gene_coef$TP = x$TIMEPOINT[match(rownames(gene_coef), x$orig.ident)]
    # ggplot(gene_coef, aes(TP, value)) +
    #   geom_point() + geom_line(aes(group=PATIENT)) +
    #   ylab("fitted coefficients")
    # mean(gene_coef$value[gene_coef$TP == "VeryLate"]) - mean(gene_coef$value[gene_coef$TP == "Late"])

    fit.cont <- contrasts.fit(fit, contrasts=mycontr)
    fit.cont <- eBayes(fit.cont, robust=T, trend=F)
    fdr <- p.adjust(fit.cont$p.value[,1], method="BH")

    # Get mean cell type proportions and relative risk for output
    # If no confounding variable included in design matrix
    if(length(mycontr)==2){
      fit.prop <- lmFit(prop, design)
    } else{ # If confounding variables included in design matrix exclude them
      new.des <- design[, mycontr!=0]
      fit.prop <- lmFit(prop, new.des)
    }

    out <- data.frame(
      PropMean = fit.prop$coefficients,
      LogFC = log2(fit.prop$coefficients[, grp1] / fit.prop$coefficients[, grp2]),
      Tstatistic=fit.cont$t[,1],
      P.Value=fit.cont$p.value[,1],
      FDR=fdr
    )
    tt = topTable(fit.cont, n = Inf)
    out$LogFC_cont = tt$logFC[match(rownames(out), rownames(tt))]

    out = out[order(out$P.Value), ]
    out$GROUP = i
    out$CTRST = paste0(grp1, "_vs_", grp2)
    out$CLUSTER = rownames(out)
    colnames(out)[colnames(out) == paste0("PropMean.", grp1)] = "PropMean.Grp1"
    colnames(out)[colnames(out) == paste0("PropMean.", grp2)] = "PropMean.Grp2"

    ct.pr = table(x$orig.ident, x$celltype) > 1
    nbr.smpl.1 = as.character(unique(x[x$TIMEPOINT == grp1, ]$orig.ident))
    nbr.smpl.2 = as.character(unique(x[x$TIMEPOINT == grp2, ]$orig.ident))
    ct.1 = colSums(ct.pr[rownames(ct.pr) %in% nbr.smpl.1, ])
    ct.2 = colSums(ct.pr[rownames(ct.pr) %in% nbr.smpl.2, ])
    out[[paste0("Samples.", "Grp1")]] = length(nbr.smpl.1)
    out[[paste0("Samples.", "Grp2")]] = length(nbr.smpl.2)
    out[[paste0("Present.", "Grp1")]] = ct.1[match(rownames(out), names(ct.1))]
    out[[paste0("Present.", "Grp2")]] = ct.2[match(rownames(out), names(ct.2))]

    l[[i]] = out
  }
  do.call("rbind", l)
}
