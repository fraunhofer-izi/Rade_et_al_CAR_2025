bm.tbl = read.csv("data/genes_table.csv")
bm.tbl = bm.tbl[bm.tbl$gene_biotyp == "protein_coding", ]
bm.ftrs = unique(bm.tbl$external_gene_name)

go.immu = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
pat = "Glycolytic|oxidative_phos|aerobic|atp_metabolic|CD4|CD8|\\_T\\_|immune|T_CELL|cytotox|cytokine|tumor|MHC|ANTIGEN|TNF|LEUKOCYTE|INTERFERON|PROTEOLYSIS|INNATE|APOPTOT|DEATH|KILL|cytolysis|KAPPAB"

go.immu = go.immu[grepl(pat, go.immu$gs_name, ignore.case = T), ]
go.immu.pathways = unique(go.immu$gs_name)
go.immu.ftrs = unique(go.immu$human_gene_symbol)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plot (pre filtering) with cutoffs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_vln_plot_cell = function(
    obj = se.meta,
    .features = "nFeature_RNA",
    .group.by = "orig.ident",
    .split.by = "STUDY",
    plot_title = "Genes Per Cell",
    x_axis_label = NULL,
    y_axis_label = "Features",
    low_cutoff = NULL,
    high_cutoff = NULL
){
  library(ggplot2)
  library(ggthemes)

  if(!is.list(obj)) {
    df = obj@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    df
  }
  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      df = x@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    })
    df = do.call("rbind", l)
    df
  }
  if(class(obj) == "data.frame") {
    df = obj %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    df
  }

  ggplot(data = df, mapping = aes(x = .data[[.group.by]], y = .data[[.features]], groups = .data[[.split.by]])) +
    geom_violin(
      size = .1,
      width = 1,
      scale = "area",
      na.rm = TRUE
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, colour = "#0077BB", size = .2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2))
    ) +
    facet_grid(. ~ .data[[.split.by]], scales = "free", space='free') +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "#BB5566", size = .3) +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
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
# Boxplot: Median values per sample (post filtering)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_median_plot = function(
    obj = se.meta,
    .var = "nFeature_RNA",
    .dotsize = 1,
    .group_by = "TIMEPOINT_WK_MO_BINNED",
    .color_by = "STUDY",
    plot_title = "Median Genes/Cell per Sample",
    y_axis_label = "Median Genes",
    y.max = NA
){

  library(Seurat)
  library(ggplot2)
  library(ggthemes)
  library(scCustomize)
  library(dplyr)

  calc.stat.tbl = function(se.obj) {
    if(.var == "cells_per_sample") {
      stat.tbl = table(se.obj@meta.data[["orig.ident"]]) %>%
        data.frame() %>%
        dplyr::rename(!!"orig.ident" := .data[["Var1"]], Number_of_Cells = .data[["Freq"]])
    } else {
      stat.tbl <- scCustomize::Median_Stats(se.obj, "orig.ident", median_var = .var, default_var = FALSE) %>%
        dplyr::slice(-n()) %>%
        droplevels()
    }
    pheno = se.obj@meta.data[!duplicated(se.obj$orig.ident), ]
    stat.tbl[[.group_by]] = pheno[[.group_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl[[.color_by]] = pheno[[.color_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl
  }

  if(!is.list(obj)) {
    df = calc.stat.tbl(obj)
  }

  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      calc.stat.tbl(x)
    })
    df = do.call("rbind", l)
  }

  if(class(obj) == "data.frame") {
    if(.var == "cells_per_sample") {
      stat.tbl = table(obj[["orig.ident"]]) %>%
        data.frame() %>%
        dplyr::rename(!!"orig.ident" := .data[["Var1"]], Number_of_Cells = .data[["Freq"]])
    } else {
      stat.tbl <- obj %>% group_by(.data[["orig.ident"]]) %>%
        summarise_at(vars(one_of(.var)), median)
      colnames(stat.tbl) <- c("orig.ident", paste0("Median_",  .var))
    }

    pheno = obj[!duplicated(obj), ]
    stat.tbl[[.group_by]] = pheno[[.group_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl[[.color_by]] = pheno[[.color_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    df = stat.tbl
  }

  ggplot(data = df, mapping = aes(x = .data[[.group_by]], y = .data[[colnames(df)[2]]], fill = .data[[.color_by]])) +
    geom_boxplot(fill = "white", outlier.colour = NA, lwd = .3) +
    geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = .dotsize, colour = NA) +
    scale_fill_manual(values = ggthemes::tableau_color_pal("Tableau 10")(10)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
    ) +
    ggtitle(plot_title) +
    ylab(y_axis_label) +
    xlab("") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y.max))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for features (assay)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_features = function(
    .obj = NULL,
    features = NULL,
    .reduc = "tsne",
    pl.points = F,
    pt.size = .3,
    plot.grid = T,
    base.size = 8,
    .x.title = "UMAP 1",
    .y.title = "UMAP 2",
    .leg.title = NULL,
    .title.size = 1,
    .assay = "RNA",
    log.scale.counts = F,
    .quantile.fltr = T,
    .na_cutoff = NULL,
    min.max = NULL,
    legend.wh = c(.3, 4),
    order = T,
    .raster.scattermore = FALSE,
    .raster.scattermore.pixel = c(512,512),
    .raster.scattermore.pointsize = 0,
    .raster = FALSE,
    .raster.dpi = 150,
    .colors = rev(MetBrewer::met.brewer("Hokusai1",n=100))
) {

  DefaultAssay(.obj) = .assay

  library(scales)
  library(cowplot)
  library(scico)

  if(!is.null(.x.title) & !is.null(.y.title)){
    if(grepl("sne", .reduc, ignore.case = T)) {.x.title = "tSNE 1"; .y.title = "tSNE 2"}
    if(grepl("umap", .reduc, ignore.case = T)) {.x.title = "UMAP 1"; .y.title = "UMAP 2"}
  }

  features = features[features %in% rownames(.obj)]
  if(length(features) == 0) { stop("None of the genes are present in the data set")}

  if(log.scale.counts == T){
    exprs.sub = .obj@assays[[.assay]]@counts[features, , drop = F]
    exprs.sub = log10(exprs.sub + 1)
  } else {
    exprs.sub = .obj@assays[[.assay]]@data[features, , drop = F]
  }

  ftr.pl = list()
  for (i in 1:length(features)) {
    ftr = features[i]
    reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
    colnames(reduc) = c("DIM_1", "DIM_2")
    reduc$EXPRS = exprs.sub[ftr, ]

    if(order == T) {
      reduc = reduc[order(reduc$EXPRS, decreasing = F), ] # For ggplot: highest values on the top
    } else {
      set.seed(1234)
      reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
    }

    ftr.exprs = reduc$EXPRS

    if(.quantile.fltr) {
      qu.max = quantile(ftr.exprs[ftr.exprs > 0], .999)
      ftr.exprs[ftr.exprs > qu.max] = qu.max
      reduc$EXPRS = ftr.exprs
    }

    if(!is.null(min.max)) {
      e.min = min.max[1]
      e.max = min.max[2]
      ftr.exprs[ftr.exprs > e.max] = e.max
    } else if(!is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs[ftr.exprs > .na_cutoff])
      e.max = max(ftr.exprs)
    } else if (is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs)
      e.max = max(ftr.exprs)
    }

    if(is.null(names(ftr))) {names(ftr) = ""}

    pl =
      ggplot(reduc, aes(x = DIM_1, y = DIM_2, color = EXPRS))
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl + geom_point(size = pt.size)
      }
    }
    pl = pl +  scale_color_gradientn(
      colors = .colors,
      na.value = "#DDDDDD",
      limits = c(e.min, e.max),
      breaks = pretty_breaks(4)
    )
    pl = pl +  guides(
      color = guide_colorbar(
        title = .leg.title, title.hjust = 0, barwidth = unit(legend.wh[1],'lines'),
        barheight = unit(legend.wh[2], 'lines'),
        ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt,
      )
    ) +
      {if(nchar(names(ftr)) == 0)ggtitle(ftr)} +
      {if(!nchar(names(ftr)) == 0)ggtitle(names(ftr))} +
      mytheme(base_size = base.size) +
      theme(
        aspect.ratio = 1,
        legend.justification = c(0,.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = 0.5, colour = "black", size = rel(.title.size))
      ) +
      xlab(.x.title) + ylab(.y.title)

    if (.raster == T) {
      pl = ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
    }

    ftr.pl[[i]] = pl
  }

  if (plot.grid == T) {
    if(length(ftr.pl) <= 4) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (4 - length(ftr.pl))))}
    if(length(ftr.pl) > 4 && length(ftr.pl) <= 8) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (8 - length(ftr.pl))))}
    plot_grid(plotlist = ftr.pl, ncol = 4, scale = .95, align = "vh")
  } else {
    ftr.pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for phenodata (meta.data)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_pheno = function(
    .obj,
    .target = NULL,
    .reduc = "umap",
    .raster.ggrastr = FALSE,
    .raster.scattermore = FALSE,
    .raster.scattermore.pixel = c(512,512),
    .raster.scattermore.pointsize = 0,
    .raster.dpi = 150,
    .sort = F,
    na_cutoff = NULL,
    .col.palette = "2",
    .col.pal.dicrete = ggthemes::tableau_color_pal()(10),
    .col.scico = "roma",
    .col.scico.d = -1,
    .col.scico.b = 0,
    .col.scico.e = 1,
    quantile.fltr = F,
    pl.points = F,
    pt.size = 1,
    legend.size = 6
) {

  reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
  colnames(reduc) = c("DIM_1", "DIM_2")
  reduc = cbind(.obj@meta.data, reduc)

  if(is.null(.target)) {
    reduc$tmp = "1"
    .target = "tmp"
    cont = F
  } else {
    cont = is.numeric(reduc[[.target]])
  }

  if(grepl("sne", .reduc, ignore.case = T)) {
    .x.title = "tSNE 1"; .y.title = "tSNE 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }
  if(grepl("umap", .reduc, ignore.case = T)) {
    .x.title = "UMAP 1"; .y.title = "UMAP 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }

  if(quantile.fltr) {
    qu.max = quantile(reduc[[.target]][ reduc[[.target]] > 0 ], .999)
    reduc[[.target]][reduc[[.target]] > qu.max] = qu.max
  }

  if(.sort == T) {
    reduc = reduc[order(reduc[[.target]], decreasing = F), ]
  } else {
    set.seed(1234)
    reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
  }

  if(!is.null(na_cutoff)) {
    e.min = min(reduc[[.target]][reduc[[.target]] > na_cutoff])
    e.max = max(reduc[[.target]])
  } else if (is.null(na_cutoff) & cont) {
    e.min = min(reduc[[.target]])
    e.max = max(reduc[[.target]])
  }

  col.cont = list(
    "1" = rev(MetBrewer::met.brewer("Hokusai1",n=100)),
    "2" = scico::scico(30, palette = .col.scico, direction = .col.scico.d, begin = .col.scico.b, end = .col.scico.e)
  )

  pl =
    ggplot(data = reduc, aes(x = DIM_1, y = DIM_2, col = .data[[.target]]))
  # scattermore::geom_scattermore(pointsize = 5, color="black")+
  # scattermore::geom_scattermore(pointsize = 4, color="white")
  if (.raster.scattermore == T) {
    pl = pl + scattermore::geom_scattermore(
      pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
    )
  } else {
    if (pl.points == F) {
      pl = pl + geom_point(shape = ".", alpha = 1)
    } else {
      pl = pl + geom_point(size = pt.size)
    }
  }
  pl = pl + theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = unit(1, "lines"),
    legend.justification = c(0,.5),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1))
  ) +
    xlab(.x.title) + ylab(.y.title)
  if(cont) {
    pl = pl +
      scale_color_gradientn(
        colors = col.cont[[.col.palette]],
        na.value = "#DDDDDD",
        limits = c(e.min, e.max),
        breaks = scales::pretty_breaks(3)
      ) +
      guides(
        color = guide_colorbar(
          title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(6, 'lines')
        )
      )
  } else {
    pl = pl +
      scale_color_manual(values = .col.pal.dicrete, na.value = "#DDDDDD") +
      guides(alpha = 'none') +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size=legend.size, shape = 16, alpha = 1)))

  }

  if (.raster.ggrastr == T) {
    ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
  } else {
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc colored by celltype
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_celltypes = function(
    obj = NULL,
    dim = "wnnUMAP",
    ct.to.pl = "celltype_short_3",
    pt.size = .1,
    leg.size = 2.5,
    leg.ncol = 1,
    col.cc = F,
    base.size = 8,
    ncol1 = 1,
    ncol2 = 1,
    ncol3 = 1,
    raster = F,
    raster.pt.size = 1
) {

  pd = get_metadata(obj)
  pd$DIM1 = pd[[paste0(dim, "_1")]]
  pd$DIM2 = pd[[paste0(dim, "_2")]]

  if(col.cc == T) {
    pd.cc = subset(pd, CellCycle == T)
    pd.cc[[ct.to.pl]] = "Cycling"
    pd = droplevels(pd[!pd$cell %in% pd.cc$cell, ])
  }

  minor.cts = names(table(pd[[ct.to.pl]])[table(pd[[ct.to.pl]]) < 10])
  pd.minor = pd[pd[[ct.to.pl]] %in% minor.cts, ]
  pd = pd[!pd[[ct.to.pl]] %in% minor.cts, ]
  pd = droplevels(pd)

  pd.lym = pd[grepl("T-Cell|B-Cell|NK|gdT|^Plasma", pd[[ct.to.pl]]), ]
  pd.lym[[ct.to.pl]] = factor(
    pd.lym[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.lym[[ct.to.pl]]])
  )

  pd.my = pd[grepl("Mono|Macrophage|other DC|cDC|pDC|Erythrocyte|Platelet", pd[[ct.to.pl]]), ]
  pd.my[[ct.to.pl]] = factor(
    pd.my[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.my[[ct.to.pl]]])
  )

  pd.other = pd[!pd[[ct.to.pl]] %in% as.character(c(unique(pd.lym[[ct.to.pl]]), unique(pd.my[[ct.to.pl]]))), ]
  if(col.cc == T) {
    pd.other = rbind(pd.other, pd.cc)
  }
  pd.other[[ct.to.pl]] = factor(
    pd.other[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.other[[ct.to.pl]]])
  )

  stopifnot(
    length(colnames(obj)) == ( nrow(pd.lym) + nrow(pd.my) + nrow(pd.other) + nrow(pd.minor) )
  )

  pl = ggplot()
  pl = pl + scattermore::geom_scattermore(data = pd, aes(x = DIM1, y = DIM2), pointsize = 6, color="black")+
    scattermore::geom_scattermore(data = pd, aes(x = DIM1, y = DIM2), pointsize = 4, color="white")
  if(raster == T) {
    pl = pl + ggrastr::geom_jitter_rast(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = .5)
    # pl = pl + scattermore::geom_scattermore(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl + guides(colour = guide_legend(
    title = "Other", title.position = "top", ncol = ncol1, order = 3, override.aes = list(shape = 16, size = leg.size)
  )) +
    scale_colour_manual(values = c(ct.col, setNames("#BBBBBB", "green")), na.value = "green") +
    ggnewscale::new_scale_color()
  if(raster == T) {
    pl = pl + ggrastr::geom_jitter_rast(data = pd.my, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = .5)
    # pl = pl + scattermore::geom_scattermore(data = pd.my, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.my, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl +  guides(colour = guide_legend(
    title = "Myeloid", title.position = "top", ncol = ncol2, order = 2, override.aes = list(shape = 16, size = leg.size)
  )) +
    scale_colour_manual(values = c(ct.col, setNames("#BBBBBB", "Cycling")), na.value = "green") +
    ggnewscale::new_scale_color()
  if(raster == T) {
    pl = pl + ggrastr::geom_jitter_rast(data = pd.lym, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".", raster.dpi = 300, scale = .5)
    # pl = pl + scattermore::geom_scattermore(data = pd.lym, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.lym, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl +  guides(colour = guide_legend(
    title = "Lymphoid", title.position = "top", ncol = ncol3, order = 1, override.aes = list(shape = 16, size = leg.size)
  )) +
    scale_colour_manual(values = c(ct.col, setNames("#FFB92D", "Cycling")), na.value = "green") +
    mytheme(base_size = base.size) +
    theme(
      legend.text = element_text(size = rel(.9)),
      legend.spacing.y = unit(1, 'mm'),
      legend.key.size = unit(3, "mm"),
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(size = 10, hjust = 1.1, face = "plain")
    ) +
    xlab("UMAP 1") + ylab("UMAP 2")

  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Box plots for cell fractions per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sc_ct_sample_fraction = function(
    inpMeta,
    label,
    group.facet = "None",
    group.color = "None",
    group.color.pal = c("#6699CC", "#997700", "#994455", "#004488", "black", "#DDAA33", "#BBBBBB"),
    dot.color = "None",
    group.psz = .5,
    base.size = 8,
    scales = "fixed",
    order.by.ave.prop = T,
    nbr.cell.cut = 100,
    filter.ct = NULL
){

  library(ggh4x)

  inpMeta = data.frame(inpMeta)

  if(dot.color == "None"){
    dot.color = NULL
  }

  cell.ident = inpMeta
  cell.ident$celltype = cell.ident[[label]]
  cell.ident = droplevels(cell.ident)
  cell.ident[[group.color]] = factor(cell.ident[[group.color]])

  cell.ident$groups_ct = paste0(
    cell.ident[[group.facet]], "_xx_", cell.ident[[group.color]], "_yy_", cell.ident$celltype
  )

  # Extract total cell number per group
  total_cells = data.frame(table(cell.ident$groups_ct)) %>%
    dplyr::rename(groups_ct = Var1, groups_ct_number = Freq)

  total_cells$celltype = gsub(".+_yy_", "", total_cells$groups_ct)
  total_cells$groups = gsub("_yy_.+", "", total_cells$groups_ct)

  l = list()
  for (i in unique(total_cells$groups)) {

    target = as.character(total_cells[total_cells$groups %in% i, ]$groups_ct)

    # Calculate and extract percents of cells per cluster per group
    g = cell.ident[cell.ident$groups_ct %in% target, ]
    perc_per_groups_ct <- prop.table(x = table(g$groups_ct, g$orig.ident), margin = 2) * 100

    # Remove sample not included in this group
    empty_columns = colSums(is.na(perc_per_groups_ct) | perc_per_groups_ct == "") == nrow(perc_per_groups_ct)
    perc_per_groups_ct = perc_per_groups_ct[, !empty_columns, drop = F]

    perc_per_groups_ct = data.frame(perc_per_groups_ct) %>%
      dplyr::rename(groups_ct = Var1, orig.ident = Var2, perc_sample = Freq)
    l[[i]] = perc_per_groups_ct
  }

  perc_per_groups_ct = do.call("rbind", l)
  rownames(perc_per_groups_ct)  = NULL

  check = perc_per_groups_ct %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(Frequency = sum(perc_sample)) %>%
    data.frame()
  stopifnot(!any(check$Frequency != "100"))

  perc_per_groups_ct = merge(perc_per_groups_ct, total_cells, by = c("groups_ct"))

  # Add groups for ggplot
  if(!is.null(dot.color)) {
    perc_per_groups_ct$dot.color = cell.ident[[dot.color]][match(perc_per_groups_ct$orig.ident, cell.ident$orig.ident)]
  } else {
    perc_per_groups_ct$dot.color = "None"
  }

  perc_per_groups_ct$group_facet = gsub("_xx_.+", "", perc_per_groups_ct$groups)
  perc_per_groups_ct$group_facet = factor(
    perc_per_groups_ct$group_facet, levels = levels(cell.ident[[group.facet]])
  )

  perc_per_groups_ct$group_color = gsub(".*_xx_", "", perc_per_groups_ct$groups)
  perc_per_groups_ct$group_color = factor(
    perc_per_groups_ct$group_color, levels = levels(cell.ident[[group.color]])
  )

  # Remove NA
  if(any(is.na(perc_per_groups_ct$group_color))) {
    print("NA is present in the data")
    perc_per_groups_ct = perc_per_groups_ct[!is.na(perc_per_groups_ct$group_color), ]
  }
  perc_per_groups_ct = droplevels(perc_per_groups_ct)

  # Remove celltypes with less than n cells in two groups
  df.sum = perc_per_groups_ct %>%
    group_by(celltype, groups, group_facet) %>%
    dplyr::slice(1)
  df.sum = df.sum %>% dplyr::group_by(celltype, group_facet) %>%
    dplyr::summarise(n = sum(groups_ct_number)) %>%
    dplyr::filter(n < nbr.cell.cut) %>%
    data.frame()
  if(nrow(df.sum) != 0) {
    df.sum$REMOVE = paste0(df.sum$celltype, "_", df.sum$group_facet)
  } else {
    df.sum$REMOVE = NULL
  }
  perc_per_groups_ct$TMP =  paste0(perc_per_groups_ct$celltype, "_", perc_per_groups_ct$group_facet)
  perc_per_groups_ct = perc_per_groups_ct[!perc_per_groups_ct$TMP %in% unique(df.sum$REMOVE), ]

  # Remove celltypes, where only one sample hast cells
  perc_per_groups_ct$TMP = perc_per_groups_ct$perc_sample > 0
  df.sum = perc_per_groups_ct %>% dplyr::group_by(celltype, group_facet) %>%
    dplyr::summarise(n = sum(TMP)) %>%
    dplyr::filter(n == 1) %>%
    data.frame()
  if(nrow(df.sum) != 0) {
    df.sum$REMOVE = paste0(df.sum$celltype, "_", df.sum$group_facet)
  } else {
    df.sum$REMOVE = NULL
  }
  perc_per_groups_ct$TMP =  paste0(perc_per_groups_ct$celltype, "_", perc_per_groups_ct$group_facet)
  perc_per_groups_ct = perc_per_groups_ct[!perc_per_groups_ct$TMP %in% unique(df.sum$REMOVE), ]
  perc_per_groups_ct$TMP = NULL

  perc_per_groups_ct$perc_sample = perc_per_groups_ct$perc_sample / 100

  if(!is.null(filter.ct)) {
    perc_per_groups_ct = perc_per_groups_ct[grepl(filter.ct, perc_per_groups_ct$celltype), ]
    perc_per_groups_ct = droplevels(perc_per_groups_ct)
  }

  if(order.by.ave.prop == T){
    lvls = perc_per_groups_ct %>%
      dplyr::group_by(celltype, group_facet) %>%
      dplyr::summarise(mean_prop = mean(perc_sample)) %>%
      dplyr::arrange(-mean_prop) %>% dplyr::select(celltype) %>% unlist
  } else {
    lvls = naturalsort(unique(as.character(perc_per_groups_ct$celltype)))
  }
  perc_per_groups_ct$celltype = factor(perc_per_groups_ct$celltype, levels = unique(lvls))

  if(length(as.character(unique(perc_per_groups_ct$group_color))) == 2) {

    set.seed(123)
    pl = ggplot(perc_per_groups_ct, aes(celltype, perc_sample))
    pl = pl + geom_boxplot(
      aes(fill = group_color), position = position_dodge2(.85, padding = .2, preserve = "single"),
      outlier.shape = NA, fatten = 1.5, linewidth = .2
    )
    if(!is.null(dot.color)) {
      pl = pl + geom_point(
        aes(col = dot.color, group = group_color), size=group.psz,
        position = position_jitterdodge(jitter.width = .2)
      ) + scale_color_manual(values = c("#DDAA33", "#BB5566", "#004488", colors_use.10))
      pl = pl + guides(color = guide_legend(title = dot.color, nrow = 1, override.aes = list(size = 2.5)))
    } else {
      pl = pl + geom_point(
        aes(col = dot.color, group = group_color), size=group.psz,
        position = position_jitterdodge(jitter.width = .2)
      ) +
        scale_color_manual(values = c("#555555")) +
        guides(color = "none")
    }
    pl = pl + ylab("Cell type proportion per sample") +
      xlab(NULL) +
      labs(fill = group.color) +
      mytheme_grid(base_size = base.size) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = rel(1), face = "plain"),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "#1A242F")
      ) +
      scale_x_discrete(drop=FALSE) +
      scale_fill_manual(values = group.color.pal, drop = F, na.value = "#BBBBBB")


    # Color
    pl = pl + guides(fill = guide_legend(nrow = 1))

    # Facetting
    if(group.facet != "None") {
      pl =
        pl + facet_wrap2(~ group_facet, ncol = 1, drop = F, scales = scales) +
        annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=.5)
    }
  } else {

    set.seed(123)
    pl = ggplot(perc_per_groups_ct, aes(group_color, perc_sample))
    pl = pl + geom_boxplot(
      aes(fill = celltype), fill = "white", position = position_dodge2(.85, padding = .2, preserve = "single"),
      outlier.shape = NA, fatten = 1.5, linewidth = .2
    )
    if(!is.null(dot.color)) {
      pl = pl + geom_point(
        aes(col = dot.color, group = celltype), size=group.psz,
        position = position_jitterdodge(jitter.width = .2)
      ) + scale_color_manual(values = colors_use.10)
      pl = pl + guides(color = guide_legend(title = dot.color, nrow = 1, override.aes = list(size = 2.5)))
    } else {
      pl = pl + geom_point(
        aes(col = dot.color, group = celltype), size=group.psz,
        position = position_jitterdodge(jitter.width = .2)
      ) +
        scale_color_manual(values = c("#555555")) +
        guides(color = "none")
    }
    pl = pl + ylab("Cell type proportion per sample") +
      xlab(NULL) +
      labs(fill = group.color) +
      mytheme_grid(base_size = base.size) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = rel(1), face = "plain"),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "#1A242F")
      ) +
      scale_x_discrete(drop=FALSE) +
      facet_wrap(~ celltype, nrow = 2)

    comparisons = combn(as.character(unique(perc_per_groups_ct$group_color)), 2, simplify = F)
    pl = pl + ggpubr::stat_compare_means(comparisons = comparisons, size = 3)
  }

  return(pl)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Celltype composition per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
comp_celltypes = function(
    p.data = NULL,
    t = "celltype_short_2",
    group.facet = "None",
    group.color = "RESPONSE_CONSENSUS_2",
    dot.color = "STUDY",
    base.size = 8,
    title.size = 10,
    group.psz = .35,
    nbr.cell.cut = 50
) {

  p.data = subset(p.data, CAR_BY_EXPRS == F)

  nbr.grps = length(as.character(unique(p.data[[group.color]])))
  if (nbr.grps == 2) {
    coord.y.max = 2
  } else {
    coord.y.max = 20
  }

  comp.pl =
    sc_ct_sample_fraction(
      inpMeta = p.data,
      label = t,
      group.facet = group.facet,
      group.color = group.color,
      dot.color = dot.color,
      nbr.cell.cut = nbr.cell.cut,
      group.psz = group.psz,
      base.size = base.size
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10),
      breaks=c(0, 0.001, 0.01, 0.1, 1)
    ) +
    coord_cartesian(ylim=c(0, coord.y.max))

  if (nbr.grps == 2) {
    # Speckle
    res.speckle = speckle::propeller(
      clusters = p.data[[t]], sample = p.data$orig.ident,
      group = p.data[[group.color]], transform = "asin"
    )

    dat = dplyr::distinct(comp.pl$data, celltype)
    res.speckle = res.speckle[rownames(res.speckle) %in% dat$celltype, ]
    res.speckle = res.speckle[dat$celltype, ]
    res.speckle$yloc = max(comp.pl$data$perc_sample) + .75
    res.speckle = add_signif(res.speckle,"P.Value", "pval_star", pval.relax = T)
    colnames(res.speckle)[1] = "celltype"
    res.speckle = droplevels(res.speckle)

    comp.pl +
      geom_text(data = res.speckle, aes(y = yloc, label = pval_star), size = 2, position = position_dodge(width = .75)) +
      guides(fill = guide_legend(title = NULL,  order = 1)) +
      guides(color = guide_legend(title = NULL, override.aes = list(shape = 16, size = 2.5))) +
      theme(
        plot.title = element_text(size = title.size, hjust = 0.5, face = "plain"),
        legend.position = "right"
      )
    # scale_color_manual(values = c("#994455", "black"))
  } else {
    comp.pl
  }

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA scatterplot with log fold change (y-axis)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_plot = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 5,
    base.size = 8,
    text.repel.size = 2.5,
    text.de.nbr.size = 2.5,
    text.cl.size = 2.5,
    ylim.extend.up = .5,
    ylim.extend.dn = .7,
    box.padding = 0.25,
    dge.nbr.up = 0,
    cluster = "celltype",
    subset.tops = NULL,
    axis.max = NULL,
    axis.min = NULL,
    legend.margin.t = -10,
    cluster.label.y = .15
){

  dgea.res$cluster = dgea.res[[cluster]]
  dgea.res.sign$cluster = dgea.res.sign[[cluster]]

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]

  cluster.de.summary = dgea.res.sign %>%
    dplyr::mutate(IsUp = avg_log2FC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(Up = sum(IsUp), Down = sum(!IsUp)) %>%
    dplyr::mutate(Down = -Down) %>%
    tidyr::gather(key = "Direction", value = "Count", -cluster)

  cluster.de.summary = cluster.de.summary[order(cluster.de.summary$cluster), ]
  cluster.de.summary$cluster = factor(cluster.de.summary$cluster, levels = unique(cluster.de.summary$cluster))

  cl.label = paste0(
    "Up: ",
    subset(cluster.de.summary, Direction == "Up")$Count,
    "\nDown: ",
    abs(subset(cluster.de.summary, Direction == "Down")$Count)
  )

  dgea.res.sign.top = dgea.res.sign
  if(!is.null(subset.tops)){
    dgea.res.sign.top = dgea.res.sign.top[dgea.res.sign.top$feature %in% subset.tops, ]
  }

  top5_up = dgea.res.sign.top %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = nbr.tops, wt = avg_log2FC)
  top5_up$GENE_CL = paste0(top5_up$feature, "_", top5_up$cluster)

  top5_down = dgea.res.sign.top %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = -nbr.tops, wt = avg_log2FC)
  top5_down$GENE_CL = paste0(top5_down$feature, "_", top5_down$cluster)

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% mutate(label_up = ifelse(GENE_CL %in% top5_up$GENE_CL, feature, ""))
  dgea.res = dgea.res %>% mutate(label_down = ifelse(GENE_CL %in% top5_down$GENE_CL, feature, ""))
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = ifelse(dgea.res$SIGNIFICANT == T, "yes", "no")

  if(is.null(axis.max)){
    axis.max = max(abs(dgea.res$avg_log2FC)) + .75
  } else {
    axis.max = axis.max
  }
  if(is.null(axis.min)){
    axis.min = -axis.max
  } else {
    axis.min = axis.min
  }

  no.cl = length(unique(dgea.res$cluster))
  top1.up = dgea.res.sign %>% group_by(cluster) %>%
    top_n(n = 1, wt = avg_log2FC) %>%
    data.frame()
  rownames(top1.up) = top1.up$cluster
  top1.up = top1.up[levels(cluster.de.summary$cluster), ]
  top1.down = dgea.res.sign %>% group_by(cluster) %>%
    top_n(n = -1, wt = avg_log2FC) %>%
    data.frame()
  rownames(top1.down) = top1.down$cluster
  top1.down = top1.down[levels(cluster.de.summary$cluster), ]

  stopifnot(identical(top1.up$cluster, top1.down$cluster))

  pos = position_jitter(width = 0.35, seed = 1234)

  data = data.frame(
    xmin = seq(.6, no.cl),
    xmax = seq(.2, no.cl) + 1.2,
    ymin = top1.down$avg_log2FC - .1,
    ymax = top1.up$avg_log2FC + .1,
    group = top1.up$cluster
  )

  set.seed(1234)
  ggplot() +
    ylim(c(axis.min, axis.max)) +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_jitter(
      data = subset(dgea.res, SIGNIFICANT == "no" & abs(avg_log2FC) > .15),
      aes(x= cluster, y= avg_log2FC, color = SIGNIFICANT), position = pos, size = .3
    ) +
    geom_jitter(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      aes(x= cluster, y= avg_log2FC, color = SIGNIFICANT), position = pos, size = .3
    ) +
    geom_rect(
      data = data,
      aes(x = group, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#DDDDDD", alpha = .3
    ) +
    geom_rect(
      data = data.frame(xmin = seq(0.5,no.cl), xmax = seq(0.5,no.cl) + 1, ymin = -cluster.label.y, ymax = cluster.label.y),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", color = "#BBBBBB", lwd = .2
    ) +
    geom_text_repel(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      mapping = aes(x = cluster, avg_log2FC, label = label_up),
      position = pos,
      size = text.repel.size,
      max.overlaps = 100,
      force = 4,
      min.segment.length = 0,
      box.padding = box.padding,
      segment.size = .1,
      direction = "y",
      ylim = c(.1, max(dgea.res$avg_log2FC) + ylim.extend.up)
    ) +
    geom_text_repel(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      mapping = aes(x = cluster, avg_log2FC, label = label_down),
      position = pos,
      size = text.repel.size,
      max.overlaps = 100,
      force = 4,
      min.segment.length = 0,
      box.padding = box.padding,
      segment.size = .1,
      direction = "y",
      ylim = c(-.1, min(dgea.res$avg_log2FC) - ylim.extend.dn)
    ) +
    mytheme(base_size = base.size) +
    theme(
      legend.position="bottom",
      axis.line.y = element_line(colour = "black"),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.margin = margin(t = legend.margin.t),
      legend.title = element_text(margin = margin(r = 3)),
      legend.text = element_text(margin = margin(l = 1))
      # legend.spacing.x = unit(0, 'mm')
    ) +
    scale_color_manual(values = c("yes" = "#6699CC", "no" = "#BBBBBB")) +
    labs(x = NULL, y = "Log2 fold change", colour = paste0("FDR <0.05"))  +
    annotate("text", x = 1:(no.cl), y = 0, label = data$group, color = "black", size = text.de.nbr.size) +
    annotate("text", x = 1:(no.cl), y = axis.max + dge.nbr.up, label = cl.label, colour = "#6699CC", size = text.cl.size) +
    guides(colour = guide_legend(override.aes = list(size=2)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot: TOP DE genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
de_tops_bubble = function(
    df,
    pl.title = NULL,
    nbr.tops = 5,
    w.o.dir = F,
    sort.by.p = T
){

  if(w.o.dir == T) {
    if(sort.by.p  == T){
      tops = df %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
        slice_head(n = nbr.tops)
    } else {
      tops = df %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(abs(avg_log2FC)) %>%
        slice_head(n = nbr.tops)
    }
  } else {
    if(sort.by.p  == T){
      tops.up =
        df %>%
        dplyr::filter(avg_log2FC > 0) %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
        slice_head(n = nbr.tops)
      tops.dwn = df %>%
        dplyr::filter(avg_log2FC < 0) %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(p_val_adj, avg_log2FC) %>%
        slice_head(n = nbr.tops)
    } else {
      tops.up =
        df %>%
        dplyr::filter(avg_log2FC > 0) %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        slice_head(n = nbr.tops)
      tops.dwn = df %>%
        dplyr::filter(avg_log2FC < 0) %>%
        dplyr::group_by(celltype) %>%
        dplyr::arrange(desc(avg_log2FC)) %>%
        slice_head(n = nbr.tops)
    }
    tops = rbind(tops.dwn, tops.up)
  }

  df.pl = df[df$feature %in% tops$feature, ]
  df.pl$feature = factor(df.pl$feature, levels = unique(tops$feature))

  thres = quantile(abs(df.pl$avg_log2FC), .99)
  df.pl$avg_log2FC[df.pl$avg_log2FC > thres & !is.infinite(df.pl$avg_log2FC)] = thres
  df.pl$avg_log2FC[df.pl$avg_log2FC < -thres & !is.infinite(df.pl$avg_log2FC)] = -thres

  ggplot(df.pl, aes(feature, celltype, fill = avg_log2FC, size = -log10(p_val_adj))) +
    geom_point(colour="black", pch=22, stroke = .2) +
    scale_size(range = c(1, 4)) +
    mytheme() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1, size = rel(.9)),
      legend.position = "bottom",
      legend.title = element_text(margin = margin(r = 2, unit = "pt"), size = rel(1)),
      legend.text = element_text(margin = margin(r = 2, t = 1, unit = "pt"), size = rel(1)),
      legend.margin = margin(t=-4, l = 5),
      plot.title = element_text(hjust = 0.5, face = "plain")
    ) +
    scale_fill_scico(
      palette = "vik", midpoint = 0, begin =  0, end = 1, direction = 1,
      limits = c(-thres, thres)
    ) +
    guides(
      fill = guide_colorbar(
        title = "logFC", order = 1,
        title.hjust = 0, title.vjust = 1.1, barwidth = unit(4, 'lines'),
        barheight = unit(.35, 'lines'), ticks.linewidth = .5/.pt
      ),
      size = guide_legend(title = "-log10(FDR)")
    ) + xlab(NULL) + ylab(NULL) +
    ggtitle(pl.title)
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_volcano = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 7,
    cl.label = NULL, # e.g.  cl.label = setNames(c("CD4", "CD8"), c("CD4 T-Cell", "CD8 T-Cell"))
    cluster.col = "cluster",
    sort.by.p = T,
    pval.column = "p_val",
    logFC.column = "avg_log2FC",
    facet.scales = "free_y",
    facet.ncol = 1,
    nudge_x = 2,
    x.axis.sym = F,
    x.axis.ext = 0,
    geom.hline = log2(1.5),
    box.padding = 0.3,
    label.padding = .12,
    label.size = 2,
    leg.title = "FDR <0.05",
    cut.y.thresh = NULL
){

  dgea.res$cluster = dgea.res[[cluster.col]]
  dgea.res.sign$cluster = dgea.res.sign[[cluster.col]]

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]

  if(sort.by.p == T) {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[pval.column]], -abs(.data[[logFC.column]])) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  } else {
    tops_up = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(.data[[logFC.column]], n = nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(.data[[logFC.column]] < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_min(.data[[logFC.column]], n = nbr.tops)
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  }

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% mutate(label_up = ifelse(GENE_CL %in% tops_up$GENE_CL, feature, ""))
  dgea.res = dgea.res %>% mutate(label_down = ifelse(GENE_CL %in% tops_down$GENE_CL, feature, ""))
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = factor(dgea.res$SIGNIFICANT, levels = c(TRUE, FALSE))

  axis.max = max(abs(dgea.res[[logFC.column]]))
  dgea.res = dgea.res %>%
    mutate(p_val = ifelse(p_val == 0, 1e-300, dgea.res[[pval.column]]))
  if(!is.null(cut.y.thresh)){
    dgea.res$p_val[dgea.res$p_val < cut.y.thresh] = cut.y.thresh
  }

  if(is.null(cl.label)) {
    # cl.label = setNames(unique(dgea.res$cluster), unique(dgea.res$cluster))
    tbl = table(
      dgea.res.sign$cluster,
      (dgea.res.sign$significant & dgea.res.sign$avg_log2FC > 0)
    )
    cl.label = setNames(
      paste0(rownames(tbl), ": ", "up: ", tbl[, 2], " | down: ", tbl[, 1]),
      rownames(tbl)
    )
  }

  g <- make_gradient(
    deg = 180, n = 500,
    cols = scico::scico(
      9, palette = 'vik', begin = .3, end = .7, direction = -1,
    )
  )
  set.seed(42)

  pl =
    ggplot(dgea.res, aes(x = .data[[logFC.column]], y = -log10(p_val))) +
    annotation_custom(
      grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT), size = .5) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT), size = .5) +
    # scattermore::geom_scattermore(data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT), pointsize = 3) +
    # scattermore::geom_scattermore(data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT), pointsize = 3) +
    theme(
      panel.spacing = unit(1.5, "lines"),
      panel.grid.minor = element_blank(),
      legend.position="bottom",
      strip.text = element_text(size = rel(1), face = "plain")
    ) +
    facet_wrap(
      ~ cluster, scales = facet.scales, labeller = labeller(cluster = cl.label),
      ncol = facet.ncol,
    ) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_up != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_up != "")$label_up,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      nudge_x = nudge_x,
      nudge_y = -2,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T & label_down != ""),
      label = subset(dgea.res, SIGNIFICANT == T & label_down != "")$label_down,
      segment.colour = "black",
      size = label.size,
      direction = "y",
      hjust = .5,
      # xlim = c(NA, 0),
      nudge_x = -nudge_x,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_vline(xintercept = -geom.hline, linetype = "dashed", linewidth = .2) +
    geom_vline(xintercept = geom.hline, linetype = "dashed", linewidth = .2) +
    scale_color_manual(values = c("TRUE" = "#555555", "FALSE" = "#BBBBBB")) +
    labs(y = "-Log10(p-value)", x = "Log2 fold change", colour = leg.title)  +
    guides(colour = guide_legend(override.aes = list(size=3)))
  if(x.axis.sym == T) {
    pl = pl + xlim(-axis.max - x.axis.ext, axis.max + x.axis.ext)
  }
  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ligand-Receptor analysis
# DB from LIANA. Parsed with iTALK
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ligand_receptor_analysis = function(
    dgea.res.sign = NULL,
    same.direction = F,
    self.interaction = F,
    font.size = 8,
    exportTable = F,
    celltype = "celltype"
) {

  dgea.res.sign$celltype = dgea.res.sign[["celltype"]]

  if (any(!"iTALK" %in% installed.packages())) {
    # Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
  }
  library(iTALK)

  if (any(!"liana" %in% installed.packages())) {
    # Sys.unsetenv("GITHUB_PAT")
    remotes::install_github('saezlab/liana')
  }
  library(liana)

  consensus_omni <- liana::select_resource("Consensus")[[1]]

  liana.db = consensus_omni %>%
    dplyr::select(
      Ligand.ApprovedSymbol = source_genesymbol,
      Receptor.ApprovedSymbol = target_genesymbol,
      Classification = category_intercell_source
    ) %>%
    data.frame()

  getLRI = function(
    res_post_splitted,
    comm_list = c('growth factor','other','cytokine','checkpoint'),
    custom.db = NULL
  ){

    comparison_grid = expand.grid(
      seq_along(res_post_splitted), seq_along(res_post_splitted)
    )
    res_post = NULL
    for (i in seq(nrow(comparison_grid))){
      for(comm_type in comm_list){
        res_cat = FindLR(
          res_post_splitted[[comparison_grid$Var1[i]]],
          res_post_splitted[[comparison_grid$Var2[i]]],
          datatype='DEG',comm_type=comm_type, database = custom.db
        )

        res_post<-rbind(res_post,res_cat)
      }
    }
    return(res_post)
  }

  if (same.direction == T) {
    res_up_splitted = dgea.res.sign %>%
      dplyr::select(gene = feature, cell_type = celltype, logFC = avg_log2FC, p.value = p_val, q.value = p_val_adj) %>%
      dplyr::filter(logFC > 0)
    res_up_splitted = split(res_up_splitted, res_up_splitted$cell_type)

    res_dn_splitted = dgea.res.sign %>%
      dplyr::select(gene = feature, cell_type = celltype, logFC = avg_log2FC, p.value = p_val, q.value = p_val_adj) %>%
      dplyr::filter(logFC < 0)
    res_dn_splitted = split(res_dn_splitted, res_dn_splitted$cell_type)

    # pre_up_LRI = getLRI(res_up_splitted)
    # pre_dn_LRI = getLRI(res_dn_splitted)
    pre_dn_LRI = getLRI(res_dn_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
    pre_up_LRI = getLRI(res_up_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
    pre_LRI = bind_rows(pre_up_LRI, pre_dn_LRI) %>% distinct
  } else {
    res_splitted = dgea.res.sign %>%
      dplyr::select(gene = feature, cell_type = celltype, logFC = avg_log2FC, p.value = p_val, q.value = p_val_adj)
    res_splitted = split(res_splitted, res_splitted$cell_type)
    pre_LRI = getLRI(res_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
  }

  lr.df = pre_LRI

  if(nrow(lr.df) == 0) {
    print("No interactions found")
  }

  lr.df$LIGAND_ID = paste0(lr.df$ligand, "_", lr.df$cell_from)
  lr.df$RECEPTOR_ID = paste0(lr.df$receptor, "_", lr.df$cell_to)
  # lr.df$LIGAND_AVE_EXPRS = dgea.res.sign$avgExpr[match(lr.df$LIGAND_ID, rownames(dgea.res.sign))]
  # lr.df$RECEPTOR_AVE_EXPRS = dgea.res.sign$avgExpr[match(lr.df$RECEPTOR_ID, rownames(dgea.res.sign))]
  # lr.df$AVE_EXPRS = rowMeans(cbind(abs(lr.df$LIGAND_AVE_EXPRS), abs(lr.df$RECEPTOR_AVE_EXPRS)))
  lr.df$AVE_LFC = rowMeans(cbind(abs(lr.df$cell_from_logFC), abs(lr.df$cell_to_logFC)))

  lr.df$LR_PAIR = paste0(lr.df$ligand, " -> ", lr.df$receptor)
  lr.df = lr.df %>% mutate(
    LR_DIR = case_when(
      cell_from_logFC > 0 & cell_to_logFC > 0 ~ "L:up | R:up",
      cell_from_logFC < 0 & cell_to_logFC < 0 ~ "L:dn | R:dn",
      cell_from_logFC > 0 & cell_to_logFC < 0 ~ "L:up | R:dn",
      cell_from_logFC < 0 & cell_to_logFC > 0 ~ "L:dn | R:up"
    )
  )
  lr.df$LR_DIR = factor(lr.df$LR_DIR, levels = c("L:up | R:up", "L:dn | R:dn", "L:up | R:dn", "L:dn | R:up"))

  if(self.interaction == F){
    lr.df = lr.df[!lr.df$cell_from == lr.df$cell_to, ]
  }

  if(nrow(lr.df) == 0){
    print("no gene found")
    stopifnot(nrow(lr.df) > 1)
  }

  max.lfc = max(c(abs(lr.df$cell_from_logFC), abs(lr.df$cell_to_logFC)))
  pl = ggplot(lr.df, aes(x = cell_to, y = LR_PAIR, fill = LR_DIR, size = AVE_LFC)) +
    geom_point(colour="black", pch=21, stroke = .3) +
    scale_size(range = c(1, 4)) +
    # facet_wrap( ~ cell_from, nrow = 1, scales = "free_x") +
    facet_grid( ~ cell_from, scales = "free_x", space = "free_x") +
    guides(color = guide_legend(
      title = "LFC Direction", ncol = 1, order = 2, override.aes = list(shape = 16, size = 3.5)
    )) +
    scale_fill_manual(values = c(
      "L:up | R:up" = "#BB5566",
      "L:dn | R:dn" = "#4477AA",
      "L:up | R:dn" = "#555555",
      "L:dn | R:up" = "#BBBBBB"
    )) +
    mytheme_grid(base_size = font.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1.25)),
      axis.title.x = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1.25)),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    ) +
    xlab("Target") + ylab("Ligand -> Receptor") +
    labs(size = "Mean LFC") +
    ggtitle("Source")

  if(exportTable == T){
    lr.df
  } else{
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot for cell identity markers
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mrkr_bubble = function(
    dgea.res = NULL,
    se.obj = NULL,
    ref.group = "celltype",
    effect.size = "avg_log2FC",
    positive.effect.size = T,
    order.by.effect.size = T,
    nbr.tops = 5,
    quantile.cut = .99,
    features = NULL,
    features.rm = NULL,
    features.grep.rm = NULL,
    font.size = 10,
    aspectRatio = 2,
    pt.size = 1,
    rm.rp = T,
    range = c(1, 4),
    export.table = F,
    rename.ftrs = F,
    rename.ftrs.col = "feature_name",
    add.ftr = NULL,
    barwidth = unit(.5, 'lines'),
    barheight = unit(4, 'lines'),
    legend.pos = "bottom"
) {

  library(scico)
  library(circlize)
  library(ComplexHeatmap)

  DefaultAssay(se.obj) = "RNA"

  dgea.res.all = dgea.res
  dgea.res = dgea.res[dgea.res$significant == T, ]

  if (!is.null(features.rm)) {
    dgea.res = dgea.res[!dgea.res$gene %in% features.rm, ]
  }

  if (!is.null(features.grep.rm)) {
    dgea.res = dgea.res[!grepl(features.grep.rm, dgea.res$gene), ]
  }

  if (rm.rp == T) {
    dgea.res = dgea.res[!grepl("^RPL|^RPS", dgea.res$gene), ]
  }

  if(is.null(features)) {

    dge.res = dgea.res[!is.na(dgea.res[[effect.size]]), ]
    if (positive.effect.size == T) {
      dge.res = dge.res[dge.res[[effect.size]] > 0, ]
    }

    if (order.by.effect.size == T) {
      dge.res = dge.res %>% dplyr::arrange(cluster, desc(abs(!!as.name(effect.size))))
    } else {
      dge.res = dge.res %>% dplyr::arrange(cluster, p_val_adj, desc(abs(!!as.name(effect.size))))
    }

    ct.keep = unique(as.character(dge.res[["cluster"]]))
    se.obj = se.obj[, se.obj@meta.data[[ref.group]] %in% ct.keep]
    se.obj@meta.data = droplevels(se.obj@meta.data)

    tops = dge.res %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = nbr.tops) %>% data.frame()

    if(!is.null(add.ftr)){
      tops = rbind(
        tops,
        dgea.res.all[dgea.res.all$gene == add.ftr[1] & dgea.res.all$cluster == add.ftr[2], ]
      )
    }

    lvls = levels(se.obj@meta.data[[ref.group]])
    lvls.order = setNames(seq(1, length(lvls)), lvls)
    tops$order = lvls.order[match(tops$cluster, names(lvls.order))]
    tops = tops[naturalsort::naturalorder(tops$order), ]
    export.tops = tops
    tops = tops$gene

  } else {
    tops = features
  }

  se.obj = se.obj[unique(tops), ]
  mat = AverageExpression(se.obj, group.by = ref.group, assays = "RNA", slot = "data")[[1]]
  mat = t(scale(t(mat)))
  df = reshape2::melt(mat)
  colnames(df) = c("GENE", "CT", "AVE")
  # df$CT = gsub("^g", "", df$CT)
  # df$CT = gsub("-", "_", df$CT)

  v.max = quantile(df$AVE, quantile.cut)
  df$AVE[df$AVE > v.max & !is.infinite(df$AVE)] = v.max
  df$AVE[df$AVE < -v.max & !is.infinite(df$AVE)] = -v.max

  df$AVE[abs(df$AVE) > v.max] = v.max

  df$PERC = dgea.res.all$pct.1[match(paste0(df$GENE, "_", df$CT), paste0(dgea.res.all$gene, "_", dgea.res.all$cluster))]
  df$PERC = df$PERC * 100
  df$PERC[df$PERC < 0.01] = NA
  if(!is.null(levels(se.obj@meta.data[[ref.group]]))){
    df$CT = factor(df$CT, levels = levels(se.obj@meta.data[[ref.group]]))
  }

  if(rename.ftrs == T) {
    ftrs = se.obj@assays$RNA@meta.features
    lvls = as.character(ftrs[[rename.ftrs.col]])[match(levels(df$GENE), rownames(ftrs))]
    df$GENE = ftrs[[rename.ftrs.col]][match(df$GENE, rownames(ftrs))]
    df$GENE = factor(df$GENE, levels = lvls)
  }

  pl =
    ggplot(df, aes(x = CT,y = GENE, size = PERC)) +
    geom_point(aes(fill = AVE), shape = 21, stroke = .3) +
    scale_size(range = range) +
    scale_fill_scico(
      palette = "vik", midpoint = 0, limits = c(-v.max,v.max),
      breaks = pretty_breaks(n = 3)
    ) +
    mytheme(base_size = font.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(.5, "lines"),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
      axis.text.y = element_text(size = rel(1)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = legend.pos,
      aspect.ratio = aspectRatio
    ) +
    guides(
      size = guide_legend(title = "Percent\nExpressed"),
      fill = guide_colorbar(
        title = "Scaled\nAverage\nExpression", order = 1,
        title.hjust = 0, barwidth = barwidth, barheight = barheight
      )
    )

  if(export.table == T) {
    export.tops
  } else {
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Heatmap: Result from UCell; ucell_enrich()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
enrich_heatmap_milestones = function(
    se.w = NULL,
    meta.vars = c("SAMPLE", "TopClones_2"),
    .slot = "custom_UCell_score",
    pl.title = NULL,
    max.value = NULL,
    lfc.thres = .1
){

  var.1 = meta.vars[[1]]
  var.2 = meta.vars[[2]]
  se.w@meta.data$VAR_1 = se.w@meta.data[[var.1]]
  se.w@meta.data$VAR_2 = se.w@meta.data[[var.2]]

  DefaultAssay(se.w) = .slot
  df = FetchData(se.w, vars = c(c("VAR_1", "VAR_2"), rownames(se.w)), layer = "data")

  df.ftrs =  df[ , -which(names(df) %in% c("VAR_1", "VAR_2"))]
  df.ftrs[df.ftrs < .1] = 0 # Methods !!!
  df = cbind(df[ , which(names(df) %in% c("VAR_1", "VAR_2"))], df.ftrs)

  df.fltr = reshape2::melt(df, ids = c("VAR_1", "VAR_2"))
  df.fltr = df.fltr %>%
    dplyr::group_by(VAR_1, VAR_2, variable) %>%
    dplyr::mutate(total = n()) %>%
    dplyr::mutate(bool = value > 0) %>%
    dplyr::mutate(bool_sum = sum(bool)) %>%
    dplyr::mutate(frac = bool_sum/total) %>%
    data.frame()
  df.fltr = df.fltr[!duplicated(paste0(df.fltr$VAR_1, df.fltr$VAR_2, df.fltr$variable)), ]
  df.fltr = df.fltr %>%
    dplyr::filter(frac < .20) %>%
    dplyr::group_by(VAR_2, variable) %>%
    dplyr::mutate(total = n())  %>%
    dplyr::filter(total == length(unique(df$VAR_1))) %>%
    data.frame()

  mat = as.data.frame(base::scale(df[ , -which(names(df) %in% c("VAR_1", "VAR_2"))]))
  df = cbind(
    df[ , which(names(df) %in% c("VAR_1", "VAR_2"))],
    mat
  )

  df.m = reshape2::melt(df, ids = c("VAR_1", "VAR_2"))

  df.m = df.m %>%
    dplyr::group_by(VAR_1, VAR_2, variable) %>%
    dplyr::summarise(ave = mean(value))

  # Remove features with negative values accros all samples
  df.m$ave[is.na(df.m$ave)] = 0
  rm.ftrs = df.m %>%
    dplyr::group_by(VAR_2, variable) %>%
    dplyr::filter(ave <= 0) %>%
    dplyr::group_by(variable) %>%
    dplyr::count()
  l = length(unique(df.m[["VAR_1"]])) * length(unique(df.m[["VAR_2"]]))
  df.m = df.m[!df.m$variable %in% rm.ftrs[rm.ftrs$n == l, ]$variable, ]

  # Cluster rows
  df.dc = df.m %>% reshape2::dcast(VAR_1 + VAR_2 ~ variable)
  colnames(df.dc)[1:2] = c("VAR_1", "VAR_2")
  full_dend <- as.dendrogram(
    hclust(dist( t(df.dc[ , -which(names(df.dc) %in% c("VAR_1", "VAR_2"))])  ))
  )
  df.m$variable = factor(df.m$variable, levels = labels(full_dend))

  df.m$ID = paste0(df.m[["VAR_2"]], "_", df.m$variable)
  r = table(df.fltr$variable)[table(df.fltr$variable) == l]
  df.m = df.m[!df.m$variable %in% names(r), ]


  obj.l = Split_Object(se.w, split.by = "nearest.cl", threads = 20)

  # # x = obj.l$"9"
  # res = lapply(obj.l, function(x){
  #   Idents(x) = var.1
  #   DefaultAssay(x) = "custom_UCell_score"
  #   lvls = levels(Idents(x))
  #   l = list()
  #   for (i in lvls) {
  #     res.wlx = FindMarkers(
  #       x, ident.1 = "non-CR", ident.2 = "CR", group.by = "BEST_RESPONSE_CONSENSUS",
  #       subset.ident = i,  fc.slot = "counts", logfc.threshold = 0
  #     )
  #     res.wlx$SIGNATURE = rownames(res.wlx)
  #     res.wlx$TIMEPOINT = i
  #     l[[i]] = res.wlx
  #   }
  #   res.wlx = do.call("rbind", l)
  #   res.wlx$p_val_adj = p.adjust(res.wlx$p_val, method = "bonferroni")
  #   res.wlx = res.wlx[order(res.wlx$SIGNATURE), ]
  #   res.wlx = res.wlx[res.wlx$p_val_adj < 0.05, ]
  #   res.wlx[abs(res.wlx$avg_log2FC) > lfc.thres, ]
  # })
  # wlx = data.table::rbindlist(res, idcol = "milestone")
  # wlx$ID = paste0(wlx$milestone, "_", wlx$SIGNATURE)
  #
  # df.m$avg_log2FC = wlx$avg_log2FC[match(df.m$ID, wlx$ID)]

  lbls = setNames(unique(gsub("\\.", " ", df.m$variable)), unique(df.m$variable))

  if(!is.null(max.value)){
    df.m$ave[df.m$ave >= max.value] = max.value
    df.m$ave[df.m$ave <= -max.value] = -max.value
  } else {
    thres = quantile(df.m$ave, .99)
    df.m$ave[df.m$ave > thres & !is.infinite(df.m$ave)] = thres
    df.m$ave[df.m$ave < -thres & !is.infinite(df.m$ave)] = -thres
    max.value = thres
  }

  ggplot(df.m) +
    geom_tile(aes(x = VAR_1, y =variable, fill=ave)) +
    # geom_point(data = subset(df.m, avg_log2FC > 0), aes(x = VAR_1, y = variable), size = 1, shape = 24) +
    # geom_point(data = subset(df.m, avg_log2FC < 0), aes(x = VAR_1, y = variable), size = 1, shape = 25) +
    geom_tile(data=df.m[is.na(df.m$ave), ], aes(x=VAR_1, y=variable, col="#BBBBBB"), fill = "#BBBBBB") +
    scale_color_manual(name="FDR<0.05", labels=NULL, values="#BBBBBB") +
    scale_fill_scico(
      palette = "vik", midpoint = 0, na.value = "#BBBBBB",
      begin = .1, end = .9, limits = c(-max.value, max.value)
    ) +
    facet_grid(~ VAR_2) +
    scale_y_discrete(labels= lbls) +
    guides(
      fill = guide_colorbar(
        title = "Scaled\nAverage\nEnrichment", order = 1,
        title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(4, 'lines'),
        ticks.linewidth = 1.5/.pt
      )
    ) + xlab(NULL) + ylab(NULL) +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1.05),
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1)),
      axis.ticks.y = element_blank(),
      legend.key.size = unit(9,"pt")
    ) +
    ggtitle(pl.title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Differential abundance analysis (DA)
# Tile plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
da_tile_pl = function(
    res.speckle,
    grp1 = "Gr_1.2",
    grp2 = "Gr_0",
    pl.title = NULL,
    alpha = 0.1,
    set.levels = NULL,
    tile.size = 4,
    stroke.size = 2
){

  if(!is.null(set.levels)){
    res.speckle$BaselineProp.clusters = factor(
      res.speckle$BaselineProp.clusters, levels = set.levels
    )
    drop.lvls = FALSE
  } else {
    drop.lvls = TRUE
  }

  ggplot(res.speckle, aes(x = BaselineProp.clusters, y = Group)) +
    geom_point(aes(shape = ifelse(P.Value >= alpha, NA, "s")), size = tile.size, stroke = stroke.size, color = "#228833") +
    geom_point(aes(color = FC), size = tile.size, shape = 15) +
    # scale_size(range = tile.size, breaks = breaks_pretty(3)) +
    scale_shape_manual(values = c(22), name = paste0("p < ", alpha), labels = c("yes", "")) +
    scale_color_scico(
      palette = "vik", midpoint = 0,
      begin = .1, end = .9, breaks = breaks_pretty(3),
      limits = c(-max(abs(res.speckle$FC)), max(abs(res.speckle$FC)))
    ) +
    mytheme(base_size = 8) +
    theme(
      panel.grid.major = element_line(colour = "grey80", linetype = 2, linewidth = .3),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
      axis.text.y = element_text(size = rel(1)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.ticks.length = unit(0.05, 'cm'),
      legend.position = "bottom",
      legend.title = element_text(margin = margin(l = 3, unit = "pt"), size = rel(1)),
      legend.margin = margin(t=-4, l = 5),
    ) +
    guides(
      size = guide_legend(title = "-log10(p-value)", title.position = "left", order = 2),
      shape = guide_legend(order = 3, title.position = "left", override.aes = list(size = 3, stroke = 1)),
      color = guide_colorbar(
        title = paste0("Log2FC (", grp1, " vs. ", grp2,")  "), order = 1,
        tticks.linewidth = .75/.pt, , frame.linewidth = 0.5/.pt,  title.vjust = .56, title.position = "left",
        barwidth = unit(4, 'lines'), barheight = unit(.35, 'lines')
      )
    ) +
    ggtitle(pl.title) +
    scale_x_discrete(drop = drop.lvls)
}
