# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Modified function from the dynverse Framework
# https://github.com/dynverse
# https://www.nature.com/articles/s41587-019-0071-9
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ti_networks = function(sds = NULL, clust.res = "RNA_snn_res.0.1"){

  # sds = sce.dm
  X = SingleCellExperiment::reducedDim(sds, "DS")
  clus.labels = as.character(sds[[clust.res]])
  clusters = unique(clus.labels)

  for (j in names(sds$slingshot@metadata$curves)) {
    sds$slingshot@metadata$curves[[j]]$dist_ind = sds$slingshot@metadata$curves[[j]]$dist_ind[
      names(sds$slingshot@metadata$curves[[j]]$dist_ind) %in% colnames(sds)
    ]
    sds$slingshot@metadata$curves[[j]]$w = sds$slingshot@metadata$curves[[j]]$w[
      names(sds$slingshot@metadata$curves[[j]]$w) %in% colnames(sds)
    ]
    sds$slingshot@metadata$curves[[j]]$lambda = sds$slingshot@metadata$curves[[j]]$lambda[
      names(sds$slingshot@metadata$curves[[j]]$lambda) %in% colnames(sds)
    ]
  }

  .dist_clusters_full <- function(c1,c2){
    mu1 <- colMeans(c1)
    mu2 <- colMeans(c2)
    diff <- mu1 - mu2
    s1 <- cov(c1)
    s2 <- cov(c2)
    return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
  }
  .dist_clusters_diag <- function(c1,c2){
    mu1 <- colMeans(c1)
    mu2 <- colMeans(c2)
    diff <- mu1 - mu2
    if(nrow(c1)==1){
      s1 <-  diag(ncol(c1))
    }else{
      s1 <- diag(diag(cov(c1)))
    }
    if(nrow(c2)==1){
      s2 <-  diag(ncol(c2))
    }else{
      s2 <- diag(diag(cov(c2)))
    }
    return(as.numeric(t(diff) %*% solve(s1 + s2) %*% diff))
  }

  # determine the distance function
  min.clus.size <- min(table(clus.labels))
  if(min.clus.size <= ncol(X)){
    message('Using diagonal covariance matrix')
    dist.fun <- function(c1,c2) .dist_clusters_diag(c1,c2)
  }else{
    message('Using full covariance matrix')
    dist.fun <- function(c1,c2) .dist_clusters_full(c1,c2)
  }

  D <- sapply(clusters,function(clID1){
    sapply(clusters,function(clID2){
      clus1 <- X[clus.labels == clID1, ,drop = FALSE]
      clus2 <- X[clus.labels == clID2, ,drop = FALSE]
      return(dist.fun(clus1, clus2))
    })
  })
  rownames(D) <- clusters
  colnames(D) <- clusters
  D = reshape2::melt(D)
  colnames(D) = c("from", "to", "length")
  D$to = as.character(D$to)
  D$from = as.character(D$from)
  rownames(D) = paste0(D$from, D$to)

  # satisfy r cmd check
  from <- to <- NULL

  # collect milestone network
  lineages <- slingLineages(sds)
  lineage_ctrl <- slingParams(sds)
  cluster_network <- lineages %>%
    purrr::map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
    unique()
  cluster_network$ID = paste0(cluster_network$from, cluster_network$to)
  cluster_network = D[cluster_network$ID, ]
  cluster_network$directed = TRUE
  rownames(cluster_network)  = NULL

  # collect clusters
  cluster <- slingClusterLabels(sds)

  # collect progressions
  lin_assign <- apply(slingCurveWeights(sds), 1, which.max)

  progressions <- purrr::map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    lin <- lineages[[l]]
    pst.full <- slingPseudotime(sds, na = FALSE)[,l]
    pst <- pst.full[ind]
    means <- sapply(lin, function(clID){
      stats::weighted.mean(pst.full, cluster[,clID])
    })
    non_ends <- means[-c(1,length(means))]
    edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]

    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1

    tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
  })

  # Conversion between milestone percentages and progressions
  check_froms <- tapply(progressions$from, progressions$cell_id, function(x) length(unique(x)) == 1)
  if (any(!check_froms)) {
    stop("In ", sQuote("progressions"), ", cells should only have 1 unique from milestone.")
  }

  check_edges <- progressions %>%
    dplyr::left_join(cluster_network, by = c("from", "to")) %>%
    dplyr::left_join(cluster_network %>% select(to = from, from = to, length2 = length), by = c("from", "to"))

  if (any(is.na(check_edges$length) & is.na(check_edges$length2))) {
    stop("All from-to combinations in ", sQuote("progressions"), " should be in ", sQuote("milestone_network"), " as well.")
  }

  # determine milestone percentages for self edges
  selfs <- progressions %>%
    dplyr::filter(from == to) %>%
    dplyr::select(cell_id, milestone_id = from) %>%
    dplyr::mutate(percentage = 1)

  progressions.2 <- progressions %>%
    dplyr::filter(from != to)

  # determine milestone percentages for 'from' milestones
  from_mls <- tapply(progressions.2$from, progressions.2$cell_id, dplyr::first, default = NA_character_)
  from_pct <- 1 - tapply(progressions.2$percentage, progressions.2$cell_id, sum, default = NA_real_)
  froms <- tibble(
    cell_id = names(from_mls) %||% character(),
    milestone_id = from_mls[cell_id] %>% unname() %>% as.character(),
    percentage = from_pct[cell_id] %>% unname() %>% as.numeric()
  )

  # determine milestone percentages for 'to' milestones
  tos <- progressions.2 %>%
    dplyr::select(cell_id, milestone_id = to, percentage)

  # return all percentages
  cluster_percentages = dplyr::bind_rows(selfs, froms, tos)

  list(
    cell_ids = colnames(sds),
    cluster_network = cluster_network,
    cluster_percentages = cluster_percentages,
    progressions = progressions
  )

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_onto_nearest_cluster <- function(ti) {
  grouping <- ti$cluster_percentages %>%
    dplyr::group_by(cell_id) %>%
    dplyr::arrange(-percentage) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(percentage = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell_id, milestone_id) %>%
    deframe()

  cell_ids <- ti$cell_ids
  ifelse(cell_ids %in% names(grouping), grouping[cell_ids], NA) %>%
    set_names(cell_ids)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
group_onto_ti_edges <- function(ti, group_template = "{from}->{to}") {
  # first map cells to largest percentage (in case of divergence regions)
  progressions <-
    ti$progressions %>%
    group_by(cell_id) %>%
    arrange(-percentage) %>%
    slice(1) %>%
    ungroup()

  # do the actual grouping
  grouping <-
    progressions %>%
    group_by(from, to) %>%
    mutate(group_id = as.character(glue::glue(group_template))) %>%
    ungroup() %>%
    select(cell_id, group_id) %>%
    deframe()

  cell_ids <- ti$cell_ids
  ifelse(cell_ids %in% names(grouping), grouping[cell_ids], NA) %>%
    set_names(cell_ids)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ti_topology_ratio = function(
    sds = NULL,
    ti = NULL,
    target = "RESPONSE_CONSENSUS_2",
    grp1 = "non-CR",
    grp2 = "CR",
    max.value = NULL,
    permuations = 100
) {

  library(ggraph)
  library(tidygraph)

  data = data.frame(colData(sds)[, colnames(colData(sds)) != "slingshot"])
  # nearest.cl = group_onto_nearest_cluster(ti)
  # stopifnot(identical(names(nearest.cl), rownames(data)))
  # data$cluster = nearest.cl

  perm = parallel::mclapply(seq(1:permuations), function(x){

    set.seed(x)
    data.rndm = data
    data.rndm[[target]] = sample(data.rndm[[target]], replace = F)

    stats.rndm = data.rndm %>%
      dplyr::group_by(cluster, .data[[target]]) %>%
      dplyr::summarise(nbr.cells = n())

    sum.cells = stats.rndm %>%
      dplyr::group_by(.data[[target]]) %>%
      dplyr::summarise(nbr.cells = sum(nbr.cells))

    stats.rndm$gr = sum.cells$nbr.cells[match(stats.rndm[[target]], sum.cells[[target]])]
    stats.rndm$freq =  stats.rndm$nbr.cells / stats.rndm$gr

    ratio.rndm = stats.rndm %>% group_by(cluster) %>%
      reframe(ratio = freq[.data[[target]] == grp1]/freq[.data[[target]] == grp2])
    ratio.rndm$ratio = log2(ratio.rndm$ratio)
    ratio.rndm = ratio.rndm[!is.infinite(ratio.rndm$ratio), ]
    ratio.rndm = ratio.rndm[!is.na(ratio.rndm$ratio), ]
    ratio.rndm

  }, mc.cores = 40)

  perm.qu = do.call("rbind", perm) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(left = quantile(ratio, 0.005), right = quantile(ratio, 0.995))

  ###

  stats = data %>%
    dplyr::group_by(cluster, .data[[target]]) %>%
    dplyr::summarise(nbr.cells = n())

  sum.cells = stats %>%
    dplyr::group_by(.data[[target]]) %>%
    dplyr::summarise(nbr.cells = sum(nbr.cells))

  stats$gr = sum.cells$nbr.cells[match(stats[[target]], sum.cells[[target]])]
  stats$freq =  stats$nbr.cells / stats$gr
  ratio = stats %>% group_by(cluster) %>%
    reframe(ratio = freq[.data[[target]] == grp1] / freq[.data[[target]] == grp2])
  ratio$ratio = log2(ratio$ratio)

  ratio = merge(ratio, perm.qu, all.x = T)
  ratio = ratio %>% mutate(
    sign = case_when(
      (ratio > 0 & ratio > right) & ratio > log2(1.25) ~ "* ",
      (ratio < 0 & ratio < left) & ratio < -log2(1.25) ~ "* ",
      TRUE ~ ""
    )
  )
  # print(ratio)

  ###

  car.df = data.frame(cluster = as.character(unique(data$cluster)))
  car.t = data %>%
    dplyr::filter(CAR_BY_EXPRS == "TRUE") %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(nbr.cells = n()) %>%
    dplyr::mutate(fraction = nbr.cells / sum(nbr.cells)) %>%
    data.frame()
  car.t = merge(car.df, car.t, all.x = T)
  car.t$fraction[is.na(car.t$fraction)] = 0
  car.t$fraction = car.t$fraction * 100
  car.t$fraction = round(car.t$fraction, 0)

  ###

  milestone_graph <- tidygraph::as_tbl_graph(ti$cluster_network)
  milestone_positions <- milestone_graph %>%
    ggraph::create_layout("kk") %>%
    dplyr::mutate(milestone_id = as.character(.data$name))

  milestone_positions$ratio = ratio$ratio[match(milestone_positions$name, ratio$cluster)]
  milestone_positions$car = car.t$fraction[match(milestone_positions$name, car.t$cluster)]
  milestone_positions$sign = ratio$sign[match(milestone_positions$name, ratio$cluster)]

  milestone_graph <- igraph::graph_from_data_frame(
    ti$cluster_network,
    vertices = milestone_positions %>% dplyr::select(-x, -y)
  ) %>%
    tidygraph::as_tbl_graph()

  ###

  max.v = max(abs(ratio$ratio))
  print(max.v)
  if(!is.null(max.value)){
    max.v = max.value
  }
  ext.axis = 20
  max.x = max(milestone_positions$x)
  max.x = max.x + (max.x / 100 * ext.axis)
  min.x = min(milestone_positions$x)
  min.x = min.x + (min.x / 100 * ext.axis)

  max.y = max(milestone_positions$y)
  max.y = max.y + (max.y / 100 * ext.axis)
  min.y = min(milestone_positions$y)
  min.y = min.y + (min.y / 100 * ext.axis)

  ###

  col.p = scico(10, palette = "vik")


  ggraph::ggraph(
    milestone_graph, "manual",
    x = milestone_positions$x, y = milestone_positions$y
  ) +
    ggraph::geom_edge_fan(width = .3) +
    ggraph::geom_edge_fan(
      aes(
        xend = .data$x + (.data$xend-.data$x)/1.5,
        yend = .data$y + (.data$yend-.data$y)/1.5
      ),
      width = .3,
      arrow = grid::arrow(type = "closed", length = unit(0.1, "cm"))
    ) +
    ggraph::geom_node_label(
      aes(fill = .data$ratio, label = .data$car),
      color = "white", size = 2.75, label.padding = unit(0.15, "lines")
    ) +
    geom_text(
      data = milestone_positions, aes(x, y, label = sign, hjust = -.3, vjust = .3),
      size = 6, color = "black"
    ) +
    scale_fill_gradientn(
      colours = c(rep(col.p[1], 2), col.p, rep(col.p[10], 2)),
      limits = c(-max.v, max.v), n.breaks = 4
    ) +
    # scico::scale_fill_scico(
    #   palette = "vik", direction = 1, midpoint = 0,
    #   limits = c(-max.v, max.v), n.breaks = 3
    # ) +
    guides(
      fill = guide_colorbar(
        title = "log2(% non-CR / % CR)", title.vjust = 1.05,
        barwidth = unit(4,'lines'), position = "bottom",
        barheight = unit(.35, 'lines'),
        ticks.linewidth = .75/.pt, frame.linewidth = 0.5/.pt
      )
    ) +
    ylim(min.y, max.y) + xlim(min.x, max.x) +
    mytheme() +
    theme(
      legend.ticks.length = unit(0.05, 'cm'),
      legend.position = "bottom",
      legend.spacing.x = unit(12, 'pt'),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ti_topology_ct_comp = function(
    sds = NULL,
    ti = NULL,
    target = "celltye",
    fill.col = til.col,
    swtch.nodes = NULL
) {

  library(ggraph)
  library(tidygraph)
  library(scatterpie)

  data = data.frame(colData(sds)[, colnames(colData(sds)) != "slingshot"])
  data$target = data[[target]]
  nearest.cl = group_onto_nearest_cluster(ti)
  stopifnot(identical(names(nearest.cl), rownames(data)))
  data$cluster = nearest.cl

  ct.comp = data %>%
    group_by(cluster, target) %>%
    dplyr::summarise(nbr.cells = n()) %>%
    dplyr::mutate(fraction = nbr.cells / sum(nbr.cells)) %>%
    data.frame()

  # reshape::cast(ct.comp %>% select(-fraction), cluster ~ celltype)
  ct.comp = reshape::cast(ct.comp %>% select(-nbr.cells), cluster ~ target)
  ct.comp[is.na(ct.comp)] = 0

  milestone_graph <- tidygraph::as_tbl_graph(ti$cluster_network)
  milestone_positions <- milestone_graph %>%
    ggraph::create_layout("kk") %>%
    dplyr::mutate(milestone_id = as.character(.data$name))

  milestone_positions = merge(milestone_positions, ct.comp, by.x = "name", by.y = "cluster", all.x = T)

  milestone_graph <- igraph::graph_from_data_frame(
    ti$cluster_network,
    vertices = milestone_positions %>% select(-.data$x, -.data$y)
  ) %>%
    tidygraph::as_tbl_graph()


  if(!is.null(swtch.nodes)){
    swtch <- function(x,i,j) {x[c(i,j)] <- x[c(j,i)]; x}
    v = milestone_positions$.ggraph.index
    milestone_positions$new.ggraph.index = swtch(
      v,
      seq(1:length(v))[v %in% swtch.nodes][1],
      seq(1:length(v))[v %in% swtch.nodes][2]
    )
  } else {
    milestone_positions$new.ggraph.index = milestone_positions$.ggraph.index
  }

  ggraph::ggraph(milestone_graph, "manual", x = milestone_positions$x, y = milestone_positions$y) +
    ggraph::geom_edge_fan(width = .3) +
    ggraph::geom_edge_fan(
      aes(
        xend = .data$x + (.data$xend-.data$x)/1.5,
        yend = .data$y + (.data$yend-.data$y)/1.5
      ),
      width = .3,
      arrow = grid::arrow(type = "closed", length = unit(.1, "cm"))
    ) +
    geom_scatterpie(
      cols = colnames(ct.comp)[2:length(colnames(ct.comp))],
      data = milestone_positions,
      colour = NA,
      pie_scale = 3.5,
      donut_radius=.6,
      bg_circle_radius = .6
    ) +
    geom_text(data = milestone_positions, aes(x, y, label = new.ggraph.index), size = 2.5, color = "black") +
    scale_fill_manual(values = fill.col) +
    mytheme() +
    theme(
      # legend.position = "none",
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cluster_nearest_milestones = function(sce.tr, se.obj){
  ti.net = ti_networks(sds = sce.tr, clust.res = "cluster")
  nearest.cl = group_onto_nearest_cluster(ti.net)
  nearest.cl = data.frame(nearest.cl)
  miles.graph = tidygraph::as_tbl_graph(ti.net$cluster_network)
  miles.pos <- miles.graph %>%
    ggraph::create_layout("kk") %>%
    dplyr::mutate(milestone_id = as.character(.data$name))
  nearest.cl$nearest.cl = miles.pos$.ggraph.index[match(nearest.cl$nearest.cl, miles.pos$milestone_id)]
  nearest.cl$nearest.cl = factor(nearest.cl$nearest.cl)
  AddMetaData(se.obj, nearest.cl)
}
