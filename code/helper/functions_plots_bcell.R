####################################################################
# Longitudinal cells
####################################################################
plot_longitudinal_cells <- function(celltypes,label="B-cells in %",remove_pre_Tec=TRUE,prod="all",test=TRUE){
  df = data.frame(PATIENT_ID=se.meta$PATIENT_ID,TIMEPOINT=se.meta$TIMEPOINT,celltype=se.meta$celltype)
  if (remove_pre_Tec){
    df=df[df$TIMEPOINT!="pre_Tec",]
  }
  df = group_by(df,PATIENT_ID, TIMEPOINT, celltype) %>%
    summarise(CellCount = n(), .groups = 'drop')

  df.cell = df %>% group_by(PATIENT_ID, TIMEPOINT) %>%
    summarise(
      TotalCells = sum(CellCount),
      ct = sum(CellCount[celltype %in% celltypes]),
      .groups = 'drop'
    )
  df.cell$perc_ct = df.cell$ct/df.cell$TotalCells*100
  df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]
  df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]

  if (prod == "ide"){
    df.cell = df.cell %>% subset(PRODUCT == "ide")
    label = paste0(label," (ide)")
  } else if (prod == "cilta"){
    df.cell = df.cell %>% subset(PRODUCT == "cilta")
    label = paste0(label," (cilta)")
  }

  df.cell = df.cell %>% droplevels()
  if(test){
    p.vals = df.cell %>% group_by(TIMEPOINT) %>% pairwise_wilcox_test(perc_ct ~ BEST_RESPONSE_CONSENSUS,exact=T)
    print(p.vals);
    p.vals = p.vals %>% add_xy_position()
  }

  longitudinal_cells = ggplot(data = df.cell, aes(x=BEST_RESPONSE_CONSENSUS, y = perc_ct, fill = BEST_RESPONSE_CONSENSUS))+
    facet_wrap(~TIMEPOINT) +
    geom_boxplot(aes(x = BEST_RESPONSE_CONSENSUS, y = perc_ct, group = BEST_RESPONSE_CONSENSUS), alpha = 1) +
    geom_point(aes(x = BEST_RESPONSE_CONSENSUS, y = perc_ct)) +
    theme_bw(base_size = 18) +
    xlab("Timepoints")+
    ylab(label)+
    mytheme(base_size = base.size) +
    theme(axis.title.x = element_blank(),
          legend.position = "top")+
    scale_fill_manual(values = c("#6699CC", "#997700")) +
    labs(fill = "") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  if (test)
    longitudinal_cells = longitudinal_cells +  stat_compare_means(aes(label = paste("p =",after_stat(p.format))))

  return(longitudinal_cells)
}


####################################################################
# Longitudinal BCMA
####################################################################
plot_longitudinal_bcma <- function(types,label="BCMA in ng/mL"){
  df.cell = pdata.elisa
  df.cell$perc_ct = df.cell[,types]
  ftr = table(df.cell$Day)
  ftr = names(ftr)[ftr>20]
  df.cell = df.cell[df.cell$Day %in% ftr,]
  df.cell$Day = ordered(df.cell$Day,levels=c("LA","Day 0","Day 7","Day 30","Day 100"))
  df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell = df.cell[!is.na(df.cell$BEST_RESPONSE_CONSENSUS),]

  p.vals = df.cell %>% group_by(Day) %>% pairwise_wilcox_test(perc_ct ~ BEST_RESPONSE_CONSENSUS, alternative = "two.sided")
  print(p.vals)
  p.vals = p.vals %>% add_xy_position()

  longitudinal_cells = ggplot(data = df.cell, aes(x=BEST_RESPONSE_CONSENSUS, y = perc_ct, fill = BEST_RESPONSE_CONSENSUS))+
    facet_wrap(~Day) +
    geom_boxplot(aes(x = BEST_RESPONSE_CONSENSUS, y = perc_ct, group = BEST_RESPONSE_CONSENSUS), alpha = 1) +
    geom_point(aes(x = BEST_RESPONSE_CONSENSUS, y = perc_ct)) +
    stat_pvalue_manual(p.vals,inherit.aes = F) +
    theme_bw(base_size = 18) +
    xlab("Timepoints")+
    ylab(label)+
    mytheme(base_size = base.size) +
    theme(axis.title.x = element_blank(),
          legend.position = "top")+
    scale_fill_manual(values = c("#6699CC", "#997700")) +
    labs(fill = "") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(longitudinal_cells)
}


####################################################################
# Longitudinal BCMA per product
####################################################################
plot_longitudinal_bcma_prod <- function(types,label="sBCMA in ng/mL",prod="all"){
  df.cell = pdata.elisa
  df.cell$perc_ct = df.cell[,types]
  ftr = table(df.cell$Day)
  ftr = names(ftr)[ftr>20]
  df.cell = df.cell[df.cell$Day %in% ftr,]
  df.cell$Day = ordered(df.cell$Day,levels=c("LA","Day 0","Day 7","Day 30","Day 100"))
  df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell = df.cell[!is.na(df.cell$BEST_RESPONSE_CONSENSUS) & !is.na(df.cell$perc_ct),]

  if (prod=="ide"){
    df.cell = df.cell %>% subset(PRODUCT == "ide")
  } else if (prod=="cilta"){
    df.cell = df.cell %>% subset(PRODUCT == "cilta")
  }

  for (i in levels(df.cell$Day)){
    if (length(unique(df.cell$BEST_RESPONSE_CONSENSUS[df.cell$Day==i]))<2)
      df.cell = df.cell %>% subset(Day != i)
  }
  df.cell = df.cell %>% droplevels()

  p.vals = df.cell %>% group_by(Day) %>% pairwise_wilcox_test(perc_ct ~ BEST_RESPONSE_CONSENSUS, alternative = "two.sided")
  print(p.vals)
  p.vals = p.vals %>% add_xy_position()

  longitudinal_cells = ggplot(data = df.cell, aes(x=BEST_RESPONSE_CONSENSUS, y = perc_ct, fill = BEST_RESPONSE_CONSENSUS))+
    facet_wrap(~Day) +
    geom_boxplot(alpha = 1) +
    geom_point() +
    stat_pvalue_manual(p.vals,inherit.aes = F) +
    theme_bw(base_size = 18) +
    xlab("Timepoints")+
    ylab(paste0(label," (",prod,")"))+
    mytheme(base_size = base.size) +
    theme(axis.title.x = element_blank(),
          legend.position = "top")+
    scale_fill_manual(values = c("#6699CC", "#997700")) +
    #scale_y_log10() +
    labs(fill = "") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(longitudinal_cells)
}

####################################################################
# Longitudinal FcRH5 per product
####################################################################
plot_longitudinal_FcRH5_prod <- function(types,label="FcRH5",prod="all"){
  df.cell = pdata.elisa
  df.cell$perc_ct = df.cell[,types]
  ftr = table(df.cell$Day)
  ftr = names(ftr)[ftr>20]
  df.cell = df.cell[df.cell$Day %in% ftr,]
  df.cell$Day = ordered(df.cell$Day,levels=c("LA","Day 0","Day 7","Day 30","Day 100"))
  df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell = df.cell[!is.na(df.cell$BEST_RESPONSE_CONSENSUS) & !is.na(df.cell$perc_ct),]

  if (prod=="ide"){
    df.cell = df.cell %>% subset(PRODUCT == "ide")
  } else if (prod=="cilta"){
    df.cell = df.cell %>% subset(PRODUCT == "cilta")
  }

  for (i in levels(df.cell$Day)){
    if (length(unique(df.cell$BEST_RESPONSE_CONSENSUS[df.cell$Day==i]))<2)
      df.cell = df.cell %>% subset(Day != i)
  }
  df.cell = df.cell %>% droplevels()

  p.vals = df.cell %>% group_by(Day) %>% pairwise_wilcox_test(perc_ct ~ BEST_RESPONSE_CONSENSUS, alternative = "two.sided")
  p.vals = p.vals %>% add_xy_position()
  print(p.vals)

  longitudinal_cells =
    ggplot(data = df.cell, aes(x=BEST_RESPONSE_CONSENSUS, y = perc_ct, fill = BEST_RESPONSE_CONSENSUS))+
    facet_wrap(~Day) +
    geom_point() +
    geom_boxplot(alpha = 1) +
    scale_y_sqrt() +
    stat_compare_means() +
    theme_bw(base_size = 18) +
    xlab("Timepoints")+
    ylab(paste0(label," (",prod,")"))+
    mytheme(base_size = base.size) +
    theme(axis.title.x = element_blank(),
          legend.position = "top")+
    scale_fill_manual(values = c("#6699CC", "#997700")) +
    labs(fill = "") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(longitudinal_cells)
}

####################################################################
# Longitudinal FcRH5 comp product
####################################################################
plot_longitudinal_FcRH5_comp_prod <- function(types,label="FcRH5"){
  df.cell = pdata.elisa
  df.cell$perc_ct = df.cell[,types]
  ftr = table(df.cell$Day)
  ftr = names(ftr)[ftr>20]
  df.cell = df.cell[df.cell$Day %in% ftr,]
  df.cell$Day = ordered(df.cell$Day,levels=c("LA","Day 0","Day 7","Day 30","Day 100"))
  df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$patient_id,se.meta$PATIENT_ID)]
  df.cell = df.cell[!is.na(df.cell$BEST_RESPONSE_CONSENSUS) & !is.na(df.cell$perc_ct),]

  for (i in levels(df.cell$Day)){
    if (length(unique(df.cell$PRODUCT[df.cell$Day==i]))<2)
      df.cell = df.cell %>% subset(Day != i)
  }
  df.cell = df.cell %>% droplevels()

  p.vals = df.cell %>% group_by(Day) %>% pairwise_wilcox_test(perc_ct ~ PRODUCT, alternative = "two.sided")
  p.vals = p.vals %>% add_xy_position()
  print(p.vals)

  longitudinal_cells =
    ggplot(data = df.cell, aes(x=PRODUCT, y = perc_ct, fill = PRODUCT))+
    facet_wrap(~Day) +
    geom_point() +
    geom_boxplot(alpha = 1) +
    scale_y_sqrt() +
    stat_compare_means() +
    theme_bw(base_size = 18) +
    xlab("Timepoints")+
    ylab(label)+
    mytheme(base_size = base.size) +
    theme(axis.title.x = element_blank(),
          legend.position = "top")+
    scale_fill_manual(values = c("#4477AA","#CC6677")) +
    labs(fill = "") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(longitudinal_cells)
}


####################################################################
# Correlation of celltypes
####################################################################
plot_corr_ct <- function(df,leg=FALSE,lab=""){
  n = names(df)
  n = gsub(pattern="BCMA",replacement="sBCMA",n)
  n = gsub(pattern="perc_bcell",replacement="B-cell",n)
  n = gsub(pattern="perc_pDC",replacement="pDC",n)
  n = gsub(pattern="perc_plasma",replacement="Plasma",n)
  names(df) = n

  corr = cor(df,method="spearman",use="pairwise.complete.obs")
  p.mat <- ggcorrplot::cor_pmat(df,method="spearman",exact=FALSE)
  pt = ggcorrplot(corr = corr, lab=T, show.legend = leg ,#method = "circle",
                  type = "lower",
                  colors = c("#4477AA","#FFFFFF","#CC6677"),
                  outline.color = "#AAAAAAAA",
                  lab_size = 5,
                  insig = "pch",
                  pch.cex = 20,
                  pch.col = "#AAAAAAAA",
                  tl.cex = 15,
                  p.mat = p.mat) +
    ggtitle(lab)
  return(pt)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get sBCMA change data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_sBCMA_data <- function(pdata){

  clin_data = pdata$pdata.clin

  # CRS
  crs_data = pdata$pdata.crs
  crs_data$CRS_GRADE = factor(crs_data$CRS_GRADE, levels = c(0, 1, 2))
  crs_data$CRS_GRADE = ordered(crs_data$CRS_GRADE)

  # ELISA
  elisa_data = pdata$pdata.elisa
  elisa_data$SAMPLE_ID = paste0(elisa_data$PATIENT_ID, "_1")

  # elisa_data = elisa_data %>%
  #   dplyr::mutate(`BCMA [ng/mL]` = as.numeric(ifelse(`BCMA [ng/mL]` == ">> Y range", max(`BCMA [ng/mL]`), `BCMA [ng/mL]`)))

  sbcma_day_0 = elisa_data %>%
    dplyr::filter(DAY == "Day 0") %>%
    dplyr::mutate(BCMA_NG_ML = as.numeric(BCMA_NG_ML)) %>%
    dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
    dplyr::rename(sbcma_day0 = "BCMA_NG_ML")

  sbcma_la = elisa_data %>%
    dplyr::filter(DAY == "LA") %>%
    dplyr::mutate(BCMA_NG_ML = as.numeric(BCMA_NG_ML)) %>%
    dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
    dplyr::rename(sbcma_la = "BCMA_NG_ML")

  sbcma_day_30 = elisa_data %>%
    dplyr::filter(DAY == "Day 30") %>%
    dplyr::mutate(BCMA_NG_ML = as.numeric(BCMA_NG_ML)) %>%
    dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
    dplyr::rename(sbcma_day30 = "BCMA_NG_ML")

  # sBCMA fold change
  sbcma_change = crs_data[c("SAMPLE_ID", "CRS_GRADE", "CRP_MAX", "CRP_BASE")] %>%
    dplyr::left_join(clin_data[,c("SAMPLE_ID", "BEST_RESPONSE_CONSENSUS")]) %>%
    dplyr::left_join(sbcma_day_0) %>%
    dplyr::left_join(sbcma_day_30) %>%
    dplyr::mutate(change = sbcma_day30 / sbcma_day0)


  # Add IL-6
  elisa_data$IL_6_PG_ML[elisa_data$IL_6_PG_ML %in% c(">> std range",">> Y range")] = max(as.numeric(elisa_data$IL_6_PG_ML),na.rm=T)
  il6_day_0 = elisa_data %>%
    dplyr::filter(DAY == "Day 0") %>%
    dplyr::mutate(IL_6_PG_ML = as.numeric(IL_6_PG_ML)) %>%
    dplyr::select(SAMPLE_ID, IL_6_PG_ML) %>%
    dplyr::rename(il6_day0 = "IL_6_PG_ML")

  il6_day_7 = elisa_data %>%
    dplyr::filter(DAY == "Day 7") %>%
    dplyr::mutate(IL_6_PG_ML = as.numeric(IL_6_PG_ML)) %>%
    dplyr::select(SAMPLE_ID, IL_6_PG_ML) %>%
    dplyr::rename(il6_day7 = "IL_6_PG_ML")

  il6_la = elisa_data %>%
    dplyr::filter(DAY == "LA") %>%
    dplyr::mutate(IL_6_PG_ML = as.numeric(IL_6_PG_ML)) %>%
    dplyr::select(SAMPLE_ID, IL_6_PG_ML) %>%
    dplyr::rename(il6_la = "IL_6_PG_ML")

  il6_day_30 = elisa_data %>%
    dplyr::filter(DAY == "Day 30") %>%
    dplyr::mutate(IL_6_PG_ML = as.numeric(IL_6_PG_ML)) %>%
    dplyr::select(SAMPLE_ID, IL_6_PG_ML) %>%
    dplyr::rename(il6_day30 = "IL_6_PG_ML")

  il6_day_100 = elisa_data %>%
    dplyr::filter(DAY == "Day 100") %>%
    dplyr::mutate(IL_6_PG_ML = as.numeric(IL_6_PG_ML)) %>%
    dplyr::select(SAMPLE_ID, IL_6_PG_ML) %>%
    dplyr::rename(il6_day100 = "IL_6_PG_ML")

  # il6 fold change
  sbcma_change = sbcma_change %>%
    dplyr::left_join(il6_day_0) %>%
    dplyr::left_join(il6_la) %>%
    dplyr::left_join(il6_day_7) %>%
    dplyr::left_join(il6_day_30) %>%
    dplyr::left_join(il6_day_100) %>%
    dplyr::mutate(il6_diff = il6_day30 - il6_day0)

  imm.s = pdata$pdata.imm.s

    # Add CAR expansion
  max_car = imm.s %>%
    dplyr::group_by(SAMPLE_ID) %>%
    dplyr::slice_max(CD3_CAR_PERC,na_rm = TRUE,with_ties=FALSE) %>%
    dplyr::select(SAMPLE_ID, CD3_CAR_PERC)

  sbcma_change = sbcma_change %>%
    dplyr::left_join(max_car)
  sbcma_change = sbcma_change %>%
    dplyr::mutate(CRP_FC = CRP_MAX / CRP_BASE)

  return(sbcma_change)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# STARtac
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
startrac = function(obj, tp = "Late", celltype = "celltype_short_3", limit = 0){

  pd = obj@meta.data

  pd = pd[pd$TIMEPOINT == tp, ]
  pd$barcode = rownames(pd)
  pd$celltype = pd[[celltype]]
  pd = droplevels(pd)
  pd = pd %>% dplyr::select(
    "Cell_Name" = barcode, "clone.id" = CLONE_PSEUDO_ID, "patient" = orig.ident,
    "majorCluster" = celltype, "group" = BEST_RESPONSE_CONSENSUS
  )
  pd$loc = "T"
  pd$patient = as.character(pd$patient)
  suppressWarnings({
    suppressMessages({
      out <- Startrac.run(pd, proj="bcma", cores=NULL,verbose=T)
    })
  })

  Startrac::plot(out,index.type="cluster.all",byPatient=T)

  df = out@cluster.data
  df = df[df$aid != "bcma", ]
  df$NCells = NULL; df$migr = NULL
  df = reshape2::melt(df, keys = c("aid", "majorCluster"))
  df = df[df$variable == "expa", ]
  df = df[!is.na(df$value), ]
  df
}
