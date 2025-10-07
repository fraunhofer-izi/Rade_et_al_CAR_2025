# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Survival analysis with cell type fractions
# Multivariate cox-regression model
# Adjusted for CAR Product
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# file = pd.tp$`Very Late`

cox.res.multi = parallel::mclapply(pd.tp, function(file) {

  TP = as.character(file$TIMEPOINT[1])
  pd = droplevels(file)

  prop = table(pd$celltype, pd$orig.ident)
  prop = prop[rowSums(prop) > 200, ]
  cell.nbr = rowSums(prop)
  prop = prop.table(prop, margin = 2)
  prop = as.data.frame.matrix(prop)
  prop = as.matrix(prop)
  prop = standardize(prop)

  pheno = pd
  pheno = pheno[!duplicated(pheno$orig.ident), ]
  rownames(pheno) = pheno$orig.ident
  pheno = pheno[intersect(rownames(pheno), colnames(prop)), ]
  prop = prop[, intersect(rownames(pheno), colnames(prop))]
  stopifnot(identical(rownames(pheno), colnames(prop)))

  # Survival object
  survivalDOD = Surv(pheno$PFS, pheno$PROGRESSION)
  survModel = survfit(survivalDOD ~ 1, data = pheno)

  # Multivariate Cox proportional hazard model
  car.product = pheno$PRODUCT
  exprs.mat.spl = split(as.matrix(prop), rownames(prop))

  cox.multi.l = list()
  for (celltype in names(exprs.mat.spl)) {
    nfit = summary(survival::coxph(survivalDOD ~ exprs.mat.spl[[celltype]] + car.product))
    res = data.frame(
      "HR" = nfit$coef[ ,2 , drop = T],
      "logHR" = nfit$coef[, 1, drop = T],
      "SE_logHR" = nfit$coef[, 3, drop = T],
      "L95CI" = nfit$conf.int[, "lower .95", drop = T],
      "U95CI" = nfit$conf.int[, "upper .95", drop = T],
      "Pval" = nfit$coef[, 5,  drop = T]
    )
    res = cbind(Covariate = c(celltype, "CAR Product"), res)
    res$CELLTYPE = celltype
    cox.multi.l[[celltype]] = res
  }

  cox.multi = do.call("rbind", cox.multi.l)
  cox.multi$TIMEPOINT = TP
  rownames(cox.multi) = NULL
  cox.multi$Pval_adj = p.adjust(cox.multi$Pval, method = "BH")
  cox.multi

}, mc.cores = 1)

df.surv.multi = do.call("rbind", cox.res.multi)
df.surv.multi$TIMEPOINT = factor(df.surv.multi$TIMEPOINT, levels = c("LP", "Late", "Very Late"))
df.surv.multi$ID = paste0(df.surv.multi$CELLTYPE, "_", df.surv.multi$TIMEPOINT)

###

df = df.surv.multi
df = df[df$Covariate != "CAR Product", ]
df$L95CI = exp(df$logHR - df$SE_logHR * 1.65)
df$U95CI = exp(df$logHR + df$SE_logHR * 1.65)


df = df[df$CELLTYPE %in% order.ct, ]
df$CELLTYPE = factor(df$CELLTYPE, levels = rev(order.ct))

df$TIMEPOINT = gsub("Apheresis", "LP", df$TIMEPOINT)
df$TIMEPOINT = factor(df$TIMEPOINT, levels = c("LP", "Late", "Very Late"))
df$Cox_Pval_Sign = ifelse(df$Pval < 0.1, "#549C74", "#555555")

dotCOLS = c("#9c9c9c", "#6699CC", "#004488")
barCOLS = c("#9c9c9c", "#6699CC", "#004488")

l = list()
for (i in 1:length(df$Pval)) {
  if (df$Pval[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval[i] > 0.0001 & df$Pval[i] <= 0.001) {
    l[[i]] = format(round(df$Pval[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_pl = unlist(l)
df$Cox_Pval_pl = ifelse(df$Pval < 0.1, paste0(df$Cox_Pval_pl, " *"), df$Cox_Pval_pl)

l = list()
for (i in 1:length(df$Pval_adj)) {
  if (df$Pval_adj[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval_adj[i] > 0.0001 & df$Pval_adj[i] <= 0.001) {
    l[[i]] = format(round(df$Pval_adj[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval_adj[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_adj_pl = unlist(l)
df$Cox_Pval_adj_pl = ifelse(
  df$Pval_adj < 0.1, paste0(df$Cox_Pval_adj_pl, " *"), df$Cox_Pval_adj_pl
)

hr.multi.1.pl =
ggplot(df, aes(x = CELLTYPE, y = HR, ymin = L95CI, ymax = U95CI, col = TIMEPOINT, fill=TIMEPOINT)) +
  geom_hline(yintercept=1, lty=2, lwd = .5, col = "#555555") +
  geom_linerange(linewidth=1.5, position=position_dodge(width = 0.8)) +
  facet_wrap(~ CELLTYPE, ncol = 1, scales = "free_y") +
  geom_point(
    size=2.5, shape=21, colour="white", stroke = 1, position=position_dodge(width = 0.8)
  ) +
  geom_text(
    data = df[!df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = U95CI + 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = -0.1, vjust = 0.5, size = rel(2.5)
  ) +
  geom_text(
    data = df[df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = L95CI - 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = 1.1, vjust = 0.5, size = rel(2.5)
  ) +
  coord_flip() +
  scale_colour_manual(values = alpha(barCOLS, .7))+
  scale_fill_manual(values = dotCOLS)+
  scale_y_log10(
    name= "Hazard ratio (90% CI)",
    limits = c(min(df$L95CI), max(df$U95CI) + 20)
  ) +
  xlab(NULL) +
  mytheme() +
  theme(
    strip.text.x  = element_blank(),
    # panel.spacing  = unit(.25, "lines"),
    legend.title=element_blank(),
    legend.margin = margin(t=-3),
    # legend.box.margin=margin(-5,0,5,0),
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    legend.text = element_text(margin = margin(l = 3, unit = "pt"), size = rel(1))
  ) +
  scale_x_discrete(expand = c(0.5, 0)) +
  ggtitle("Adjusted for CAR Product")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Survival analysis with cell type fractions
# Multivariate cox-regression model
# Adjusted for CAR Product + EMD + sBCMA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# file = pd.tp$`Very Late`

clin_data = readRDS(
  "publication/clinicial_data/clinical_table_DF_2024_10_28.Rds"
)$pdata.clin

# EMD
clin_data$EMD_BEFORE_CART = tidyr::replace_na(clin_data$EMD_BEFORE_CART, 0)
clin_data$EMD_BEFORE_CART = ifelse(clin_data$EMD_BEFORE_CART == 1, "grp1", "grp0")
clin_data$EMD_BEFORE_CART #grp0=no grp1=yes

# sBCMA groups
elisa_data = readRDS(
  "publication/clinicial_data/clinical_table_DF_2024_10_28.Rds"
)$pdata.elisa
elisa_data$SAMPLE_ID = paste0(elisa_data$PATIENT_ID, "_1")

sbcma_day_0 = elisa_data %>%
  dplyr::filter(DAY == "Day 0") %>%
  dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
  dplyr::rename(sbcma_day0 = "BCMA_NG_ML")

median_sbcma_day0 = median(
  as.numeric(subset(elisa_data, DAY == "Day 0")$BCMA_NG_ML), na.rm = T
)

clin_data = clin_data %>%
  dplyr::left_join(sbcma_day_0) %>%
  dplyr::mutate(
    sbcma_day0_bin = factor(ifelse(sbcma_day0 > median_sbcma_day0, "High", "Low"), levels = c("Low", "High"))
  )

cox.res.multi = parallel::mclapply(pd.tp, function(file) {

  TP = as.character(file$TIMEPOINT[1])
  pd = droplevels(file)

  prop = table(pd$celltype, pd$orig.ident)
  prop = prop[rowSums(prop) > 200, ]
  cell.nbr = rowSums(prop)
  prop = prop.table(prop, margin = 2)
  prop = as.data.frame.matrix(prop)
  prop = as.matrix(prop)
  prop = standardize(prop)

  pheno = pd
  pheno = pheno[!duplicated(pheno$orig.ident), ]
  rownames(pheno) = pheno$orig.ident

  pheno$EMD = clin_data$EMD_BEFORE_CART[match(
    pheno$PATIENT_ID, clin_data$PATIENT_ID
  )]
  pheno = pheno[!is.na(pheno$EMD), ]

  pheno$SBCMA = clin_data$sbcma_day0_bin[match(
    pheno$PATIENT_ID, clin_data$PATIENT_ID
  )]
  pheno = pheno[!is.na(pheno$SBCMA), ]

  pheno = pheno[intersect(rownames(pheno), colnames(prop)), ]
  prop = prop[, intersect(rownames(pheno), colnames(prop))]
  stopifnot(identical(rownames(pheno), colnames(prop)))

  # Survival object
  survivalDOD = Surv(pheno$PFS, pheno$PROGRESSION)
  survModel = survfit(survivalDOD ~ 1, data = pheno)

  # Multivariate Cox proportional hazard model
  car.product = pheno$PRODUCT
  emd = pheno$EMD
  sbcma = pheno$SBCMA
  exprs.mat.spl = split(as.matrix(prop), rownames(prop))

  cox.multi.l = list()
  for (celltype in names(exprs.mat.spl)) {
    nfit = summary(survival::coxph(survivalDOD ~ exprs.mat.spl[[celltype]] + car.product + emd + sbcma))
    res = data.frame(
      "HR" = nfit$coef[ ,2 , drop = T],
      "logHR" = nfit$coef[, 1, drop = T],
      "SE_logHR" = nfit$coef[, 3, drop = T],
      "L95CI" = nfit$conf.int[, "lower .95", drop = T],
      "U95CI" = nfit$conf.int[, "upper .95", drop = T],
      "Pval" = nfit$coef[, 5,  drop = T]
    )
    res = cbind(Covariate = c(celltype, "CAR Product", "EMD", "sBCMA"), res)
    res$CELLTYPE = celltype
    cox.multi.l[[celltype]] = res
  }

  cox.multi = do.call("rbind", cox.multi.l)
  cox.multi$TIMEPOINT = TP
  rownames(cox.multi) = NULL
  cox.multi$Pval_adj = p.adjust(cox.multi$Pval, method = "BH")
  cox.multi

}, mc.cores = 1)

df.surv.multi = do.call("rbind", cox.res.multi)
df.surv.multi$TIMEPOINT = factor(df.surv.multi$TIMEPOINT, levels = c("LP", "Late", "Very Late"))
df.surv.multi$ID = paste0(df.surv.multi$CELLTYPE, "_", df.surv.multi$TIMEPOINT)

###

df = df.surv.multi
df = df[df$Covariate != "CAR Product", ]
df = df[df$Covariate != "EMD", ]
df = df[df$Covariate != "sBCMA", ]
df$L95CI = exp(df$logHR - df$SE_logHR * 1.65)
df$U95CI = exp(df$logHR + df$SE_logHR * 1.65)


df = df[df$CELLTYPE %in% order.ct, ]
df$CELLTYPE = factor(df$CELLTYPE, levels = rev(order.ct))

df$TIMEPOINT = gsub("Apheresis", "LP", df$TIMEPOINT)
df$TIMEPOINT = factor(df$TIMEPOINT, levels = c("LP", "Late", "Very Late"))
df$Cox_Pval_Sign = ifelse(df$Pval < 0.1, "#549C74", "#555555")

dotCOLS = c("#9c9c9c", "#6699CC", "#004488")
barCOLS = c("#9c9c9c", "#6699CC", "#004488")

l = list()
for (i in 1:length(df$Pval)) {
  if (df$Pval[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval[i] > 0.0001 & df$Pval[i] <= 0.001) {
    l[[i]] = format(round(df$Pval[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_pl = unlist(l)
df$Cox_Pval_pl = ifelse(df$Pval < 0.1, paste0(df$Cox_Pval_pl, " *"), df$Cox_Pval_pl)

l = list()
for (i in 1:length(df$Pval_adj)) {
  if (df$Pval_adj[i] <= 0.0001) {
    l[[i]] = formatC(df$Pval[i], format = "e", digits= 2)
  } else if (df$Pval_adj[i] > 0.0001 & df$Pval_adj[i] <= 0.001) {
    l[[i]] = format(round(df$Pval_adj[i], digits=5), nsmall = 5)
  } else {
    l[[i]] = format(round(df$Pval_adj[i], digits=3), nsmall = 3)
  }
}
df$Cox_Pval_adj_pl = unlist(l)
df$Cox_Pval_adj_pl = ifelse(
  df$Pval_adj < 0.1, paste0(df$Cox_Pval_adj_pl, " *"), df$Cox_Pval_adj_pl
)


hr.multi.2.pl =
ggplot(df, aes(x = CELLTYPE, y = HR, ymin = L95CI, ymax = U95CI, col = TIMEPOINT, fill=TIMEPOINT)) +
  geom_hline(yintercept=1, lty=2, lwd = .5, col = "#555555") +
  geom_linerange(linewidth=1.5, position=position_dodge(width = 0.8)) +
  facet_wrap(~ CELLTYPE, ncol = 1, scales = "free_y") +
  geom_point(
    size=2.5, shape=21, colour="white", stroke = 1, position=position_dodge(width = 0.8)
  ) +
  geom_text(
    data = df[!df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = U95CI + 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = -0.1, vjust = 0.5, size = rel(2.5)
  ) +
  geom_text(
    data = df[df$CELLTYPE %in% levels(df$CELLTYPE)[1:1], ],
    aes(
      y = L95CI - 0.01, label = paste0(Cox_Pval_pl, " | ", Cox_Pval_adj_pl)
    ),
    color = "black", position=position_dodge(width = .8),
    hjust = 1.1, vjust = 0.5, size = rel(2.5)
  ) +
  coord_flip() +
  scale_colour_manual(values = alpha(barCOLS, .7))+
  scale_fill_manual(values = dotCOLS)+
  scale_y_log10(
    name= "Hazard ratio (90% CI)",
    limits = c(min(df$L95CI), max(df$U95CI) + 20)
  ) +
  xlab(NULL) +
  mytheme() +
  theme(
    strip.text.x  = element_blank(),
    # panel.spacing  = unit(.25, "lines"),
    legend.title=element_blank(),
    legend.margin = margin(t=-3),
    # legend.box.margin=margin(-5,0,5,0),
    legend.position = "bottom",
    axis.ticks.y = element_blank(),
    legend.text = element_text(margin = margin(l = 3, unit = "pt"), size = rel(1))
  ) +
  scale_x_discrete(expand = c(0.5, 0)) +
  ggtitle("Adjusted for CAR Product, EMD and sBCMA")
