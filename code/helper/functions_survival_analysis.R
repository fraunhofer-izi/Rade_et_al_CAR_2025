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
    ...
) {

  surv.obj <<- obj

  formula = as.formula(paste0("Surv(PFS_M, PROGRESSION)~", group))
  survModel = survfit(formula, data=surv.obj)
  survModel$call$formula = formula

  surv.pl = ggsurvplot(
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
    risk.table.fontsize = 2.5,
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

  surv.pl$plot = surv.pl$plot + annotate("richtext", x = x_anno, y = y_anno, hjust = 0, vjust=0, label = anno, size = 2.3, label.color=NA, fill=NA)

  surv.pl$table = surv.pl$table +
    theme(
      panel.border = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank()
    ) +
    coord_cartesian(xlim=c(0, cutoff)) +
    scale_x_continuous(expand = c(0, 2), breaks = seq(0, cutoff, by = x.breaks))

  # surv.gg = surv.pl$plot + inset_element(surv.pl$table, left = 0, bottom = 0, right = 1, top =.32)

  surv.pl$plot + inset_element(surv.pl$table, left = -.01, bottom = -0.03, right = 1.0075, top =.3)

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
    format(round(tabtext[, 1, drop = F], digits=2)),
    target
  )

  coi = log(summary(nfit)$conf.int[, 3:4])
  coi = matrix(coi, ncol=2)
  coi = format(round(coi, digits=2))

  lhr = paste0(tabtext[, 1], " [", coi[, 1], ", ", coi[, 2], "]")
  tabtext[, 1] = lhr
  rownames(tabtext) = leg.txt
  colnames(tabtext) = c("  logHR [95% CI]  ", "  p-value  ")

  ggtexttable(tabtext, theme = ttheme("light", base_size = 8))

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Table functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pvalue <- function(x, ...) {
  x <- x[-length(x)]  # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform an ANOVA
    # p <- summary(aov(y ~ g))[[1]][["Pr(>F)"]][1]
    # For numeric variables, perform Wilcoxon's rank-sum test
    p <- wilcox.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- fisher.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           # "Median (SD)"=sprintf("%s (&plusmn; %s)", MEDIAN, SD),
                                                           "Median (IQR)"=sprintf("%s [%s - %s]", MEDIAN, Q1, Q3),
                                                           "Range" = sprintf("%s - %s", MIN, MAX)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}
