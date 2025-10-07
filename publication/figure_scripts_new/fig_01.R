.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "cowplot", "patchwork", "readxl", "clinfun", "tidyr", "survival", "survminer",
  "ggtext", "prodlim", "table1", "rstatix", "ggarrow", "ggstar", "car", "scico",
  "openxlsx"
)
.bioc_packages = c()

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

source("code/helper/styles.R")
source("code/helper/functions.R")
theme_set(mytheme(9))
palette_response = c(CR = "#7B9AB6", `VGPR/PR` = "#E9C54E" , `SD/PD`="#9B740A")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
pdata = readRDS("publication/clinicial_data/clinical_table_DF_2024_10_28.Rds")

response_groups_str = "'nCR'='CR'; 'VGPR'='VGPR/PR'; 'PR'='VGPR/PR'; 'SD'='SD/PD'; 'PD'='SD/PD'"
response_levels = c("CR", "VGPR/PR", "SD/PD")

clin_data = pdata$pdata.clin

clin_data = clin_data %>%
  mutate(
    STATUS_BEFORE_CAR = factor(car::recode(STATUS_BEFORE_CAR, response_groups_str), levels = response_levels),
    BEST_RESPONSE = factor(car::recode(BEST_RESPONSE, response_groups_str), levels = response_levels),
    RESPONSE_1_M = factor(car::recode(RESPONSE_1_M, response_groups_str), levels = response_levels),
    RESPONSE_3_M = factor(car::recode(RESPONSE_3_M, response_groups_str), levels = response_levels),
    RESPONSE_6_M = factor(car::recode(RESPONSE_6_M, response_groups_str), levels = response_levels),
    RESPONSE_12_M = factor(car::recode(RESPONSE_12_M, response_groups_str), levels = response_levels),
  )

# EMD
clin_data$EMD_BEFORE_CART = tidyr::replace_na(clin_data$EMD_BEFORE_CART, 0)
clin_data$EMD_BEFORE_CART = ifelse(clin_data$EMD_BEFORE_CART == 1, "grp1", "grp0")

# Refractoriness
clin_data$REFRACTORINESS_BIN = factor(
  ifelse(clin_data$REFRACTORINESS_GROUP=="TCexposed", "TCexposed", "TCRRMM/PentaRRMM"),
  levels=c("TCexposed", "TCRRMM/PentaRRMM")
)

elisa_data = pdata$pdata.elisa
elisa_data$SAMPLE_ID = paste0(elisa_data$PATIENT_ID, "_1")

imm.s = pdata$pdata.imm.s
imm.s$CAR_RATIO = imm.s$CD4_CAR_PERC / imm.s$CD8_CAR_PERC

# START ### Supp Table
imm.s.supps = imm.s
imm.s.supps = imm.s.supps %>% dplyr::select(
  PATIENT_ID, DAY, CD3_CAR_PERC, CD3_CAR_ABS, CD4_CAR_PERC, CD4_CAR_ABS,
  CD8_CAR_PERC, CD8_CAR_ABS
)
imm.s.supps$PATIENT_ID = gsub("Patient0", "P", imm.s.supps$PATIENT_ID)
imm.s.supps$DAY[imm.s.supps$DAY == "Leukapheresis"] = "LP"

sheet = "Suppl_Table_1"
xlsx.filename = "publication/supplementary_info_new/table_fc_car_subsets.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  imm.s.supps,
  startRow = 1, startCol = 1)
saveWorkbook(wb, xlsx.filename, overwrite = T)
# END ### Supp Table

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1b-c: Swim plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
swim_data = clin_data %>%
  dplyr::mutate(Product = PRODUCT) %>%
  dplyr::mutate(
    sample_id = factor(
      gsub("atient0", "", gsub("_1", "", SAMPLE_ID)),
      levels = (lapply(clin_data %>% arrange(OS) %>% pull(SAMPLE_ID), function(x) gsub("atient0", "", gsub("_1", "", x))) %>% unlist)
    )
  )

# get response at each timepoint
d30 = swim_data %>%
  dplyr::mutate(
    start = -50,
    end = case_when(
      PFS<30 & OS < 30 ~ OS,
      PFS<45 & OS > 30 ~ PFS,
      .default = 30)
  ) %>%
  dplyr::select(sample_id, Product, STATUS_BEFORE_CAR, start, end)
colnames(d30)[3] = "response"

d90 = swim_data %>%
  dplyr::filter(PFS > 30) %>%
  dplyr::mutate(
    response = as.character(RESPONSE_1_M),
    start = ifelse(PFS<45, PFS, 30),
    end = case_when(
      PFS<90 & OS < 90 ~ OS,
      PFS<90 & OS > 90 & RESPONSE_1_M == "SD/PD" ~ OS,
      PFS<110 & OS > 90 ~ PFS,
      .default = 90)
  ) %>%
  dplyr::mutate(response = factor(
    dplyr::case_when(
      PFS<30 ~ "SD/PD",
      PFS<90 ~ RESPONSE_1_M,
      .default = response),
    levels = response_levels)
  ) %>%
  dplyr::select(sample_id, Product, response, start, end)

d180 = swim_data %>%
  dplyr::filter(PFS > 90) %>%
  dplyr::mutate(
    response = as.character(RESPONSE_3_M),
    start = ifelse(PFS<110, PFS, 90),
    end = case_when(
      PFS<180  ~ PFS,
      # PFS<180 & OS < 180 & ~ OS,
      .default = 180)
  ) %>%
  dplyr::mutate(response = factor(
    dplyr::case_when(
      PFS<90 ~ "SD/PD",
      PFS<180 ~ RESPONSE_3_M,
      .default = response),
    levels = response_levels)
  ) %>%
  dplyr::select(sample_id, Product, response, start, end)

d365 = swim_data %>%
  dplyr::filter(PFS > 180) %>%
  dplyr::mutate(
    response = as.character(RESPONSE_6_M),
    start = ifelse(PFS<180, PFS, 180),
    end = case_when(
      PFS<365 & OS < 365 ~ OS,
      PFS<365 & OS > 365 ~ PFS,
      PFS>365 & OS > 365 & PROGRESSION == 1 ~ PFS,
      .default = 365)
  ) %>%
  dplyr::mutate(response = factor(
    dplyr::case_when(
      PFS<180 ~ "SD/PD",
      PFS<365 ~ RESPONSE_6_M,
      is.na(response) ~ RESPONSE_6_M,
      .default = response),
    levels = response_levels)
  ) %>%
  dplyr::select(sample_id, Product, response, start, end)

lfu = swim_data %>%
  dplyr::mutate(
    response = RESPONSE_12_M,
    start = ifelse(PROGRESSION == 1, PFS, 365),
    end = OS
  ) %>%
  dplyr::filter(start<end, sample_id != "P64") %>%
  dplyr::mutate(response = factor(
    dplyr::case_when(
      PROGRESSION == 1 ~ "SD/PD",
      is.na(response) ~ RESPONSE_6_M,
      .default = response),
    levels = response_levels)
  ) %>%
  dplyr::select(sample_id, Product, response, start, end)

segment_df = rbind(
  d30, d90, d180, d365, lfu
)

# Add events
death_df = swim_data[,c("sample_id", "Product", "OS", "DEATH")]
colnames(death_df)[3:4] = c("time", "status")
death_df$event = "Death"

progression_df = swim_data[,c("sample_id", "Product", "PFS", "PROGRESSION")]
colnames(progression_df)[3:4] = c("time", "status")
progression_df$event = "Progression"

event_df = rbind(death_df, progression_df) %>%
  dplyr::filter(status == 1) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(difference_pr_de = dplyr::lag(time) - time) %>%
  dplyr::filter(is.na(difference_pr_de) | difference_pr_de != 0) %>%
  ungroup() %>%
  data.frame()

# arrow_df
arrow_df = swim_data %>%
  dplyr::filter(DEATH != 1, PROGRESSION != 1) %>%
  data.frame()

ann_x <- 0.55
ann_y_1 <- 0.2
ann_y_2 <- 0.15

swim_segments_ide =
  ggplot() +
  geom_segment(
    data=subset(segment_df, Product=="ide"),
    aes(x=start, xend=end, y=sample_id, yend=sample_id, color=response),
    alpha=1, linewidth = 2
  ) +
  geom_arrow_segment(
    data = subset(arrow_df,Product=="ide"),
    aes(x = OS+15, xend = OS+62, y = sample_id, colour = BEST_RESPONSE),
    linewidth=.5, alpha=1, show.legend = F
  ) +
  geom_point(
    data = subset(event_df,Product=="ide"),
    aes(y=sample_id, x = time, shape=event), size=1.25
  ) +
  # geom_star(
  #   data = subset(single_cell_tp, Product=="ide"),
  #   aes(y=sample_id, x = TP, fill = "single-cell sample"), size=.75
  # ) +
  scale_fill_manual(values=c("black")) +
  scale_shape_manual(name = "Event", values = c(Progression = 21, Death = 17)) +
  scale_color_manual(name ="IMWG response criteria", values = palette_response) +
  scale_x_continuous(
    breaks = c(-100, -50, seq(0, 900, 200)),
    labels = c(" ", "LA", seq(0, 900, 200)), expand = c(0.02,0.05)
  ) +
  xlab("Days after CAR-T cell infusion") +
  annotate(
    "text", label="Ide-cel (n=34)", hjust=0, size = 3, fontface = "bold",
    x = max(subset(segment_df, Product=="ide")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="ide")$sample_id))*.25
  ) +
  annotate(
    "text", label="ORR: 67%", hjust=0, size = 3,
    x = max(subset(segment_df, Product=="ide")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="ide")$sample_id))*ann_y_1
  ) +
  annotate(
    "text", label="CR rate: 38%", hjust=0, size = 3,
    x = max(subset(segment_df, Product=="ide")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="ide")$sample_id))*ann_y_2
  ) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face="bold", size=rel(1)),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    # legend.title = element_text(margin = margin(b = 2)),
    plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1)),
    panel.background = element_rect(fill = "#F7F5EE"),
    legend.text = element_text(margin = margin(l = 3, unit = "pt"))
  ) +
  guides(
    color = guide_legend(title.position = "left", order=1),
    shape=guide_legend(title.position="left", order=2,  override.aes = list(size = 1.5, stroke = 1)),
    fill=guide_legend(size=3.5, title = "", order=3)
  )
  # ggtitle("Ide-cel (n=34)")

swim_segments_cilta =
  ggplot() +
  geom_segment(
    data=subset(segment_df, Product=="cilta"),
    aes(x=start, xend=end, y=sample_id, yend=sample_id, color=response),
    alpha=1, linewidth = 2.5
  ) +
  geom_arrow_segment(
    data = subset(arrow_df,Product=="cilta"),
    aes(x = OS+15, xend = OS+62, y = sample_id, colour = BEST_RESPONSE),
    linewidth=.5, alpha=1, show.legend = F
  ) +
  geom_point(
    data = subset(event_df,Product=="cilta"),
    aes(y=sample_id, x = time, shape=event), size=1.25
  ) +
  # geom_star(
  #   data = subset(single_cell_tp, Product=="cilta"),
  #   aes(y=sample_id, x = TP, fill = "single-cell sample"), size=.75
  # ) +
  scale_fill_manual(values=c("black")) +
  scale_shape_manual(name = "Event", values = c(Progression = 21, Death = 17)) +
  scale_color_manual(name ="IMWG response criteria", values = palette_response) +
  scale_x_continuous(
    breaks = c(-100, -50, seq(0, 900, 200)),
    labels = c(" ", "LA", seq(0, 900, 200)), expand = c(0.02,0.05)
  ) +
  # facet_wrap(facets = vars(PRODUCT_LAB)) +
  xlab("Days after CAR-T cell infusion") +
  annotate(
    "text", label="Cilta-cel (n=27)", hjust=0, size = 3, fontface = "bold",
    x = max(subset(segment_df, Product=="cilta")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="cilta")$sample_id))*.25
  ) +
  annotate(
    "text", label="ORR: 93%", hjust=0, size = 3,
    x = max(subset(segment_df, Product=="cilta")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="cilta")$sample_id))*ann_y_1
  ) +
  annotate(
    "text", label="CR rate: 78%", hjust=0, size = 3,
    x = max(subset(segment_df, Product=="cilta")$end)*ann_x,
    y = length(unique(subset(segment_df, Product=="cilta")$sample_id))*ann_y_2
  ) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face="bold", size=rel(1)),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1)),
    panel.background = element_rect(fill = "#F7F5EE"),
    legend.text = element_text(margin = margin(l = 3, unit = "pt"))

  ) +
  guides(
    color = guide_legend(title.position = "left", order=1),
    shape=guide_legend(title.position="left", order=2,  override.aes = list(size = 1.5, stroke = 1)),
    fill=guide_legend(size=3.5, title = "", order=3)
  )
  # ggtitle("Cilta-cel (n=27)")

swim_segments =
  plot_grid(
  plot_grid(
    plot_grid(NULL, swim_segments_ide + theme(legend.position = "none"), nrow = 2, rel_heights = c(.05, 1)),
    NULL,
    plot_grid(NULL, swim_segments_cilta + theme(legend.position = "none"), nrow = 2, rel_heights = c(.05, 1)),
    nrow=1, rel_widths = c(1, 0.125, 1),
    labels = c("B", "", "C"), label_fontface = "bold",
    label_size = 12, vjust = 1.1
  ),
  cowplot::get_legend(swim_segments_ide),
  nrow=2, rel_heights = c(6,.5)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1d-f ) Kaplan Meier Analyses
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
surv_plot_product = surv_plot(
  obj = clin_data,
  group = "PRODUCT",
  legend.labs = c("cilta" = "Cilta-cel", "ide" = "Ide-cel"),
  leg.title = "CAR T cell product",
  col.pal = c("#006EAE", "#647314")
)

clin_data$PRODUCT_C = factor(clin_data$PRODUCT, levels = c("ide", "cilta"))
product_hr_table = cox_table(
  clin_data,
  "PRODUCT_C",
  "Product (cilta-cel)"
)

survfit(Surv(PFS_M, STATUS) ~ BEST_RESPONSE_CONSENSUS, clin_data)

surv_plot_response = surv_plot(
  obj = clin_data,
  group = "BEST_RESPONSE_CONSENSUS",
  legend.labs = c("CR" = "CR", "non-CR" = "non-CR"),
  leg.title = "Response"
)

clin_data$BEST_RESPONSE_CONSENSUS_CR = factor(
  clin_data$BEST_RESPONSE_CONSENSUS, levels = c( "non-CR", "CR")
)
response_hr_table = cox_table(
  clin_data,
  "BEST_RESPONSE_CONSENSUS_CR",
  "Response (CR)"
)

survfit(Surv(PFS_M, STATUS)~REFRACTORINESS_BIN, clin_data)

table(clin_data$REFRACTORINESS_BIN)
surv_plot_refr = surv_plot(
  obj = clin_data,
  group = "REFRACTORINESS_BIN",
  legend.labs = c("TCexposed" = "TCexposed", "TCRRMM/PentaRRMM" = "TCRRMM/PentaRRMM"),
  leg.title = "Refractoriness",
  col.pal = c("#000000", "#969696")
)

refr_hr_table = cox_table(
  clin_data,
  "REFRACTORINESS_BIN",
  "TCRRMM/PentaRRMM"
)

pl_2 =
  plot_grid(
  plot_grid(
    surv_plot_product,
    NULL,
    surv_plot_response,
    NULL,
    surv_plot_refr,
    nrow = 1, rel_widths = c(1,0.1,1,0.1,1),
    labels = c("D", "", "E", "", "F"), label_fontface = "bold", label_size = 12, vjust = 2
  ),
  plot_grid(
    NULL,
    product_hr_table,
    NULL,
    response_hr_table,
    NULL,
    refr_hr_table,
    NULL,
    nrow = 1, rel_widths = c(.07, 1, 0.15, 1, 0.1, 1, .0)
  ),
  nrow=2, rel_heights = c(3,1)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1g) COX-PH model adjusted for product:
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sbcma_la = elisa_data %>%
  dplyr::filter(DAY == "LA") %>%
  dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
  dplyr::rename(sbcma_la = "BCMA_NG_ML")

sbcma_day_0 = elisa_data %>%
  dplyr::filter(DAY == "Day 0") %>%
  dplyr::select(SAMPLE_ID, BCMA_NG_ML) %>%
  dplyr::rename(sbcma_day0 = "BCMA_NG_ML")

median_sbcma_la = median(
  as.numeric(subset(elisa_data, DAY == "Day 0")$BCMA_NG_ML), na.rm = T
)

max_exp = imm.s %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::slice_max(CD3_CAR_PERC, with_ties = F) %>%
  dplyr::rename(max_expansion = "CD3_CAR_PERC") %>%
  dplyr::select(SAMPLE_ID, max_expansion)

clin_data = clin_data %>%
  dplyr::left_join(sbcma_day_0) %>%
  dplyr::left_join(max_exp) %>%
  dplyr::mutate(
    sbcma_day0_bin = factor(ifelse(sbcma_day0 > median_sbcma_la, "High", "Low"), levels = c("Low", "High")),
    status_before_cart_num = rev(as.numeric(STATUS_BEFORE_CAR))
  )

covariates <- c(
  "status_before_cart_num", "EMD_BEFORE_CART", "HR_CYTOGEN", "sbcma_day0_bin"
)
clin_data[[covariates[2]]] = factor(clin_data[[covariates[2]]])

cox.model <- sapply(covariates, function(x) {
  res = coxph(
    as.formula(paste('Surv(PFS, STATUS)~', x, "+ PRODUCT")),
    data = clin_data
  )
  res = summary(res)
  res = c(
    "logHR" = signif(res$coef[1, , drop = F][1], digits = 2),
    "SE_logHR" = res$coef[1, , drop = F][3],
    "L90CI" = res$coef[1, 1, drop = T] - res$coef[1, 3, drop = T] * 1.65,
    "U90CI" = res$coef[1, 1, drop = T] + res$coef[1, 3, drop = T] * 1.65,
    "Pval" = signif(res$coef[1, , drop = F][5], digits = 2)
  )
  res
}) %>% t()

cox.model = as.data.frame(cox.model, check.names = FALSE) %>%
  tibble::rownames_to_column("category") %>%
  dplyr::arrange(logHR)
cox.model$category = factor(cox.model$category, levels = cox.model$category)
cox.model$Pval = ifelse(cox.model$Pval < 0.1, paste0("p=", cox.model$Pval), "")

y.labels = setNames(
  c(
    "Remission\nbefore CAR T", "EMD before\ninfusion",
    "High-risk\ncytogenetics", paste0("sBCMA day 0\n(>",round(median_sbcma_la,1)," ng/mL)")
  ),
  covariates
)

cox.model.HR = cox.model
cox.model.HR$logHR = exp(cox.model.HR$logHR)
cox.model.HR$L90CI = exp(cox.model.HR$L90CI)
cox.model.HR$U90CI = exp(cox.model.HR$U90CI)

forest.pl =
ggplot(
  cox.model.HR,
  aes(x=logHR, y=(category), xmin=L90CI, xmax=U90CI, color = logHR)
) +
  geom_vline(
    aes(xintercept=1), color="grey50", linetype="dashed", alpha=.6, lwd = .3
  ) +
  geom_linerange(linewidth=2.25) +
  geom_point(size=3.5, shape=16, stroke = 1) +
  geom_point(size=4, shape=21, colour="white", stroke = 1) +
  geom_text(
    aes(x = logHR, label = Pval), color = "black",
    hjust = .5, vjust = -1.3, size = 2.5
  ) +
  theme(
    #aspect.ratio = .7,
    legend.position = "none",
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = rel(1)),
  ) +
  scale_color_scico(palette = "vikO", midpoint = 1, begin = .2, end = .8) +
  xlab("Hazard ratio (90% CI)") +
  scale_x_log10() +
  scale_y_discrete(labels = y.labels) +
  ggtitle("Cox Regression")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Remission status at CAR-T infusion
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rem.pl = clin_data %>%
  dplyr::group_by(PRODUCT, STATUS_BEFORE_CAR, BEST_RESPONSE) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::mutate(perc = n/sum(n))
rem.pl$PRODUCT = stringr::str_to_title(rem.pl$PRODUCT)
rem.pl = ggplot(rem.pl, aes(x=STATUS_BEFORE_CAR, fill=BEST_RESPONSE, y=perc)) +
  geom_bar(position = "fill", stat="identity") +
  geom_text(aes(label=n), position = position_stack( vjust=.5), size = 3, show.legend = F) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.margin = margin(t=-5),
    panel.spacing = unit(.5, "lines"),
    axis.title.y = element_text(vjust = + 2)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  facet_grid(. ~ PRODUCT) +
  xlab("Remission before CAR T") +
  ylab("Proportion") +
  scale_fill_manual(name = "Response\nafter\nCAR-T\ninfusion", values = palette_response)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1h-i) CAR expansion
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plot_expansion = function(
    obj.imm = imm.s,
    obj.pd = clin_data,
    group = "PRODUCT",
    palette = c(cilta = "#004488", ide = "#44AA99"),
    labels = c("cilta-cel", "ide-cel"),
    leg.pos = "bottom"
) {

  pdata.exp = obj.pd
  pdata.exp$GROUP = pdata.exp[[group]]

  imm.exp = obj.imm %>%
    dplyr::select(SAMPLE_ID, DAY, CD3_CAR_PERC, CAR_RATIO) %>%
    dplyr::filter(DAY != "Leukapheresis") %>%
    dplyr::left_join(pdata.exp[,c("SAMPLE_ID", "GROUP")]) %>%
    dplyr::mutate(CAR_RATIO = ifelse(DAY == "Day 0", 0, CAR_RATIO)) %>%
    tidyr::pivot_longer(c("CD3_CAR_PERC", "CAR_RATIO")) %>%
    dplyr::mutate(
      DAY = factor(DAY, levels = c("Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))
    )

  p.wlx = imm.exp %>%
    dplyr::filter(DAY != "Day 0") %>%
    dplyr::group_by(DAY, name) %>%
    rstatix::wilcox_test(value ~ GROUP) %>%
    rstatix::add_significance(p.col = "p", output.col = "p_signif")
    p.wlx$p_signif[p.wlx$p_signif == "ns"] = ""

  l = list()
  for (i in 1:length(p.wlx$p)) {
    if (p.wlx$p[i] <= 0.0001) {
      l[[i]] = formatC(p.wlx$p[i], format = "e", digits= 2)
    } else if (p.wlx$p[i] > 0.0001 & p.wlx$p[i] <= 0.001) {
      l[[i]] = format(round(p.wlx$p[i], digits=5), nsmall = 5)
    } else {
      l[[i]] = format(round(p.wlx$p[i], digits=3), nsmall = 3)
    }
  }
  p.wlx$p_label = unlist(l)
  p.wlx$p_label = paste0("p =\n", p.wlx$p_label)
  p.wlx$p_label = ifelse(p.wlx$p < .1, p.wlx$p_label, "")

  a = imm.exp %>% dplyr::filter(name == "CD3_CAR_PERC", !is.na(value))
  tmp = ggplot(a, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess")
  max.v = max(ggplot_build(tmp)[[1]][[1]]$ymax)
  min.v = min(ggplot_build(tmp)[[1]][[1]]$ymin)

  pl.perc = ggplot(a, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_rect(aes(xmin=1, xmax=2, ymin=-Inf, ymax = Inf), alpha=0.025, fill="#DDDDDD", color = NA) +
    geom_smooth(aes(group=GROUP), linewidth=.75, method = "loess") +
    coord_cartesian(ylim=c(min.v, max.v+ 5), xlim = c(1.35, 4.75)) +
    geom_text(
      data = subset(p.wlx, name == "CD3_CAR_PERC"),
      aes(x=DAY, y=Inf, label = p_label), inherit.aes=F, size=2.5, vjust=1.25, hjust=.5
    ) +
    scale_color_manual(name = " ", values = palette, labels = labels) +
    scale_fill_manual(name = " ", values = palette, labels = labels) +
    ylab("% CD3+CAR+") +
    theme(
      legend.position = "bottom",
      axis.title.y = element_text(vjust = + 2),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      legend.key.spacing.x = unit(10, "pt")
    ) +
    ggtitle("CD3+CAR+")

  ###

  b = imm.exp %>% dplyr::filter(name == "CAR_RATIO", !is.na(value))
  tmp = ggplot(b, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_smooth(aes(group=GROUP), linewidth=.5, method = "loess")
  max.v = max(ggplot_build(tmp)[[1]][[1]]$ymax)
  min.v = min(ggplot_build(tmp)[[1]][[1]]$ymin)

  pl.ratio =
    ggplot(b, aes(x=DAY, y=value, color=GROUP, fill=GROUP)) +
    geom_rect(aes(xmin=1, xmax=2, ymin=-Inf, ymax = Inf), alpha=0.025, fill="#DDDDDD", color = NA) +
    geom_smooth(aes(group=GROUP), linewidth=.75, method = "loess") +
    coord_cartesian(ylim=c(min.v, max.v + .1), xlim = c(1.35, 4.75)) +
    geom_text(
      data = subset(p.wlx, name == "CAR_RATIO"),
      aes(x=DAY, y=Inf, label = p_label), inherit.aes=F, size=2.5, vjust=1.25, hjust=.5
    ) +
    scale_color_manual(name = " ", values = palette, labels = labels) +
    scale_fill_manual(name = " ", values = palette, labels = labels) +
    ylab("CD4+CAR+/CD8+CAR+") +
    theme(
      # aspect.ratio = 1,
      legend.position = "bottom",
      axis.title.y = element_text(vjust = + 2),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, hjust=1, vjust=1),
      plot.title = element_text(hjust = 0.5, size = rel(1)),
      legend.key.spacing.x = unit(10, "pt")
    ) +
    ggtitle("CAR ratio")

  plot_grid(
    plot_grid(
      pl.perc + theme(legend.position = "none"),
      NULL,
      pl.ratio + theme(legend.position = "none"),
      ncol = 3, rel_widths = c(1, .1, 1),
      labels = c("I", "", "J"), label_fontface = "bold", label_size = 12,
      vjust = 1.1
    ),
    plot_grid(ggpubr::get_legend(pl.perc)),
    nrow = 2, rel_heights = c(1, .1)
  )
}

pl.exp.response = plot_expansion(
  group = "BEST_RESPONSE_CONSENSUS",
  palette = c("CR"="#7B9AB6", "non-CR"="#9B740A"),
  labels = c("non-CR" = "non-CR", "CR"="CR"),
  leg.pos = "bottom"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final pot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

o.pl = ggdraw() +
  cowplot::draw_image(
    "publication/figures_main_new/fig_01_overview.PNG", scale = .855)

ggsave2(
  "publication/figures_main_new/fig_01.png",
  plot_grid(
    plot_grid(
      o.pl + theme(plot.margin = unit(c(0,0,28,-30), "pt")),
      swim_segments, NULL,
      nrow = 1, rel_widths = c(1.22, 1.5, .01),
      labels = c("A", "", ""), label_fontface = "bold", label_size = 12, vjust = 1.1
    ),
    NULL,
    pl_2,
    NULL,
    plot_grid(
      rem.pl,
      NULL,
      plot_grid(forest.pl, NULL, nrow = 2, rel_heights = c(1, .15)),
      NULL,
      pl.exp.response,
      NULL,
      ncol = 6, rel_widths = c(.3, .01, .245, .025, .43, 0),
      labels = c("G", "", "H", "", ""), label_fontface = "bold", label_size = 12, vjust = 1.1
    ),
    nrow = 5, rel_heights = c(.374, -.01, .39, .03, .29)
  ),
  height = 165,
  width = 165,
  units="mm",
  dpi=400,
  scale=1.5,
  device = png, type = "cairo",
  bg= "white"
)


ggsave2(
  "publication/figures_main_new/fig_01.pdf",
  plot_grid(
    plot_grid(
      o.pl + theme(plot.margin = unit(c(0,0,28,-30), "pt")),
      swim_segments, NULL,
      nrow = 1, rel_widths = c(1.22, 1.5, .01),
      labels = c("A", "", ""), label_fontface = "bold", label_size = 12, vjust = 1.1
    ),
    NULL,
    pl_2,
    NULL,
    plot_grid(
      rem.pl,
      NULL,
      plot_grid(forest.pl, NULL, nrow = 2, rel_heights = c(1, .15)),
      NULL,
      pl.exp.response,
      NULL,
      ncol = 6, rel_widths = c(.3, .01, .245, .025, .43, 0),
      labels = c("G", "", "H", "", ""), label_fontface = "bold", label_size = 12, vjust = 1.1
    ),
    nrow = 5, rel_heights = c(.374, -.01, .39, .03, .29)
  ),
  height = 165,
  width = 165,
  units="mm",
  dpi=400,
  scale=1.5,
  bg= "white",
  device = cairo_pdf
)


