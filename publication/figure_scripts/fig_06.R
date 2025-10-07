# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Install and load packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "scCustomize","harmony","beeswarm","ggplot2","reshape2","tibble",
  "ggpubr","cowplot","tidyr","dplyr","rstatix","ggtext","conflicted","UCell",
  "scico","GGally","stringr","ggcorrplot","scRepertoire","yaml","tidycmprsk","ggsurvfit",
  "survival","parameters","forcats", "readxl"
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

if (any(!"Startrac" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("Japrin/sscVis")
  devtools::install_github("ncborcherding/Startrac")
}
library(Startrac)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load functions and styles
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("code/helper/functions.R")
source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions_plots_bcell.R")
theme_set(mytheme())

cust.th = theme(
  legend.spacing.y = unit(2, 'mm'),
  legend.key.size = unit(4, "mm"),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  plot.title = element_blank()
)

base.size = 15

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read datasets
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
bcell_df = readRDS("./data/bcell_df.RDS")

se.meta = readRDS(paste0(manifest$meta_pub, "integration/05_seurat_harmony_all_new.Rds"))
clono_vdj_b = readRDS(paste0(manifest$meta_pub, "integration/05_vdj_b_new.Rds"))
se.meta = Seurat::AddMetaData(se.meta, clono_vdj_b)

pdata = readRDS("publication/clinicial_data/clinical_table_DF_2024_10_28.Rds")

if(is.null(se.meta@meta.data$CLONE_PSEUDO_ID)){
  pd = se.meta@meta.data
  pd$barcode = rownames(pd)
  df = pd %>%
    dplyr::group_by(CTstrict) %>%
    dplyr::count(name = "cloneFreqAll") %>%
    dplyr::arrange(desc(cloneFreqAll))
  df$pseudo_id = paste0("Clone_", seq(1:nrow(df)))
  se.meta$CLONE_PSEUDO_ID = df$pseudo_id[match(se.meta$CTstrict, df$CTstrict)]
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# B-cell clones only B- and Plasma cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.bcell = subset(se.meta, subset = celltype_short_3 == "B-Cell")
se.plasmacell = subset(se.meta, subset = celltype_short_3 == "Plasma cell")
se.bcell_comb = subset(se.meta, subset = celltype_short_3 %in% c("B-Cell","Plasma cell"))

### Integration Analysis of B-cell subsets
dims.use.rna = 15

# RNA Integration
se.bcell_comb = integration(
  obj = se.bcell_comb,
  no.ftrs = 1000,
  threads = 1,
  .nbr.dims = dims.use.rna,
  run.integration = T,
  harmony.group.vars = c("STUDY")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc colored by celltype
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd = get_metadata(se.bcell_comb)
pd$celltype=ordered(pd$celltype,levels = c(
  "B naive",
  "B naive ISG+",
  "B Transitional",
  "B effector CD27-",
  "B effector CD27+",
  "B memory",
  "B memory CD95+",
  "Plasma cell"
))

b_cell.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = celltype)) +
  scattermore::geom_scattermore(pointsize = 8, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 2.5, pixels = c(512,512)) +
  mytheme(base_size = base.size) +
  theme(
    legend.position = "right",
    legend.margin = margin(t=-5, b = 9),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 6)
  )) +
   scale_color_manual(
     values = bcell.col
   )
#b_cell.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc colored by clonesize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd$cloneSize = as.character(pd$cloneSize)
pd$cloneSize = ifelse(
  grepl("Hyper", pd$cloneSize), "Hyperexpanded (100 > X)", pd$cloneSize
)
pd$cloneSize = gsub(" \\(", "\\\n(", pd$cloneSize)
pd$cloneSize = forcats::fct_relevel(
  as.factor(pd$cloneSize), "Single\n(0 < X <= 1)", after = Inf
)

clono.col = setNames(
  c("#003285", "#5AB2FF", "#A0DEFF", "#F08A5D", "#B83B5E", "#BBBBBB"),
  c(rev(levels(pd$cloneSize)), "NA")
)

clonesz.pl =
  ggplot(data = pd, aes(x = umap_1, y = umap_2, col = cloneSize)) +
  scattermore::geom_scattermore(pointsize = 8, color="black")+
  scattermore::geom_scattermore(pointsize = 6, color="white") +
  scattermore::geom_scattermore(pointsize = 2.5, pixels = c(512,512)) +
  mytheme(base_size = base.size) +
  theme(
    legend.position = "right",
    legend.margin = margin(t=-5, b = 9),
    legend.title = element_text(margin = margin(r = 3)),
    legend.spacing.y = unit(2, 'mm'),
    legend.key.size = unit(4, "mm"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_blank()
  ) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 6)
  ))+
  scale_color_manual(values = clono.col)
# clonesz.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Candidate genes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ftrs_to_plot=c(
  "FCRL5","CD38","GPRC5D","LAMP5","TNFRSF17"
)
names(ftrs_to_plot) = ftrs_to_plot
#names(ftrs_to_plot)[names(ftrs_to_plot) == "TNFRSF17"] = "BCMA"

global.max = quantile(unlist(FetchData(se.bcell_comb, ftrs_to_plot)), .99)

ftrs.l = list()
pd = get_metadata(se.bcell_comb)

for (i in 1:length(ftrs_to_plot)) {
  exprs = FetchData(se.bcell_comb, vars = ftrs_to_plot[i], layer = "data")
  pd$EXPRS = exprs[[1]]
  pd$EXPRS[pd$EXPRS >= global.max] = global.max
  c = c(scico::scico(30, palette = "navia", direction = -1))

  pl =
    ggplot(pd, aes(x = umap_1, y = umap_2, color = EXPRS)) +
    scattermore::geom_scattermore(pointsize = 6.5, color="black")+
    scattermore::geom_scattermore(pointsize = 4, color="white") +
    scattermore::geom_scattermore(pointsize = 3, pixels = c(512,512)) +
    scale_color_gradientn(
      colors =c[1:(length(c)-4)],
      na.value = "white", limits = c(0, global.max),
      breaks = scales::pretty_breaks(3)
    ) +
    mytheme(base_size = base.size) +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title  = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, colour = "black", size = rel(1))
    ) +
    ggtitle(names(ftrs_to_plot[i])) +
    guides(
      color = guide_colorbar(
        title = "RNA Exprs ", title.vjust = 1.1,
        barwidth = unit(7.5, 'lines'), barheight = unit(.6, 'lines'),
        ticks.linewidth = 2/.pt, ticks = T, frame.colour="black",
        frame.linewidth = 0.5/.pt
      )
    )
  ftrs.l[[i]] = pl
}

legend = ggpubr::get_legend(ftrs.l[[1]],)
ftrs.l = lapply(ftrs.l, function(x){
  x = x + theme(legend.position='none', panel.border = element_blank())
})

ftrs.pl = plot_grid(plotlist = ftrs.l, nrow = 2, scale = .98)
ftrs.pl = plot_grid(ftrs.pl, ggdraw(legend), rel_heights = c(1, .15), nrow = 2)
#ftrs.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# B-cells and plasma cells in patients
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.bcell_comb$PATIENT_ID_SHORT = gsub("Patient0", "P", se.bcell_comb$SAMPLE_ID)
se.bcell_comb$PATIENT_ID_SHORT = gsub("_1$", "", se.bcell_comb$PATIENT_ID_SHORT)
se.bcell_comb@meta.data = droplevels(se.bcell_comb@meta.data)

df = table(se.bcell_comb$celltype,se.bcell_comb$PATIENT_ID_SHORT) %>%
  reshape2::melt(varnames = c("CellType", "Patient"), value.name = "CellCount")

df <- df %>%
  group_by(Patient) %>%
  mutate(TotalCells = sum(CellCount)) %>%
  ungroup()

df <- df %>%
  mutate(Patient = reorder(Patient, -TotalCells))

df <- df %>%
  mutate(CellType = ordered(CellType, levels =
                              c(
                                "B naive",
                                "B naive ISG+",
                                "B Transitional",
                                "B effector CD27-",
                                "B effector CD27+",
                                "B memory",
                                "B memory CD95+",
                                "Plasma cell"
                              )))


levels(df$Patient)[df$Patient[df$CellType=="Plasma cell" & df$CellCount>=50]] = paste0(
  "* ",levels(df$Patient)[df$Patient[df$CellType=="Plasma cell" & df$CellCount>=50]]
)
df = df[df$TotalCells>=100,]

bar.pl = ggplot(df, aes(x = Patient, y = CellCount, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Cells of B-cell lineage per Patient",
    x = "Patients with \u2265 100 cells",
    y = "Number of Cells",
    fill = "Cell Type"
  ) +
  mytheme(base_size = base.size) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(
    values = bcell.col
  ) +
  theme(
    legend.position = "inside",
    legend.title = element_blank(),
    legend.text = element_text(size = rel(.9)),
    legend.position.inside = c(0.8, 0.55)
  )
#bar.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Shannon diversity B-cells vs plasma cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##############################################################
# B-cell analysis
# Clone Pseudo IDs
se.bcell_comb_sh =subset(se.bcell, subset = cloneSize %in% levels(se.bcell$cloneSize))

star.lp.bc = startrac(obj = se.bcell_comb_sh, tp = "LP", celltype = "celltype_short_3")
star.l.bc = startrac(obj = se.bcell_comb_sh, tp = "Late", celltype = "celltype_short_3")
star.vl.bc = startrac(obj = se.bcell_comb_sh, tp = "Very Late", celltype = "celltype_short_3")

##############################################################
# Plasma cell analysis
# Clone Pseudo IDs
se.bcell_comb_sh =subset(se.plasmacell, subset = cloneSize %in% levels(se.plasmacell$cloneSize))

star.lp.pc = startrac(obj = se.bcell_comb_sh, tp = "LP", celltype = "celltype_short_3")
star.l.pc = startrac(obj = se.bcell_comb_sh, tp = "Late", celltype = "celltype_short_3")
star.vl.pc = startrac(obj = se.bcell_comb_sh, tp = "Very Late", celltype = "celltype_short_3")

star = rbind(star.lp.bc, star.l.bc, star.vl.bc,star.lp.pc, star.l.pc, star.vl.pc)

star = merge(
  star, se.bcell_comb_sh@meta.data[!duplicated(se.bcell_comb_sh$orig.ident), ],
  by.x = "aid", by.y = "orig.ident"
)

shannon.pl =
  ggplot(star, aes(majorCluster, value, fill = majorCluster)) +
  geom_boxplot(outlier.colour = NA, position=position_dodge(0.8), fatten = 1, alpha = .5, lwd = 1) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), stroke = .5, dotsize = 1.1) +
  # ylim(0, max(pd.cl$freq + .1)) +
  stat_compare_means(
    aes(group = majorCluster, label = paste0("p = ",..p.format..)), label.x = 1.5, label.y = 1.1, hjust = .5, vjust = 1,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#CCBB44","#7C2529")) +
  scale_x_discrete(labels = c("B-cell"="B-cell", "Plasma cell"="PC")) +
  xlab(NULL) + ylab("Clonality") +
  mytheme(base_size = base.size) +
  theme(
    panel.spacing = unit(.5, "lines"),
    axis.title.y = element_text(vjust = + 2),
    legend.margin=margin(t=-5),
    legend.box.spacing = unit(0, "pt"),
    legend.position = "none"
  ) +
  facet_wrap( ~ TIMEPOINT)
#shannon.pl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dynamics BMCA cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# B-cells
bcell.dyn.pl =plot_longitudinal_cells(c("B naive",
  "B naive ISG+",
  "B Transitional",
  "B effector CD27-",
  "B effector CD27+",
  "B memory",
  "B memory CD95+"), test = F) +
  theme(
    legend.margin = margin(b = -12),
    legend.key.size = unit(7, "mm"),
    legend.position = "bottom"
  )
bcell.dyn.pl = bcell.dyn.pl + stat_compare_means(
  aes(label = paste("p =",after_stat(p.format))),
  label.x = 1.5, hjust = .5, vjust = 1, size = 4.5
)

#pDC
pdc.dyn.pl = plot_longitudinal_cells(c("pDC"),label = "pDCs in %",test = F)
pdc.dyn.pl = pdc.dyn.pl +
  theme(
    legend.margin = margin(b = -12),
    legend.key.size = unit(7, "mm"),
    legend.position = "bottom"
  ) +
  stat_compare_means(
    aes(label = paste("p =",after_stat(p.format))),
    label.x = 1.5, hjust = .5, vjust = 1, size = 4.5
  )

# Plasma cells
plasma.dyn.pl = plot_longitudinal_cells(c("Plasma cell"),label = "Plasma cells in %",test = F)
plasma.dyn.pl = plasma.dyn.pl +
  theme(
    legend.margin = margin(b = -12),
    legend.key.size = unit(7, "mm"),
    legend.position = "bottom"
  ) +
  scale_y_sqrt() +
  stat_compare_means(
    aes(label = paste("p =",after_stat(p.format))),
    label.x = 1.5, hjust = .5, vjust = 1, size = 4.5
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dynamics clinical B-cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bcell_df_form = bcell_df %>%
  mutate(cumulative_days = floor(as.numeric(difftime(Datum, date_cart, units = "days")))) %>%
  mutate(diff_lp = floor(as.numeric(difftime(Datum, date_leukapheresis, units = "days")))) %>%
  rstatix::filter(!is.na("date_cart"), cumulative_days<1000, diff_lp>-10)
bcell_df_form$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(bcell_df_form$patient_id,se.meta$PATIENT_ID)]
bcell_df_form$PRODUCT = se.meta$PRODUCT[match(bcell_df_form$patient_id,se.meta$PATIENT_ID)]
bcell_df_form = bcell_df_form[!is.na(bcell_df_form$BEST_RESPONSE_CONSENSUS),]

tps = c(seq(-10,200,10))
bcell_df_inter = data.frame(patient_id=character(), Wert=numeric(), cumulative_days=numeric(), BEST_RESPONSE_CONSENSUS=character(), PRODUCT=character())
for (i in unique(bcell_df_form$patient_id)){
  sel = bcell_df_form$patient_id == i
  interpol = approx(x=bcell_df_form$cumulative_days[sel],y=bcell_df_form$Wert[sel],xout=tps,ties = mean)
  df = data.frame(
    patient_id=rep(i,length(interpol[[1]])),
    Wert=interpol$y,
    cumulative_days=interpol$x,
    BEST_RESPONSE_CONSENSUS=bcell_df_form$BEST_RESPONSE_CONSENSUS[bcell_df_form$patient_id == i][[1]],
    PRODUCT=bcell_df_form$PRODUCT[bcell_df_form$patient_id == i][[1]]
  )
  bcell_df_inter = rbind(bcell_df_inter,df)
}

p.vals = bcell_df_inter %>%
  dplyr::group_by(cumulative_days) %>%
  pairwise_wilcox_test(Wert ~ BEST_RESPONSE_CONSENSUS, alternative = "two.sided")
#print(p.vals)
p.vals$y.position = max(bcell_df_inter$Wert,na.rm = T)*1.05
p.vals$xmin=tps
p.vals$xmax=tps

mean_bcell_df_inter <- bcell_df_inter %>%
  dplyr::group_by(cumulative_days, BEST_RESPONSE_CONSENSUS) %>%
  dplyr::summarise(mean_value = mean(Wert,na.rm=T), .groups = "drop")

bcell.clin.dyn.pl = bcell_df_inter  %>% drop_na(Wert) %>%
  ggplot(aes(x=cumulative_days, y=Wert,color=BEST_RESPONSE_CONSENSUS)) +
  geom_line(aes(group = patient_id),alpha=0.3) +
  scale_color_manual(values = c("#7B9AB6", "#9B740A"), name="") +
  geom_point(alpha=0.3)+
  stat_pvalue_manual(p.vals,inherit.aes = F) +
  geom_line(data = mean_bcell_df_inter, aes(x = cumulative_days, y = mean_value, group = BEST_RESPONSE_CONSENSUS),
            linewidth = 1.5, linetype = "solid") +
  mytheme_grid(base_size = 15) +
  geom_hline(aes(yintercept=109), color="darkblue", linetype="dashed", alpha=.6) +
  geom_hline(aes(yintercept=411), color="firebrick", linetype="dashed", alpha=.6) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_markdown(),
    legend.position = "none",
    legend.title = element_text()
  ) +
  labs(
    x="Days after CAR-T cell infusion",
    y="B cells/&mu;l"
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Longitudinal_sBCMA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pdata.elisa = pdata$pdata.elisa

df.cell = pdata.elisa
ftr = table(df.cell$DAY)
ftr = names(ftr)[ftr>20]
df.cell = df.cell[df.cell$DAY %in% ftr,]
df.cell$DAY = ordered(df.cell$DAY,levels=c("LA","Day 0","Day 7","Day 30","Day 100"))
df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]
df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]

df.cell.day = df.cell[!is.na(df.cell$BEST_RESPONSE_CONSENSUS) & !is.na(df.cell$BCMA_NG_ML),]
df.cell.day$DAY = gsub("Day ","",df.cell.day$DAY)
df.cell.day$DAY[df.cell.day$DAY=="LA"] = "-5"
df.cell.day$DAY = as.numeric(df.cell.day$DAY)

p.vals = df.cell %>% group_by(DAY) %>% pairwise_wilcox_test(BCMA_NG_ML ~ BEST_RESPONSE_CONSENSUS, alternative = "two.sided")
print(p.vals)
p.vals$xmin = sort(unique(df.cell.day$DAY))
p.vals$xmax = sort(unique(df.cell.day$DAY))
p.vals$y.position = log10(max(df.cell$BCMA_NG_ML,na.rm = T)*1.15)
p.vals$BEST_RESPONSE_CONSENSUS = "CR"

options(scipen=6)
bcma.dyn.pl = ggplot(
  data = df.cell.day,
  aes(x=DAY, y=BCMA_NG_ML, fill=BEST_RESPONSE_CONSENSUS,col=BEST_RESPONSE_CONSENSUS, group = BEST_RESPONSE_CONSENSUS)
  ) +
  geom_point() +
  geom_smooth(method="loess", formula = 'y ~ x')+
  xlab("Days after CAR-T infusion") +
  ylab("sBCMA [ng/mL]")+
  mytheme(base_size = base.size) +
  scale_fill_manual(values = c("#7B9AB6", "#9B740A")) +
  scale_color_manual(values = c("#7B9AB6", "#9B740A")) +
  stat_pvalue_manual(data =p.vals, size = 5) +
  scale_y_log10(limits = c(0.01, 10000)) +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15,0.15)
  )
#bcma.dyn.pl


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Interactions BCMA vs BCMA-celltypes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = data.frame(PATIENT_ID=se.meta$PATIENT_ID,TIMEPOINT=se.meta$TIMEPOINT,celltype=se.meta$celltype)
df=df[df$TIMEPOINT!="pre_Tec",]

df = group_by(df,PATIENT_ID, TIMEPOINT, celltype) %>%
  summarise(CellCount = n(), .groups = 'drop')

df.cell = df %>% group_by(PATIENT_ID, TIMEPOINT) %>%
  summarise(
    TotalCells = sum(CellCount),
    ct_bcell = sum(CellCount[celltype %in% c("B naive","B naive ISG+","B Transitional","B effector CD27-", "B effector CD27+", "B memory", "B memory CD95+")]),
    ct_pDC = sum(CellCount[celltype %in% "pDC"]),
    ct_plasma = sum(CellCount[celltype %in% "Plasma cell"]),
    .groups = 'drop'
  )
df.cell$perc_bcell = df.cell$ct_bcell/df.cell$TotalCells*100
df.cell$perc_pDC = df.cell$ct_pDC/df.cell$TotalCells*100
df.cell$perc_plasma = df.cell$ct_plasma/df.cell$TotalCells*100


df.cell$BEST_RESPONSE_CONSENSUS = se.meta$BEST_RESPONSE_CONSENSUS[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]
df.cell$PRODUCT = se.meta$PRODUCT[match(df.cell$PATIENT_ID,se.meta$PATIENT_ID)]

pdata.elisa.merge = pdata.elisa
pdata.elisa.merge$Day[pdata.elisa.merge$DAY=="LA"] = "LP"
pdata.elisa.merge$Day[pdata.elisa.merge$DAY=="Day 30"] = "Late"
pdata.elisa.merge$Day[pdata.elisa.merge$DAY=="Day 100"] = "Very Late"
pdata.elisa.merge$mtc = paste(pdata.elisa.merge$PATIENT_ID, pdata.elisa.merge$Day,sep="_")
df.cell$mtc = paste(df.cell$PATIENT_ID,df.cell$TIMEPOINT,sep="_")

mtc = match(df.cell$mtc, pdata.elisa.merge$mtc)
df.cell$BCMA = pdata.elisa.merge$BCMA_NG_ML[mtc]
df.cell = df.cell %>% droplevels()

# LP
tp = "LP"
df = df.cell %>% subset(TIMEPOINT == tp) %>% dplyr::select(perc_bcell,perc_pDC,perc_plasma,BCMA)
cor.lp.pt = plot_corr_ct(df,lab=tp,leg = FALSE)

# Late
tp = "Late"
df = df.cell %>% subset(TIMEPOINT == tp) %>% dplyr::select(perc_bcell,perc_pDC,perc_plasma,BCMA)
cor.late.pt = plot_corr_ct(df,lab=tp,leg = FALSE)

# Very Late
tp = "Very Late"
df = df.cell %>% subset(TIMEPOINT == tp) %>% dplyr::select(perc_bcell,perc_pDC,perc_plasma,BCMA)
cor.verylate.pt = plot_corr_ct(df,lab=tp,leg = FALSE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create correlations plots for sBCMA change vs. Expansion, CRP, and CRS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sbcma_change = get_sBCMA_data(pdata)

bcma.exp.pt = sbcma_change %>% drop_na(change, CD3_CAR_PERC) %>%
  ggplot(aes(x=log(change,2) , y = CD3_CAR_PERC)) +
  geom_point(aes(fill=CRS_GRADE), shape=21, size = 3) +
  geom_smooth(method="lm", formula = 'y ~ x', color="grey20", linewidth = .3) +
  stat_cor(method = "spearman", size=5.5) +
  guides(fill = guide_legend(title = "CRS grade")) +
  ylab("Max CAR+/CD3+") + xlab("logFC sBCMA (day30 vs day0)")  +
  mytheme(base_size = base.size) +
  # theme(aspect.ratio = .618) +
  scale_fill_manual(
    values = c("#6699CC", "#EE99AA", "#994455")
  ) +
  guides(
    fill = guide_legend(
      title = "CRS grade", override.aes = list(size = 3)
    )
  )

bcma.crp.pt = sbcma_change %>% drop_na(change, CRP_MAX) %>%
  ggplot(aes(x=log(change,2) , y=CRP_MAX)) +
  geom_point(aes(fill=CRS_GRADE), shape=21, size = 3) +
  geom_smooth(method="lm", formula = 'y ~ x', color="grey20", linewidth = .3) +
  stat_cor(method = "spearman", size=5.5) +
  guides(fill = guide_legend(title = "CRS grade")) +
  ylab("max CRP") + xlab("logFC sBCMA (day30 vs day0)") +
  mytheme(base_size = base.size) +
  theme(
    # aspect.ratio = .618,
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c("#6699CC", "#EE99AA", "#994455")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# B-cell reconstitution
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
missing_timepoints = bcell_df_form %>%
  dplyr::filter(cumulative_days > -1) %>%
  group_by(patient_id) %>%
  slice_min(cumulative_days) %>%
  dplyr::filter(cumulative_days>50) %>%
  pull(patient_id)

bcell_df_form = bcell_df_form %>%
  dplyr::filter(patient_id %in% se.meta@meta.data$PATIENT_ID, !patient_id %in% missing_timepoints)

# bcell_df_form%>%
#   pull(patient_id) %>%
#   unique() %>%
#   length()
# # Total 55

# Time point of B cell aplasia
bcell_min_tp = bcell_df_form %>%
  group_by(patient_id) %>%
  dplyr::filter(patient_id %in% se.meta@meta.data$PATIENT_ID, cumulative_days > -1) %>%
  slice_min(order_by = Wert) %>%
  slice_min(cumulative_days) %>%
  mutate(DATE_MIN_B = Datum) %>%
  dplyr::select(patient_id,cumulative_days, Wert, DATE_MIN_B)

# Only look at time points after B cell aplasia
bcell_df_form = bcell_df_form %>%
  left_join(bcell_min_tp[,c("patient_id", "DATE_MIN_B")]) %>%
  dplyr::filter(Datum >= DATE_MIN_B) %>%
  mutate(
    status = factor(ifelse(Wert>0, 1, 0)),
    PRODUCT = factor(PRODUCT, levels = c("ide", "cilta"))
  )

# Calculate first tp of B cell recovery after possible aplasia
events =  bcell_df_form %>%
  dplyr::filter(status==1) %>%
  group_by(patient_id) %>%
  slice_min(order_by = cumulative_days)

#nrow(events)
#View(events %>% dplyr::select(patient_id, status, cumulative_days, Wert))

censored = bcell_df_form %>%
  dplyr::filter(!patient_id %in% events$patient_id) %>%
  group_by(patient_id) %>%
  slice_max(cumulative_days)
#View(censored %>% dplyr::select(patient_id, status, cumulative_days, Wert))

bcell_rec = rbind(events, censored)
#nrow(bcell_rec)

bcell_rec_pl = cuminc(Surv(cumulative_days, status) ~ PRODUCT, bcell_rec) %>%
  ggcuminc() +
  add_confidence_interval() +
  # add_risktable() +
  add_censor_mark(shape=124, size=3) +
  scale_ggsurvfit(x_scales = list(limits = c(0,320), breaks = seq(0, 300, by = 100))) +
  mytheme(15) +
  theme(legend.position = "inside", legend.position.inside = c(0.1,0.85), plot.title = element_text(face="bold", hjust=.5, size=rel(1.25))) +
  xlab("Days after CAR-T infusion") +
  ylab("B cell recovery") +
  scale_color_manual(values = c(cilta = "#004488", ide = "#44AA99")) +
  scale_fill_manual(values = c(cilta = "#004488", ide = "#44AA99")) +
  guides(colour = guide_legend(
    #legend.position = "inside", legend.position.inside = c(0.1,0.85),
    #title.position="top",
    #title.hjust = 0,
    # keywidth = 1.5,
    # keyheight = .7,
    override.aes = list(size = .6))
  )
#ggtitle("B cell recovery")

coxph(Surv(cumulative_days, as.numeric(status)) ~ PRODUCT, bcell_rec)
crr_mod = coxph(Surv(cumulative_days, as.numeric(status)) ~ PRODUCT + BEST_RESPONSE_CONSENSUS, data=bcell_rec)

tbl <-
  crr_mod %>%
  gtsummary::tbl_regression(exponentiate = TRUE) %>%
#  gtsummary::add_global_p() %>%
  add_n(location = "level")
#tbl

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CAL-1 cellines
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cart_data = read_xlsx(
  "data/Auswerungen Cal-1 killing CAR-T and T cell Teclistamab.xlsx",
  sheet = 2
)

bispec_data = read_xlsx(
  "data/Auswerungen Cal-1 killing CAR-T and T cell Teclistamab.xlsx",
  sheet = 3
)

meta_cart = data.frame(
  name = c("Cilta-cel_2", "Cilta-cel_5", "Cilta-cel_8", "Ide-cel_66"),
  label = c( "P95 (2% Cilta)", "P49 (5% Cilta)", "P49 (8% Cilta)", "P16 (66% Ide)")
)

cart_summary = cart_data %>%
  tidyr::pivot_longer(cols = c(2:6)) %>%
  dplyr::group_by(replicate) %>%
  dplyr::mutate(killing_rate = 1- (value / value[name=="control"]), `% CAR-T` = gsub(".+_", "", name)) %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(
    n=n(),
    mean_killing = mean(killing_rate),
    sd = sd(killing_rate) ,
    se = sd / sqrt(n),
    ci= se * qt(.95/2 + .5, n - 1)
  ) %>%
  dplyr::filter(name != "control") %>%
  dplyr::left_join(meta_cart, by = join_by(name))

cal.kill.pl = cart_summary %>%
  mutate(
    r = seq(1,4,1),
    label = factor(
      cart_summary$label,
      levels = c("P95 (2% Cilta)", "P49 (5% Cilta)","P49 (8% Cilta)", "P16 (66% Ide)")
  )) %>%
  ggplot(aes(x=label, y=100*mean_killing, fill=r)) +
  geom_bar(stat = "identity", color="black") +
  geom_errorbar(aes(ymin=100*(mean_killing-se), ymax = 100*(mean_killing+se)), width=.25) +
  mytheme_grid(base_size = base.size) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1.05),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  ylab("% CAL-1 killing (24h)") +
  scale_y_continuous(limits = c(-21,100), breaks = seq(0,100,20)) +
  scico::scale_fill_scico(palette = "vik", begin=.2, end=.9) +
  guides(fill = "none")

meta_bispec = data.frame(
  group = c("control", "Tec Only", "HD-1", "HD-1", "HD-2", "HD-2", "HD-3", "HD-3"),
  label = c("control", "Tec", "HD-1 T", "HD-1 T+Tec", "HD-2 T", "HD-2 T+Tec", "HD-3 T", "HD-3 T+Tec"),
  name = c("control", "tec", "donor-1", "donor-1_tec", "donor-2", "donor-2_tec", "donor-3", "donor-3_tec"),
  condition = c("CAL-1 only", "Tec only", "T-cells", "T-cells + Tec", "T-cells", "T-cells + Tec", "T-cells", "T-cells + Tec")
)


df.b = bispec_data %>%
  tidyr::pivot_longer(cols = c(2:9)) %>%
  dplyr::left_join(meta_bispec, by="name") %>%
  dplyr::group_by(replicate) %>%
  dplyr::mutate(killing_rate = 1- (value / value[name=="control"])) %>%
  dplyr::filter(name !="control") %>%
  dplyr::mutate(killing_rate = killing_rate * 100)

# tmp = df.b[df.b$group == "HD-1", ]
# tmp = tmp[order(tmp$condition), ]
# wilcox.test(tmp$killing_rate ~ tmp$condition, paired = T)

stat.test <- df.b[df.b$group != "Tec Only", ] %>%
  dplyr::group_by(group) %>%
  rstatix::wilcox_test(killing_rate ~ condition, alternative = "less") %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "condition")

facet.title <- c("HD-1", "HD-2", "HD-3", "Tec")
names(facet.title) <- c("HD-1", "HD-2", "HD-3", "Tec Only")

cal.bispec.pl =
  ggplot(df.b, aes(x=condition, y=killing_rate)) +
  geom_boxplot(aes(fill = condition), outlier.colour = NA, position=position_dodge(0.8), fatten = 1, alpha = .7, lwd = 0.3) +
  geom_dotplot(aes(fill = condition), binaxis='y', stackdir='center', position=position_dodge(0.8), stroke = .5, dotsize = .6) +
  facet_grid(
    ~ group, scale="free_x", space="free_x",
    labeller = labeller(group = facet.title)
  ) +
  # stat_pvalue_manual(stat.test, label = "p", y.position = 105) +
  scale_y_continuous(limits = c(-21,100), breaks = seq(0,100,25)) +
  scale_fill_manual(
    values = c("#004488", "#DDAA33", "#BB5566")
  ) +
  mytheme(base_size = base.size) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank() ,
    legend.position = "bottom",
    panel.spacing = unit(.5, "lines"),
    legend.box.spacing = unit(0, "pt"),
  ) +
  ylab("% CAL-1 killing (24h)") + labs(fill = NULL)

tec.img = ggdraw() + cowplot::draw_image("publication/figures_main/fig_06_tec.png", scale = 1)
car.img = ggdraw() + cowplot::draw_image("publication/figures_main/fig_06_car.png", scale = 1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create the plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="publication/figures_main/fig_06.png",
  plot_grid(
    plot_grid(
      plot_grid(b_cell.pl,NULL,clonesz.pl,nrow=3,rel_heights = c(1.1, .075, 1.1)),
      NULL, plot_grid(
        plot_grid(NULL, shannon.pl, nrow = 2, rel_heights = c(.04, 1)),
        NULL, bar.pl, nrow=3, rel_heights = c(1.0, .075, 1.3),
        labels = c("","C",""), label_fontface = "bold", label_size = 24,hjust = -0.2
      ), NULL,
      plot_grid(NULL, ftrs.pl, nrow = 2, rel_heights = c(.025, 1)),
      ncol = 5, rel_widths = c(0.8,0.1,1,0.1,1.2),
      labels = c("A", "", "B", "", "D"), label_fontface = "bold", label_size = 24, hjust = -0.2
    ),
    plot_grid(NULL),
    plot_grid(bcell.dyn.pl,NULL,pdc.dyn.pl,NULL,plasma.dyn.pl,NULL,ncol=6,rel_widths = c(1,0.1,1,0.1,1,0.01),
              labels = c("E", "", "F", "", "G"), label_fontface = "bold", label_size = 24, hjust = -0.2
    ),
    plot_grid(NULL),
    plot_grid(
      NULL,
      NULL,
      plot_grid(car.img, NULL, nrow = 2, rel_heights = c(.8, .2)),
      NULL,
      cal.kill.pl + theme(plot.margin = margin(t=34, b=28,l=5,r=5, unit = "pt")),
      NULL,
      plot_grid(tec.img, NULL, nrow = 2, rel_heights = c(.8, .2)),
      NULL,
      cal.bispec.pl + theme(plot.margin = margin(t=5, b=80,l=5,r=5, unit = "pt")),
      NULL,
      ncol = 10, rel_widths = c(0.05, -0.04, .4, 0.05, .3, .15, .4, 0.05, .55, .01),
      # ncol = 5, rel_widths = c(0.05, .45, .15, 1, 1.9),
      labels = c("H", "", "", "", "", "", "I", "", "", ""), label_fontface = "bold", label_size = 24, hjust = -0.2
    ),
    NULL,
    plot_grid(bcell.clin.dyn.pl,NULL,bcell_rec_pl,NULL,bcma.dyn.pl,NULL,ncol=6,rel_widths = c(1,0.1,1,0.1,1,0.01),
              labels = c("J", "", "K", "", "L"), label_fontface = "bold", label_size = 24, hjust = -0.2, vjust = 0
    ),
    plot_grid(NULL),
    plot_grid(
      plot_grid(cor.lp.pt, cor.late.pt, cor.verylate.pt, ncol=3),
      NULL,
      bcma.exp.pt,
      NULL,
      bcma.crp.pt,
      NULL,
      ncol = 6, rel_widths = c(1, 0.1, 0.58, 0.025, 0.42, 0.01),
      labels = c("M", "", "N", "", "O", ""), label_fontface = "bold", label_size = 24, hjust = -0.2,vjust = 0
    ),
    nrow = 9, rel_heights = c(1.4, .1, 0.6, .1, .9, 0,  0.8, .1, 0.8)
  ),
  width = 165, height = 180, dpi = 300, bg = "white", units = "mm", scale = 3,
  device = png, type = "cairo"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create the plot as pdf
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="publication/figures_main/fig_06.pdf",
  plot_grid(
    plot_grid(
      plot_grid(b_cell.pl,NULL,clonesz.pl,nrow=3,rel_heights = c(1.1, .075, 1.1)),
      NULL, plot_grid(
        plot_grid(NULL, shannon.pl, nrow = 2, rel_heights = c(.04, 1)),
        NULL, bar.pl, nrow=3, rel_heights = c(1.0, .075, 1.3),
        labels = c("","C",""), label_fontface = "bold", label_size = 18,hjust = -0.2
      ), NULL,
      plot_grid(NULL, ftrs.pl, nrow = 2, rel_heights = c(.025, 1)),
      ncol = 5, rel_widths = c(0.8,0.1,1,0.1,1.2),
      labels = c("A", "", "B", "", "D"), label_fontface = "bold", label_size = 18, hjust = -0.2
    ),
    plot_grid(NULL),
    plot_grid(bcell.dyn.pl,NULL,pdc.dyn.pl,NULL,plasma.dyn.pl,NULL,ncol=6,rel_widths = c(1,0.1,1,0.1,1,0.01),
              labels = c("E", "", "F", "", "G"), label_fontface = "bold", label_size = 18, hjust = -0.2
    ),
    plot_grid(NULL),
    plot_grid(
      NULL,
      NULL,
      plot_grid(car.img, NULL, nrow = 2, rel_heights = c(.8, .2)),
      NULL,
      cal.kill.pl + theme(plot.margin = margin(t=34, b=28,l=5,r=5, unit = "pt")),
      NULL,
      plot_grid(tec.img, NULL, nrow = 2, rel_heights = c(.8, .2)),
      NULL,
      cal.bispec.pl + theme(plot.margin = margin(t=5, b=80,l=5,r=5, unit = "pt")),
      NULL,
      ncol = 10, rel_widths = c(0.05, -0.04, .4, 0.05, .3, .15, .4, 0.05, .55, .01),
      # ncol = 5, rel_widths = c(0.05, .45, .15, 1, 1.9),
      labels = c("H", "", "", "", "", "", "I", "", "", ""), label_fontface = "bold", label_size = 18, hjust = -0.2
    ),
    NULL,
    plot_grid(bcell.clin.dyn.pl,NULL,bcell_rec_pl,NULL,bcma.dyn.pl,NULL,ncol=6,rel_widths = c(1,0.1,1,0.1,1,0.01),
              labels = c("J", "", "K", "", "L"), label_fontface = "bold", label_size = 18, hjust = -0.2, vjust = 0
    ),
    plot_grid(NULL),
    plot_grid(
      plot_grid(cor.lp.pt, cor.late.pt, cor.verylate.pt, ncol=3),
      NULL,
      bcma.exp.pt,
      NULL,
      bcma.crp.pt,
      NULL,
      ncol = 6, rel_widths = c(1, 0.1, 0.58, 0.025, 0.42, 0.01),
      labels = c("M", "", "N", "", "O", ""), label_fontface = "bold", label_size = 18, hjust = -0.2,vjust = 0
    ),
    nrow = 9, rel_heights = c(1.4, .1, 0.6, .1, .9, 0,  0.8, .1, 0.8)
  ),
  width = 165, height = 180, dpi = 300, bg = "white", units = "mm", scale = 3,
  device = cairo_pdf
)
