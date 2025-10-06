library(dplyr)
library(yaml)

options(scipen=999)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq: CellRanger count
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nodes = 1
ntasks = 1
ttime = "66:00:00"
mail = "FAIL"
mem = 200000
cpu = 30

manifest = yaml.load_file("manifest.yaml")

cellranger = manifest$ukl_b4$data_dl
ref.fasta = paste0(manifest$refs, "reference_sources/Homo_sapiens.GRCh38.dna.primary_assembly.fa")

out.dir = paste0(manifest$ukl_b4$work, "souporcell/output/")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sbatch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_subscript = function(path, job_id, sample.paths){
  file.create(path)

  write("#!/bin/bash", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(paste("#SBATCH -J", job_id), file = path, append = TRUE)
  write(paste("#SBATCH --nodes", nodes), file = path, append = TRUE)
  write(paste("#SBATCH --ntasks", ntasks), file = path, append = TRUE)
  write(paste("#SBATCH --time", ttime), file = path, append = TRUE)
  write(paste("#SBATCH --cpus-per-task", cpu), file = path, append = TRUE)
  write(paste("#SBATCH --mem", mem), file = path, append = TRUE)
  write(paste("#SBATCH --exclude=ribnode[009,020,006,010,012,017]"), file = path, append = TRUE)
  write(paste("#SBATCH -e", paste0(job_id, ".e")), file = path, append = TRUE)
  write(paste("#SBATCH -o", paste0(job_id, ".o")), file = path, append = TRUE)
  write("#SBATCH --mail-type=END,FAIL", file = path, append = TRUE)

  write("", file = path, append = TRUE)
  write("ml Singularity", file = path, append = TRUE)
  write("export SINGULARITY_BINDPATH=/mnt/ribolution/", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(
    paste0(
      "singularity exec ", manifest$base$workdata, "/cohorts/souporcell_latest.sif souporcell_pipeline.py",
      " -i ", sample.paths$bam.PATH ,
      " -b ", sample.paths$barcode.PATH,
      " -f ", ref.fasta,
      " -t ", 30,
      " -o ", paste0(out.dir, sample.paths$SAMPLE),
      " -k ", 2
  ), file = path, append = TRUE)

  write("", file = path, append = TRUE)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sample paths
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bam.files = list.files(
  path = cellranger, pattern = "sample_alignments.bam$", full.names = T,
  recursive = T
)
barcode.files = list.files(
  path = cellranger, pattern = "barcodes.tsv.gz$",
  full.names = T, recursive = T
)
barcode.files = barcode.files[grepl("count/sample_filtered_feature_bc_matrix", barcode.files)]

df = data.frame(
  SAMPLE = basename(dirname(dirname(bam.files))),
  bam.PATH = bam.files,
  barcode.PATH = barcode.files
)

stopifnot(identical(
  basename(dirname(dirname(dirname(df$barcode.PATH)))),
  basename(dirname(dirname(df$bam.PATH)))
))
stopifnot(identical(
  df$SAMPLE,
  basename(dirname(dirname(df$bam.PATH)))
))

# k = c("multi_T3-A", "multi_T3-B", "multi_T3-D")
# df = df[df$SAMPLE %in% k, ]
# table(df$SAMPLE)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# submit to Slurm
# There is a link with this script under this path: "out.dir".
# The script was executed there.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (sample in unique(df$SAMPLE)) {
  job_id = paste0("souporcell_", sample)
  print(job_id)
  sample.paths = subset(df, SAMPLE == sample)
  write_subscript(paste0(out.dir, job_id, ".slurm"), job_id, sample.paths)

  cmd = paste0("sbatch ", paste0(out.dir, job_id, ".slurm"))
  system(cmd)
}
