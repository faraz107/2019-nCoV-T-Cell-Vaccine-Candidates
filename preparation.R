
library(fs)
library(here)
library(Biostrings)
library(purrr)
library(seqinr)

# Reference data -------

ORGANISM_NAME <- "Wuhan seafood market pneumonia virus"
biomartr::is.genome.available(db = "refseq", organism = ORGANISM_NAME)
biomartr::getGFF(db = "refseq", organism = ORGANISM_NAME, gunzip = TRUE, reference = FALSE)

suppressMessages({
  temp <- read_lines(file = here("_ncbi_downloads", "annotation", 
                                 paste0("Wuhan_seafood_market_pneumonia_virus_genomic_refseq.gff")))
  NCBI_nCov <- temp[which(str_detect(temp, "sequence-region"))] %>% word(start = 2, sep = " ")
})

gb_nCov <- genbankr::readGenBank(file = here("Data", "NCBI", "nCov.gb"))

gb_cds <- cds(gb_nCov)

translation_info <- data.frame(gene = gb_cds$gene, 
                               start = gb_cds %>% IRanges::start(),
                               end = gb_cds %>% IRanges::end(), 
                               stringsAsFactors = FALSE)

translation_indices <- translation_info[4:13,]

translation_indices %<>% bind_rows(data.frame(gene = "orf1a", start = 266, end = 13468, stringsAsFactors = FALSE))
translation_indices %<>% bind_rows(data.frame(gene = "orf1b", start = 13468, end = 21555, stringsAsFactors = FALSE))

translation_indices %<>% arrange(start)

# Analysis --------

dir_create(path = here("Data", "Filtered"))

raw_ls <- dir_ls(path = here("Data", "Raw"))

exclude_ls <- c("EPI_ISL_402121", "EPI_ISL_402126", "EPI_ISL_403928", "EPI_ISL_403931",
                "EPI_ISL_402120", "EPI_ISL_406959", "EPI_ISL_406960")

filtered_ls <- vector(mode = "character")

for (f in raw_ls) {
  if (sum(str_detect(string = f, pattern = exclude_ls)) == 0) {
    filtered_ls <- c(filtered_ls, f)
    temp <- str_split(string = f, pattern = "/")[[1]]
    file_copy(path = f, new_path = here("Data", "Filtered", temp[length(temp)]), overwrite = TRUE)
  }
  
}

genomes <- list()
g_names <- list()

for (f in filtered_ls) {
    g <- read.fasta(file = f, seqtype = "DNA", whole.header = TRUE, as.string = TRUE)
    genomes <- c(genomes, as.character(g))
    g_names <- c(g_names, names(g))
  }

seqinr::write.fasta(sequences = genomes, 
                    names = g_names, 
                    as.string = FALSE, 
                    file.out = here("mafft-linux", "input"))

# setwd(here("mafft-linux"))
# system("mafft input > output")
# 
# dir_create(path = here("Data", "MSAs"))
# 
# file_copy(path = here("mafft-linux", "output"), 
#           new_path = here("Data", "MSAs", "genome.fasta"), overwrite = TRUE)
# 
# setwd(here())

# Using UGENE to remove ensure alignment is along refernce Wuhan-Hu-1 sequence

data_file <- here("Data", "MSAs", "genome_referenced.fasta")
align_MSA <- read.alignment(file = data_file, format = "fasta")

DNA_Mat <- toupper(as.matrix.alignment(x = align_MSA))

for (l in seq.int(1,length(translation_indices$gene))) {
prepare.proteinMSA(DNA_Mat = DNA_Mat, 
                   PROTEIN = translation_indices$gene[l], 
                   START = translation_indices$start[l], 
                   END = translation_indices$end[l])
}



  
