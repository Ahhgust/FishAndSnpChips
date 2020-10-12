#!/usr/local/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailing=TRUE)

if (length(args)<2) {
    stop("I need 1 .csv file and one index (column) to work!")
}

file <- args[[1]]
index <- as.integer(args[[2]]) - 1 # 0-based index (internally), 1-based (externally)

if (is.na(index)) {
    stop("Problem parsing index")
}

annos <- read_tsv("keystone.annos.GC.10.20.100.encodeHighsignal.chromHMM.bed.hg38.gz", na=".", col_names=FALSE, progress=FALSE, col_types=cols())
colnames(annos) <- c("Chr", "Pos0", "Position", "foo", "GC10", "GC20", "GC100", "EncodeHigh", "ChromState")

# the annotations refer to the SNPs in the snp-chip. Apparently some SNPs are duplicated
# which lead to dup records in the annotations. This should fix the dups in the annos
annos <- distinct(annos)


header <- read_csv(file,
                   progress=FALSE, col_types=cols(), n_max=1)

# file format:
#cols 1:10 <- info on site;


ncols <- length(colnames(header))

# ok. col_types can take a string (1-char per column)
# ? means keep it, guess the type, _ means skip it
# this constructs said string
tibble(Cols=1:ncols,
       Incl=Cols < 6 | between(Cols, (index*4) + 11, (index*4) + 14),
       Char=ifelse(Incl==TRUE, "?", "_")) %>%
           pull(Char) %>%
           paste0(collapse="") -> coltypes
    
# each person is assocaited w/ 4 columns: in the order, Genotype, Score, Theta, R (*4 in above formula)


# pull out the position info (cols 1-5, and the person's info 
tib <- read_csv(file,
                col_types=coltypes,
                guess_max=2612359,
                progress=FALSE)


collies <- colnames(tib)

nc <- ncol(tib)

# gives the sample name (who) + ".R"
sname <- collies[nc]
# remove .R
sname <- substr(sname, 0, nchar(sname)-2)

colnames(tib)[ nc ] <- "R"
colnames(tib)[ nc -1 ] <- "Theta"
colnames(tib)[ nc -2 ] <- "Score"
colnames(tib)[ nc -3 ] <- "Genotype"
                
# remove non-standard chromosomes; stick to just the autosomes...
filter(tib,
       Chr %in% 1:22) -> tib

# how many no-calls are there? (consider a min of 1 as we're taking logs)
nocallCount <- max(sum( tib$Genotype=="NC"), 1)

# go with UCSC style chromosome naming...
tib$Chr <- paste0("chr", tib$Chr)

combined <- inner_join(
    tib,
    annos,
    by=c("Chr", "Position"))

combined$LogNocall <- log10(nocallCount)

mutate(combined,
       EncodeHigh=ifelse( is.na(EncodeHigh), "NA", EncodeHigh),
       ChromState=ifelse( is.na(ChromState), "NA", ChromState),
       Genotype=factor(Genotype, levels=c("AA", "AB", "BB", "NC")),
       S=tan( (Theta*pi)/2),
       X=R/(S+1),
       Y=(R*S)/(S+1)) -> combined

combined <- drop_na(combined)


mod <- readRDS("ElasticNet.withnocalls.genotyping.rds")

predict(mod, newdata=combined, type='prob', s= mod$lambdaOpt) -> pres   

likeToGL <- function(x,y,z) {
    sprintf("%.6f,%0.6f,%0.6f", log10(1-x), log10(1-y), log10(1-z))
}

# takes three genotype likelihoods; computes quality of site.
likeToQual <- function(x) {
    ifelse(x>0.9999,
           99,
           ifelse(x < 0.0001, 1,
                  round(-10*log10(1- x ))
                  )
           )
}

combined %>% 
    select(Chr, Position, Name, Genotype) %>%
    mutate(Chr=factor(Chr, levels=paste0("chr", 1:22))) %>%
    bind_cols(pres) %>%
    arrange(Chr, Position) %>%
    mutate(
        MaxLike=pmax(
            AA,
            AB,
            BB),
        MaxGeno=case_when(
            AA==MaxLike ~ "0/0",
            AB==MaxLike ~ "0/1",
            BB==MaxLike ~ "1/1",
            TRUE        ~ "?"
        ))    -> all # TODO: actually two cases; AA can mean ref/ref or alt/alt. Need more info!!


all %>%
    transmute(
       `#CHROM`=Chr,
       POS=Position,
       ID=Name,
       REF="A",
       ALT="B",
       QUAL=likeToQual(MaxLike),
       FILTER=ifelse(QUAL>19, "PASS", "q20"),
       INFO=paste("OriginalGeno", Genotype,sep="="),
       FORMAT="GT:GL",
       !!sname:=paste(MaxGeno,
                      likeToGL(AA, AB, BB),
                      sep=":")
    ) -> asvcf

filename <- paste(sname, "vcf", sep=".")

# todo; flesh out header...
cat("##fileformat=VCFv4.3\n##fileDate=2020\n##source=illumina2vcf\n", file=filename)

write_tsv(asvcf, path=filename, append=TRUE)




