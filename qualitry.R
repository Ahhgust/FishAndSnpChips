#!/usr/local/bin/Rscript


suppressPackageStartupMessages(library(tidyverse))

theme_set(theme_bw(base_size=16))







library(scales)
library(glmnet)

# Piggy piggy. read in the whole file.
tib <- read_csv("keystone192.2.csv", guess_max=2612359)

test <- head(tib)

tib %>% select(Name, Chr, Position, contains("GType")) %>%
    gather(key="Individual", value="Genotype",  -Name, -Chr, -Position) -> genotypes

tib %>% select(Name, Chr, Position, contains(".Score")) %>%
    gather(key="Individual", value="Score",  -Name, -Chr, -Position) -> scores

tib %>% select(Name, Chr, Position, contains(".Theta")) %>%
    gather(key="Individual", value="Score",  -Name, -Chr, -Position) -> thetas

tib %>% select(Name, Chr, Position, contains(".R")) %>%
    gather(key="Individual", value="Score",  -Name, -Chr, -Position) -> rs

rm(tib)
gc()

# this is clobbering the ID
genotypes %>% separate(Individual, into=c("ID"), sep="[.]", extra="drop") -> genotypes

scores %>% separate(Individual, into=c("ID"), sep="[.]", extra="drop") -> scores

thetas %>% separate(Individual, into=c("ID"), sep="[.]", extra="drop") -> thetas
rs %>% separate(Individual, into=c("ID"), sep="[.]", extra="drop") -> rs


tib.long <- bind_cols(genotypes,
                      select(scores, Score),
                      select(thetas, Score) %>% rename(Theta=Score),
                      select(rs, Score) %>% rename(R=Score))
                      



tib.long %>%
    filter( grepl("_[0-9]+ng", ID)) %>%
    separate(ID, into=c("SampPre", "SampPost", "ng"), sep="_") ->
    dilutionSeries

dilutionSeries %>%
    write_tsv("long.parsed.tsv.gz")

dilutionSeries <- read_tsv("long.parsed.tsv.gz", guess_max=62696569)


isCorrect <- function(ng, genotype, idx) {
    truth <- genotype[ which(ng=='200ng') ]
    obs <- genotype[ idx ]
    return(obs==truth)
}

filter(dilutionSeries,
       ng=='200ng')  %>%
    select(-ng) %>%
    rename(CorrectGenotype=Genotype,
           CorrectScore=Score) -> correct


left_join(
    filter(dilutionSeries, ng != '200ng'),
    correct,
    by=c("Name", "Chr", "Position", "SampPre", "SampPost")) ->
        dSeries

select(dSeries, -Theta.y, -R.y) %>%
    rename(Theta=Theta.x, R=R.x) %>%
    filter(Genotype != "NC",
           CorrectGenotype != "NC",
        Chr %in% 1:22) -> dSeries.autosomes

nrow(dSeries.autosomes)

filter(dSeries, Chr %in% 1:22) %>%
    group_by(SampPre, SampPost, ng) %>%
    summarize(NumNocall=sum(Genotype=="NC")) ->nocallCount

#ggplot(nocallCount, aes(x=ng, y=NumNocall)) + geom_jitter() 

annos <- read_tsv("keystone.annos.GC.10.20.100.encodeHighsignal.chromHMM.bed", na=".", col_names=FALSE, progress=FALSE, col_types=cols())

colnames(annos) <- c("Chr", "Pos0", "Position", "foo", "GC10", "GC20", "GC100", "EncodeHigh", "ChromState")

#lift <- read_tsv("keystoneAsBed.hg19.bed", col_names=FALSE, progress=FALSE, col_types=cols())

#colnames(lift) <- c("Chr", "Pos0", "PosHg38", "foo")

#separate(lift, foo, into=c("ChrHg19", "Pos0too", "PosHg19")) -> lift

#filter(lift,
 #      Chr==ChrHg19) %>%
  #  select(Chr, PosHg19, PosHg38) -> lift

#mutate(lift,
#       PosHg19=as.numeric(PosHg19) ) -> lift


#inner_join(annos,
 #          lift, by=c("Chr"="Chr", "Position"="PosHg38")) -> annos


nocallCount$ngF <- factor(nocallCount$ng, levels=c("1ng", "5ng", "10ng", "25ng", "50ng", "100ng", "150ng"))

ggplot(nocallCount, aes(x=ngF, y=log(NumNocall))) + geom_jitter()


select(annos,
       -Pos0,
       -foo
       ) %>%
    mutate(Chr=substr(Chr, 4, 999)) %>% # remove chr prefix
    inner_join(
        select( dSeries.autosomes, -CorrectScore),
        by=c("Chr", "Position")
        )  -> autosWithAnnos

inner_join(
    autosWithAnnos,
    nocallCount,
    by=c("SampPre", "SampPost", "ng")) -> autosWithAnnosAndNocalls


mutate(autosWithAnnosAndNocalls,
       EncodeHigh=ifelse( is.na(EncodeHigh), "NA", EncodeHigh),
       ChromState=ifelse( is.na(ChromState), "NA", ChromState),
       Genotype=factor(Genotype, levels=c("AA", "AB", "BB")),
       CorrectGenotype=factor(CorrectGenotype, levels=c("AA", "AB", "BB")),
       S=tan( (Theta*pi)/2),
       X=R/(S+1),
       Y=(R*S)/(S+1),
       LogNocall=log10(NumNocall)) -> autosTransformed
       
       
# TODO: polar coordinates (Theta, R)
# vs cartesian coordinates (X,Y)

xtrain <- model.matrix(
    ~ . -1, select(
                autosTransformed,
                GC10,GC20, GC100, Score, Theta, R, LogNocall,
                EncodeHigh, ChromState, Genotype, Score)
)


mod <- cv.glmnet(x=xtrain,y = autosTransformed$CorrectGenotype, family =   "multinomial", alpha=0.5, nfolds=3)

preds <- predict(mod, newx=xtrain, s="lambda.min", type="response")
predClass <- predict(mod, newx=xtrain, s="lambda.min", type="class") 

sum(predClass==autosWithAnnosAndNocalls$CorrectGenotype)/nrow(autosWithAnnosAndNocalls)

sum(autosWithAnnosAndNocalls$Genotype==autosWithAnnosAndNocalls$CorrectGenotype)/nrow(autosWithAnnosAndNocalls)


mod <- cv.glmnet(x=xtrain,y = autosTransformed$CorrectGenotype, family =   "multinomial", alpha=1, nfolds=3)

preds <- predict(mod, newx=xtrain, s="lambda.min", type="response")
predClass <- predict(mod, newx=xtrain, s="lambda.min", type="class") 

sum(predClass==autosWithAnnosAndNocalls$CorrectGenotype)/nrow(autosWithAnnosAndNocalls)

sum(autosWithAnnosAndNocalls$Genotype==autosWithAnnosAndNocalls$CorrectGenotype)/nrow(autosWithAnnosAndNocalls)



