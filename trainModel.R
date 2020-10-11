#!/usr/local/bin/Rscript


suppressPackageStartupMessages(library(tidyverse))

theme_set(theme_bw(base_size=16))


library(scales)
library(glmnet)
library(caret)

# let us do 5-fold cross validation in parallel
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)


dilutionSeries <- read_tsv("long.parsed.tsv.gz", guess_max=5367227)


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
    filter(#Genotype != "NC",
           CorrectGenotype != "NC",
        Chr %in% 1:22) -> dSeries.autosomes

nrow(dSeries.autosomes)

filter(dSeries, Chr %in% 1:22) %>%
    group_by(SampPre, SampPost, ng) %>%
    summarize(NumNocall=sum(Genotype=="NC")) ->nocallCount

#ggplot(nocallCount, aes(x=ng, y=NumNocall)) + geom_jitter() 

annos <- read_tsv("keystone.annos.GC.10.20.100.encodeHighsignal.chromHMM.bed", na=".", col_names=FALSE, progress=FALSE, col_types=cols())

colnames(annos) <- c("Chr", "Pos0", "Position", "foo", "GC10", "GC20", "GC100", "EncodeHigh", "ChromState")


nocallCount$ngF <- factor(nocallCount$ng, levels=c("1ng", "5ng", "10ng", "25ng", "50ng", "100ng", "150ng"))



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
       Genotype=factor(Genotype, levels=c("AA", "AB", "BB", "NC")),
       CorrectGenotype=factor(CorrectGenotype, levels=c("AA", "AB", "BB")),
       S=tan( (Theta*pi)/2),
       X=R/(S+1),
       Y=(R*S)/(S+1),
       LogNocall=log10(NumNocall)) -> autosTransformed
       


set.seed(123)


autosTransformed$Position <-  as.integer(autosTransformed$Position)

training <- filter(autosTransformed,
                   row_number() %% 10 <7)

test <- filter(autosTransformed,
                   row_number() %% 10 > 6)

training <- training[ complete.cases(training),]
test <- test[ complete.cases(test),]
    
trct <- trainControl(method="cv", number=5, returnResamp="all",
                     classProbs=TRUE, summaryFunction=mnLogLoss)

grid <- expand.grid(alpha=0:10/10,
                    lambda = 1/2**(18:10))

# train the elastic net over the search grid
# consider every variable and their cross-products
cvMod <- train(
    CorrectGenotype ~ (Genotype + GC10 + GC20 + GC100 + EncodeHigh +ChromState + Genotype + Score + Theta + R + LogNocall + X + Y)^2,
    data=data.frame(training),
    method = "glmnet", 
    trControl = trct,
    metric = "logLoss",
    tuneGrid = grid)


bestParams <- cvMod$bestTune

predict(cvMod, newdata=test, type='prob', s=  cvMod$lambdaOpt) -> pres

mutate(pres,
       MaxLike=pmax(
           AA,
           AB,
           BB),
       MaxGeno=case_when(
           AA==MaxLike ~ "AA",
           AB==MaxLike ~ "AB",
           BB==MaxLike ~ "BB",
           TRUE        ~ "?"
           )
       ) %>%
    bind_cols(test) -> testWithAnswers

# naive calling
mutate(testWithAnswers,
       Correct=
           as.character(CorrectGenotype)==
       as.character(Genotype)) %>%
    group_by(ng) %>%
    count(Correct)


mutate(testWithAnswers,
       Correct=as.character(CorrectGenotype)==
           as.character(MaxGeno)) %>%
    group_by(ng) %>%
    count(Correct)


saveRDS(cvMod, file="ElasticNet.withnocalls.genotyping.rds")

# x is the probability that the basecall is correct
# converts to a phred-based call.
likeToPL <- function(x) {
    pmin(
        round( -10*log10(1-x) ),
        99)
}

likeToGL <- function(x) {
    log10(1-x)
}


