#import raw data and library
library(dplyr)

#setwd to where the downloaded repository is 
plasmid_sequence <-read.delim('data.human TET2 sequence single strand.txt', header = FALSE, stringsAsFactors = FALSE)

#assign value to b: b is the position in a gene where we want primers;
#here is an example: we want primers around 1680bp, so b=1680, and you can use any numbers

b <- 1680

#----------------------------------------------------------------------------------------------------------------
#clean up raw data and convert to data frame format (from 5' to 3')

plasmid_sequence <- strsplit(plasmid_sequence[1,1],"")
plasmid_sequence <- data.frame(plasmid_sequence[1], row.names = NULL, 
                               check.rows = FALSE, check.names = TRUE, 
                               fix.empty.names = TRUE, stringsAsFactors=FALSE)

#based on the plasmid sequence, generate all primers around b position, and named "primer_precursors". 
#a = primer length -1

a <- 17
primer_precursor_17 <- NULL
for (i in 0:a){
  primer_precursor_17 <- rbind(primer_precursor_17, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_17 <- data.frame(primer_precursor_17, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 18
primer_precursor_18 <- NULL
for (i in 0:a){
  primer_precursor_18 <- rbind(primer_precursor_18, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_18 <- data.frame(primer_precursor_18, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 19
primer_precursor_19 <- NULL
for (i in 0:a){
  primer_precursor_19 <- rbind(primer_precursor_19, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_19 <- data.frame(primer_precursor_19, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 20
primer_precursor_20 <- NULL
for (i in 0:a){
  primer_precursor_20 <- rbind(primer_precursor_20, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_20 <- data.frame(primer_precursor_20, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 21
primer_precursor_21 <- NULL
for (i in 0:a){
  primer_precursor_21 <- rbind(primer_precursor_21, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_21 <- data.frame(primer_precursor_21, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 22
primer_precursor_22 <- NULL
for (i in 0:a){
  primer_precursor_22 <- rbind(primer_precursor_22, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_22 <- data.frame(primer_precursor_22, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 23
primer_precursor_23 <- NULL
for (i in 0:a){
  primer_precursor_23 <- rbind(primer_precursor_23, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

primer_precursor_23 <- data.frame(primer_precursor_23, row.names = NULL, 
                                  check.rows = FALSE, check.names = TRUE, 
                                  fix.empty.names = TRUE, stringsAsFactors=FALSE)

#combine all primers generated with a length range from 18bp to 22 bp
primer_precursors <- bind_rows(primer_precursor_17, primer_precursor_18, primer_precursor_19, 
                               primer_precursor_20, primer_precursor_21, primer_precursor_22, 
                               primer_precursor_23, .id = NULL)


#---------------------------------------------------------------------------------------------------------------
#make complememt plamisd sequence, and extract all primers in the complemented version of this plasmid sequence (from 3' to 5')

plasmid_sequence <-read.delim('data.human TET2 sequence single strand.txt', header = FALSE, stringsAsFactors = FALSE)
plasmid_sequence <- strsplit(plasmid_sequence[1,1],"")
plasmid_sequence <- paste(unlist(plasmid_sequence), collapse='')
complement_plasmid_sequence <- chartr("ATGC","TACG", plasmid_sequence)
complement_plasmid_sequence <- as.character(complement_plasmid_sequence)
complement_plasmid_sequence <- strsplit(complement_plasmid_sequence,"")
complement_plasmid_sequence <- data.frame(complement_plasmid_sequence[1], row.names = NULL, 
                                          check.rows = FALSE, check.names = TRUE, 
                                          fix.empty.names = TRUE, stringsAsFactors=FALSE)

a <- 17
complement_primer_precursor_17 <- NULL
for (i in 0:a){
  complement_primer_precursor_17 <- rbind(complement_primer_precursor_17, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
complement_primer_precursor_17 <- data.frame(complement_primer_precursor_17, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
#reverse the columns of this dataframe, so that primers are from 5' to 3'
complement_primer_precursor_17 = complement_primer_precursor_17[ ,order(ncol(complement_primer_precursor_17):1)]
colnames(complement_primer_precursor_17) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 'X18')

a <- 18
complement_primer_precursor_18 <- NULL
for (i in 0:a){
  complement_primer_precursor_18 <- rbind(complement_primer_precursor_18, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
complement_primer_precursor_18 <- data.frame(complement_primer_precursor_18, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_18 = complement_primer_precursor_18[ ,order(ncol(complement_primer_precursor_18):1)]
colnames(complement_primer_precursor_18) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 'X18', 'X19')

a <- 19
complement_primer_precursor_19 <- NULL
for (i in 0:a){
  complement_primer_precursor_19 <- rbind(complement_primer_precursor_19, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
complement_primer_precursor_19 <- data.frame(complement_primer_precursor_19, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_19 = complement_primer_precursor_19[ ,order(ncol(complement_primer_precursor_19):1)]
colnames(complement_primer_precursor_19) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 'X18', 'X19','X20')

a <- 20
complement_primer_precursor_20 <- NULL
for (i in 0:a){
  complement_primer_precursor_20 <- rbind(complement_primer_precursor_20, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

complement_primer_precursor_20 <- data.frame(complement_primer_precursor_20, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_20 = complement_primer_precursor_20[ ,order(ncol(complement_primer_precursor_20):1)]
colnames(complement_primer_precursor_20) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 'X18', 'X19','X20','X21')

a <- 21
complement_primer_precursor_21 <- NULL
for (i in 0:a){
  complement_primer_precursor_21 <- rbind(complement_primer_precursor_21, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
complement_primer_precursor_21 <- data.frame(complement_primer_precursor_21, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_21 = complement_primer_precursor_21[ ,order(ncol(complement_primer_precursor_21):1)]
colnames(complement_primer_precursor_21) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 
                                              'X18', 'X19','X20','X21','X22')

a <- 22
complement_primer_precursor_22 <- NULL
for (i in 0:a){
  complement_primer_precursor_22 <- rbind(complement_primer_precursor_22, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}

complement_primer_precursor_22 <- data.frame(complement_primer_precursor_22, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_22 = complement_primer_precursor_22[ ,order(ncol(complement_primer_precursor_22):1)]
colnames(complement_primer_precursor_22) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 
                                              'X18', 'X19','X20','X21','X22','X23')

a <- 23
complement_primer_precursor_23 <- NULL
for (i in 0:a){
  complement_primer_precursor_23 <- rbind(complement_primer_precursor_23, 
                                          complement_plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
complement_primer_precursor_23 <- data.frame(complement_primer_precursor_23, row.names = NULL, 
                                             check.rows = FALSE, check.names = TRUE, 
                                             fix.empty.names = TRUE, stringsAsFactors=FALSE)
complement_primer_precursor_23 = complement_primer_precursor_23[ ,order(ncol(complement_primer_precursor_23):1)]
colnames(complement_primer_precursor_23) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 
                                              'X12', 'X13', 'X14', 'X15', 'X16', 'X17', 
                                              'X18', 'X19','X20','X21','X22','X23','X24')

complement_primer_precursors <- bind_rows(complement_primer_precursor_17, complement_primer_precursor_18, 
                                          complement_primer_precursor_19, complement_primer_precursor_20, 
                                          complement_primer_precursor_21, complement_primer_precursor_22, 
                                          complement_primer_precursor_23, .id = NULL)




#####################################################################################################################
##Emma Code for primer precursor, GC
#We will first compute the GC content of the primer precursors generated from 18bp-24bp
#We will use the R package,  seqinr  which is frequently used by biologists for sequence manipulation

install.packages("seqinr")
library(seqinr)

#calculate the gc content of each row in 18bp primers
pp_18_GC <-apply(as.matrix(primer_precursor_17), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 19bp primers
pp_19_GC <-apply(as.matrix(primer_precursor_18), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 20bp primers
pp_20_GC <-apply(as.matrix(primer_precursor_19), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 21bp primers
pp_21_GC <-apply(as.matrix(primer_precursor_20), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 22bp primers
pp_22_GC <-apply(as.matrix(primer_precursor_21), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 23bp primers
pp_23_GC <-apply(as.matrix(primer_precursor_22), MARGIN = 1, FUN = GC)

#calculate the gc content of each row in 24bp primers
pp_24_GC <-apply(as.matrix(primer_precursor_23), MARGIN = 1, FUN = GC)

#merge data so that the potential primers have their associated GC
primer_precursor_17["GC Content"] <- pp_18_GC
primer_precursor_18["GC Content"] <- pp_19_GC
primer_precursor_19["GC Content"] <- pp_20_GC
primer_precursor_20["GC Content"] <- pp_21_GC
primer_precursor_21["GC Content"] <- pp_22_GC
primer_precursor_22["GC Content"] <- pp_23_GC
primer_precursor_23["GC Content"] <- pp_24_GC

#remove primers that are not within 0.4-0.6(40%-60%) range for GC content
#we can do this through dplyr
install.packages('tidyverse')
library(dplyr)

#if the GC content >= 0.4 and <=0.6, it will create a new file called pp_(bplength)bp_pass_GC_threshold
##output are primers that contain 0.4-0.6 (or, 40%-60%) GC content

#18bp
pp_18bp_pass_GC_threshold <- filter(primer_precursor_17, primer_precursor_17$`GC Content` >= 0.4, primer_precursor_17$`GC Content`<= 0.6)

#19bp
pp_19bp_pass_GC_threshold <- filter(primer_precursor_18, primer_precursor_18$`GC Content` >= 0.4, primer_precursor_18$`GC Content`<= 0.6)

#20bp
pp_20bp_pass_GC_threshold <- filter(primer_precursor_19, primer_precursor_19$`GC Content` >= 0.4, primer_precursor_19$`GC Content`<= 0.6)

#21bp
pp_21bp_pass_GC_threshold <- filter(primer_precursor_20, primer_precursor_20$`GC Content` >= 0.4, primer_precursor_20$`GC Content`<= 0.6)

#22bp
pp_22bp_pass_GC_threshold <- filter(primer_precursor_21, primer_precursor_21$`GC Content` >= 0.4, primer_precursor_21$`GC Content`<= 0.6)

#23bp
pp_23bp_pass_GC_threshold <- filter(primer_precursor_22, primer_precursor_22$`GC Content` >= 0.4, primer_precursor_22$`GC Content`<= 0.6)

#24bp
pp_24bp_pass_GC_threshold <- filter(primer_precursor_23, primer_precursor_23$`GC Content` >= 0.4, primer_precursor_23$`GC Content`<= 0.6)

################################################################################################################################################################################################################################
#We will now compute the GC content of the complementary primer precursors generated from 18bp-24bp
#workflow is similar to primer precursors

library(seqinr)
#calculate the gc content of each row in 18bp primers
cpp_18_GC <-apply(as.matrix(complement_primer_precursor_17), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 19bp primers
cpp_19_GC <-apply(as.matrix(complement_primer_precursor_18), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 20bp primers
cpp_20_GC <-apply(as.matrix(complement_primer_precursor_19), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 21bp primers
cpp_21_GC <-apply(as.matrix(complement_primer_precursor_20), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 22bp primers
cpp_22_GC <-apply(as.matrix(complement_primer_precursor_21), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 23bp primers
cpp_23_GC <- apply(as.matrix(complement_primer_precursor_22), MARGIN = 1, FUN  = GC)

#calculate the gc content of each row in 24bp primers
cpp_24_GC <-apply(as.matrix(complement_primer_precursor_23), MARGIN = 1, FUN  = GC)

#merge data so that the potential primers have their associated GC
complement_primer_precursor_17["GC Content"] <- cpp_18_GC
complement_primer_precursor_18["GC Content"] <- cpp_19_GC
complement_primer_precursor_19["GC Content"] <- cpp_20_GC
complement_primer_precursor_20["GC Content"] <- cpp_21_GC
complement_primer_precursor_21["GC Content"] <- cpp_22_GC
complement_primer_precursor_22["GC Content"] <- cpp_23_GC
complement_primer_precursor_23["GC Content"] <- cpp_24_GC

#remove primers that are not within 0.4-0.6(40%-60%) range for GC content
#we can do this through dplyr

library(dplyr)

#if the GC content >= 0.4 and <=0.6, it will create a new file called pp_(bplength)bp_pass_GC_threshold
##output are primers that contain 0.4-0.6 (or, 40%-60%) GC content

#18bp
cpp_18bp_pass_GC_threshold <- filter(complement_primer_precursor_17, complement_primer_precursor_17$`GC Content` >= 0.4, complement_primer_precursor_17$`GC Content`<= 0.6)

#19bp
cpp_19bp_pass_GC_threshold <- filter(complement_primer_precursor_18, complement_primer_precursor_18$`GC Content` >= 0.4, complement_primer_precursor_18$`GC Content`<= 0.6)

#20bp
cpp_20bp_pass_GC_threshold <- filter(complement_primer_precursor_19, complement_primer_precursor_19$`GC Content` >= 0.4, complement_primer_precursor_19$`GC Content`<= 0.6)

#21bp
cpp_21bp_pass_GC_threshold <- filter(complement_primer_precursor_20, complement_primer_precursor_20$`GC Content` >= 0.4, complement_primer_precursor_20$`GC Content`<= 0.6)

#22bp
cpp_22bp_pass_GC_threshold <- filter(complement_primer_precursor_21, complement_primer_precursor_21$`GC Content` >= 0.4, complement_primer_precursor_21$`GC Content`<= 0.6)

#23bp
cpp_23bp_pass_GC_threshold <- filter(complement_primer_precursor_22, complement_primer_precursor_22$`GC Content` >= 0.4, complement_primer_precursor_22$`GC Content`<= 0.6)

#24
cpp_24bp_pass_GC_threshold <- filter(complement_primer_precursor_23, complement_primer_precursor_23$`GC Content` >= 0.4, complement_primer_precursor_23$`GC Content`<= 0.6)


#############################################################################################################################################
#now to calculate the Tm of the primer precursors

###calculate the Tm of a primer through the TmCalculator package
install.packages("TmCalculator")
library(TmCalculator)
#
##there are three different ways to calcualte Tm, based off of GC/based off of other nucleotides, and the Wallace Method (good for primer 14-20bp)

##we will use the Tm_NN function which takes into account the bases as well as ionic concentrations of the PCR reaction
#you can change the input of the pcr reaction ionic concentrations using this method
##we will need to convert into strings in order for this function to work

install.packages("Biostrings")
library(Biostrings)


####create empty vector and then run function over the rows
#In the for loop the colon creates a sequence of numbers from 1 to 15, and loops over i = 1, i= 2, .... 

#18bp
pp_18bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_18bp_pass_GC_threshold)){
  pp_18bp_pass_GC_threshold_Tm <- c(pp_18bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_18bp_pass_GC_threshold[i,]))))
}

#19bp
pp_19bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_19bp_pass_GC_threshold)){
  pp_19bp_pass_GC_threshold_Tm <- c(pp_19bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_19bp_pass_GC_threshold[i,]))))
}

#20bp
pp_20bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_20bp_pass_GC_threshold)){
  pp_20bp_pass_GC_threshold_Tm <- c(pp_20bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_20bp_pass_GC_threshold[i,]))))
}

#21bp
pp_21bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_21bp_pass_GC_threshold)){
  pp_21bp_pass_GC_threshold_Tm <- c(pp_21bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_21bp_pass_GC_threshold[i,]))))
}

#22bp
pp_22bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_22bp_pass_GC_threshold)){
  pp_22bp_pass_GC_threshold_Tm <- c(pp_22bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_22bp_pass_GC_threshold[i,]))))
}

#23bp
pp_23bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_23bp_pass_GC_threshold)){
  pp_23bp_pass_GC_threshold_Tm <- c(pp_23bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_23bp_pass_GC_threshold[i,]))))
}

#24bp
pp_24bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(pp_24bp_pass_GC_threshold)){
  pp_24bp_pass_GC_threshold_Tm <- c(pp_24bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(pp_24bp_pass_GC_threshold[i,]))))
}


#add Tm generated to pp_#bp_pass_GC_threshold
pp_18bp_pass_GC_threshold <-cbind(pp_18bp_pass_GC_threshold, pp_18bp_pass_GC_threshold_Tm)
pp_19bp_pass_GC_threshold <-cbind(pp_19bp_pass_GC_threshold, pp_19bp_pass_GC_threshold_Tm)
pp_20bp_pass_GC_threshold <-cbind(pp_20bp_pass_GC_threshold, pp_20bp_pass_GC_threshold_Tm)
pp_21bp_pass_GC_threshold <-cbind(pp_21bp_pass_GC_threshold, pp_21bp_pass_GC_threshold_Tm)
pp_22bp_pass_GC_threshold <-cbind(pp_22bp_pass_GC_threshold, pp_22bp_pass_GC_threshold_Tm)
pp_23bp_pass_GC_threshold <-cbind(pp_23bp_pass_GC_threshold, pp_23bp_pass_GC_threshold_Tm)
pp_24bp_pass_GC_threshold <-cbind(pp_24bp_pass_GC_threshold, pp_24bp_pass_GC_threshold_Tm)

#rename new column generated "Tm"
colnames(pp_18bp_pass_GC_threshold)[20] <- "Tm"
colnames(pp_19bp_pass_GC_threshold)[21] <- "Tm"
colnames(pp_20bp_pass_GC_threshold)[22] <- "Tm"
colnames(pp_21bp_pass_GC_threshold)[23] <- "Tm"
colnames(pp_22bp_pass_GC_threshold)[24] <- "Tm"
colnames(pp_23bp_pass_GC_threshold)[25] <- "Tm"
colnames(pp_24bp_pass_GC_threshold)[26] <- "Tm"

#filter out the Tm that are not within 50-60 degrees

library(dplyr)

#use the dplyr filter function to filter melting temperatures between 50-60 degrees.The primers that meet this criteria will be 
#placed into a new dataframe called pp_#bp_pass_GC_and_Tm_threshold


#18bp
pp_18bp_pass_GC_and_Tm_threshold <- filter(pp_18bp_pass_GC_threshold, pp_18bp_pass_GC_threshold$Tm >=50, pp_18bp_pass_GC_threshold$Tm <=60)

#19bp
pp_19bp_pass_GC_and_Tm_threshold <- filter(pp_19bp_pass_GC_threshold, pp_19bp_pass_GC_threshold$Tm >=50, pp_19bp_pass_GC_threshold$Tm <=60)

#20bp
pp_20bp_pass_GC_and_Tm_threshold <- filter(pp_20bp_pass_GC_threshold, pp_20bp_pass_GC_threshold$Tm >=50, pp_20bp_pass_GC_threshold$Tm <=60)

#21bp
pp_21bp_pass_GC_and_Tm_threshold <- filter(pp_21bp_pass_GC_threshold, pp_21bp_pass_GC_threshold$Tm >=50, pp_21bp_pass_GC_threshold$Tm <=60)

#22bp
pp_22bp_pass_GC_and_Tm_threshold <- filter(pp_22bp_pass_GC_threshold, pp_22bp_pass_GC_threshold$Tm >=50, pp_22bp_pass_GC_threshold$Tm <=60)

#23bp
pp_23bp_pass_GC_and_Tm_threshold <- filter(pp_23bp_pass_GC_threshold, pp_23bp_pass_GC_threshold$Tm >=50, pp_23bp_pass_GC_threshold$Tm <=60)

#24bp
pp_24bp_pass_GC_and_Tm_threshold <- filter(pp_24bp_pass_GC_threshold, pp_24bp_pass_GC_threshold$Tm >=50, pp_24bp_pass_GC_threshold$Tm <=60)


#############################################################################################
#calculate the Tm of the complementary primer precursors
#workflow is similar to primer precursors Tm calculation


library(TmCalculator)
#
##there are three different ways to calcualte Tm, based off of GC/based off of other nucleotides, and the Wallace Method (good for primer 14-20bp)

##we will use the Tm_NN function which takes into account the bases as well as ionic concentrations of the PCR reaction
#you can change the input of the pcr reaction ionic concentrations using this method
##we will need to convert into strings in order for this function to work


library(Biostrings)

####create empty vector and then run function over the rows
#In the for loop the colon creates a sequence of numbers from 1 to 15, and loops over i = 1, i= 2, .... 

#18bp
cpp_18bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_18bp_pass_GC_threshold)){
  cpp_18bp_pass_GC_threshold_Tm <- c(cpp_18bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_18bp_pass_GC_threshold[i,]))))
}

#19bp
cpp_19bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_19bp_pass_GC_threshold)){
  cpp_19bp_pass_GC_threshold_Tm <- c(cpp_19bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_19bp_pass_GC_threshold[i,]))))
}

#20bp
cpp_20bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_20bp_pass_GC_threshold)){
  cpp_20bp_pass_GC_threshold_Tm <- c(cpp_20bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_20bp_pass_GC_threshold[i,]))))
}

#21bp
cpp_21bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_21bp_pass_GC_threshold)){
  cpp_21bp_pass_GC_threshold_Tm <- c(cpp_21bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_21bp_pass_GC_threshold[i,]))))
}

#22bp
cpp_22bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_22bp_pass_GC_threshold)){
  cpp_22bp_pass_GC_threshold_Tm <- c(cpp_22bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_22bp_pass_GC_threshold[i,]))))
}

#23bp
cpp_23bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_23bp_pass_GC_threshold)){
  cpp_23bp_pass_GC_threshold_Tm <- c(cpp_23bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_23bp_pass_GC_threshold[i,]))))
}

#24bp
cpp_24bp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(cpp_24bp_pass_GC_threshold)){
  cpp_24bp_pass_GC_threshold_Tm <- c(cpp_24bp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(cpp_24bp_pass_GC_threshold[i,]))))
}



#add Tm generated to pp_#bp_pass_GC_threshold
cpp_18bp_pass_GC_threshold <-cbind(cpp_18bp_pass_GC_threshold, cpp_18bp_pass_GC_threshold_Tm)
cpp_19bp_pass_GC_threshold <-cbind(cpp_19bp_pass_GC_threshold, cpp_19bp_pass_GC_threshold_Tm)
cpp_20bp_pass_GC_threshold <-cbind(cpp_20bp_pass_GC_threshold, cpp_20bp_pass_GC_threshold_Tm)
cpp_21bp_pass_GC_threshold <-cbind(cpp_21bp_pass_GC_threshold, cpp_21bp_pass_GC_threshold_Tm)
cpp_22bp_pass_GC_threshold <-cbind(cpp_22bp_pass_GC_threshold, cpp_22bp_pass_GC_threshold_Tm)
cpp_23bp_pass_GC_threshold <-cbind(cpp_23bp_pass_GC_threshold, cpp_23bp_pass_GC_threshold_Tm)
cpp_24bp_pass_GC_threshold <-cbind(cpp_24bp_pass_GC_threshold, cpp_24bp_pass_GC_threshold_Tm)

#rename new column generated "Tm"
colnames(cpp_18bp_pass_GC_threshold)[20] <- "Tm"
colnames(cpp_19bp_pass_GC_threshold)[21] <- "Tm"
colnames(cpp_20bp_pass_GC_threshold)[22] <- "Tm"
colnames(cpp_21bp_pass_GC_threshold)[23] <- "Tm"
colnames(cpp_22bp_pass_GC_threshold)[24] <- "Tm"
colnames(cpp_23bp_pass_GC_threshold)[25] <- "Tm"
colnames(cpp_24bp_pass_GC_threshold)[26] <- "Tm"


#filter out the Tm that are not within 50-60 degrees


library(dplyr)

#use the dplyr filter function to filter melting temperatures between 50-60 degrees.The primers that meet this criteria will be 
#placed into a new dataframe called pp_#bp_pass_GC_and_Tm_threshold


#18bp
cpp_18bp_pass_GC_and_Tm_threshold <- filter(cpp_18bp_pass_GC_threshold, cpp_18bp_pass_GC_threshold$Tm >=50, cpp_18bp_pass_GC_threshold$Tm <=60)

#19bp
cpp_19bp_pass_GC_and_Tm_threshold <- filter(cpp_19bp_pass_GC_threshold, cpp_19bp_pass_GC_threshold$Tm >=50, cpp_19bp_pass_GC_threshold$Tm <=60)

#20bp
cpp_20bp_pass_GC_and_Tm_threshold <- filter(cpp_20bp_pass_GC_threshold, cpp_20bp_pass_GC_threshold$Tm >=50, cpp_20bp_pass_GC_threshold$Tm <=60)

#21bp
cpp_21bp_pass_GC_and_Tm_threshold <- filter(cpp_21bp_pass_GC_threshold, cpp_21bp_pass_GC_threshold$Tm >=50, cpp_21bp_pass_GC_threshold$Tm <=60)

#22bp
cpp_22bp_pass_GC_and_Tm_threshold <- filter(cpp_22bp_pass_GC_threshold, cpp_22bp_pass_GC_threshold$Tm >=50, cpp_22bp_pass_GC_threshold$Tm <=60)

#23bp
cpp_23bp_pass_GC_and_Tm_threshold <- filter(cpp_23bp_pass_GC_threshold, cpp_23bp_pass_GC_threshold$Tm >=50, cpp_23bp_pass_GC_threshold$Tm <=60)

#24bp
cpp_24bp_pass_GC_and_Tm_threshold <- filter(cpp_24bp_pass_GC_threshold, cpp_24bp_pass_GC_threshold$Tm >=50, cpp_24bp_pass_GC_threshold$Tm <=60)


#######################################
#Ange's code: Finding primers and complement primer strands that begin with G/C pairs 
########################################
library(tidyverse)

# Generate a list of primers based on Emma's code 
Allprimerprecursors<-list(pp_18bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_18bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_18bp_pass_GC_and_Tm_threshold))],
                          pp_19bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_19bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_19bp_pass_GC_and_Tm_threshold))],
                          pp_20bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_20bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_20bp_pass_GC_and_Tm_threshold))],
                          pp_21bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_21bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_21bp_pass_GC_and_Tm_threshold))],
                          pp_22bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_22bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_22bp_pass_GC_and_Tm_threshold))],
                          pp_23bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_23bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_23bp_pass_GC_and_Tm_threshold))],
                          pp_24bp_pass_GC_and_Tm_threshold[,c(-ncol(pp_24bp_pass_GC_and_Tm_threshold) + 1, -ncol(pp_24bp_pass_GC_and_Tm_threshold))])

# Exclude primer precursors that begin and end with A or T 
list_of_good_primer_precursors <- list()
for(file in Allprimerprecursors){
  whichcol <-1
  Lastcall <- file[ ,whichcol]
  NotAorT<- !(Lastcall %in% c("A","T"))
  file <-file[NotAorT, ]
  
  
  whichcol <-ncol(file)
  Lastcall <- file[ ,whichcol]
  NotAorT<- !(Lastcall %in% c("A","T"))
  file <-file[NotAorT, ]
  list_of_good_primer_precursors[[ncol(file) - 17]] <- file
  
}

good_primer_precursors <- bind_rows(list_of_good_primer_precursors)


#get the complement sequence based off the forward sequence (change A to T, etc.)
complement_primer_precursors <- as_tibble(apply(good_primer_precursors, 2, chartr, old = "ATGC", new = "TACG"))
#reverse the complement sequence so that it too is written from 5' to 3'
good_complement_primer_precursors <- complement_primer_precursors %>% select(ncol(complement_primer_precursors):1)

#Unite a final list of primer and complement primer precursors that begin and end with G/C pairs 
final_list_of_forward_primers <- good_primer_precursors %>% 
  unite("Primer", X1:X24, na.rm = TRUE, sep = "")

final_list_of_forward_complement_primers <- good_complement_primer_precursors %>% 
  unite("Complment 5-3 Primer", X24:X1, na.rm = TRUE, sep = "")

# Print a final list of good primers and complement primers that begin and end with G/C pairs 
knitr::kable(final_list_of_forward_primers)
knitr::kable(final_list_of_forward_complement_primers)