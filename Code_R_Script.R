#import raw data and library
list.of.packages <- c("dplyr", "seqinr", "TmCalculator", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)

#make sure getwd to where the downloaded repository is 
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
primer_precursor_producer <- function(a){
  primer_precursor <- NULL
  for (i in 0:a){
    primer_precursor <- rbind(primer_precursor, plasmid_sequence[(b-i):(b+a-i),], .id = NULL)}
    primer_precursor <- data.frame(primer_precursor, row.names = NULL, 
                                    check.rows = FALSE, check.names = TRUE, 
                                    fix.empty.names = TRUE, stringsAsFactors=FALSE)
  return(primer_precursor)
}

#combine all primers generated with a length range from 18bp to 22 bp
primer_precursors <- bind_rows(primer_precursor_producer(17), primer_precursor_producer(18), 
                               primer_precursor_producer(19), primer_precursor_producer(20),
                               primer_precursor_producer(21), primer_precursor_producer(22),
                               primer_precursor_producer(23), .id = NULL)

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



#############################################Emma Code
#calculate the GC content of each row in pimer_precursors df
#change NA to N, as the functions are already coded to ignore "N" functions in the default settings
######primer_precursors
#change NA to N
primer_precursors[is.na(primer_precursors)] <- "N"

#Calculate GC 
pp_GC <-apply(as.matrix(primer_precursors), MARGIN = 1, FUN = GC)

#add list to the end of the complement_primer_precursors df called "GC content"
primer_precursors["GC Content"] <- pp_GC

#filter out GC content that is not between 40-60, using dplyr
#if the GC content >= 40 and <=60, it will create a new file called primer_precursors_pass_GC_threshold
##output are primers that contain 40-60 (or, 40%-60%) GC content

primer_precursors_pass_GC_threshold <-filter(primer_precursors, primer_precursors$`GC Content` >= 40,primer_precursors$'GC Content' <= 60)


################################################################################################################################################################################################################################
#We will now compute the GC content of the complementary_primer_precursors df
#workflow is similar to primer precursors
#change NA to N
complement_primer_precursors[is.na(complement_primer_precursors)] <- "N"

#Calculate GC
cpp_GC <-apply(as.matrix(complement_primer_precursors), MARGIN = 1, FUN = GC)

#add list to the end of the complement_primer_precursors DF called "GC content"
complement_primer_precursors["GC Content"] <- cpp_GC

#filter out GC content that is not between 40-60
complement_primer_precursors_pass_GC_threshold <-filter(complement_primer_precursors, complement_primer_precursors$`GC Content` >= 40,complement_primer_precursors$`GC Content`<= 60)

################################################################################################################################################################################################################################
#now to calculate the Tm of the primer precursors

###calculate the Tm of a primer through the TmCalculator package

##there are three different ways to calcualte Tm, based off of GC/based off of other nucleotides, and the Wallace Method (good for primer 14-20bp)

##we will use the Tm_NN function which takes into account the bases as well as ionic concentrations of the PCR reaction
#you can change the input of the pcr reaction ionic concentrations using this method
##we will need to convert into strings in order for this function to work

####create empty vector and then run function over the rows
#In the for loop the colon creates a sequence of numbers from 1 to 15, and loops over i = 1, i= 2, .... 

#calculate Tm
pp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(primer_precursors_pass_GC_threshold)){
  pp_pass_GC_threshold_Tm <- c(pp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(primer_precursors_pass_GC_threshold[i,]))))
}

#add list to the end of the complement_primer_precursors_pass_GC_threshold called "Tm"
primer_precursors_pass_GC_threshold["Tm"] <-pp_pass_GC_threshold_Tm

#filter out GC content that is not between 0.4-0.6
primer_precursors_pass_GC_and_Tm_threshold <-filter(primer_precursors_pass_GC_threshold, primer_precursors_pass_GC_threshold$Tm >= 50, primer_precursors_pass_GC_threshold$Tm <=60)


################################################################################################################################################################################################################################
#now to calculate the Tm of the complementary_primer_precursors
#workflow is similar to primer_precursors


#calculate Tm
cpp_pass_GC_threshold_Tm <- c()
for (i in 1:nrow(complement_primer_precursors_pass_GC_threshold)){
  cpp_pass_GC_threshold_Tm <- c(cpp_pass_GC_threshold_Tm, Tm_NN(as.matrix(c2s(complement_primer_precursors_pass_GC_threshold[i,]))))
}

#add list to the end of the complement_primer_precursors_pass_GC_threshold called "Tm"
complement_primer_precursors_pass_GC_threshold["Tm"] <-cpp_pass_GC_threshold_Tm

#filter out GC content that is not between 0.4-0.6 
complement_primer_precursors_pass_GC_and_Tm_threshold <-filter(complement_primer_precursors_pass_GC_threshold, complement_primer_precursors_pass_GC_threshold$Tm >= 50, complement_primer_precursors_pass_GC_threshold$Tm <=60)


#######################################
#Ange's code: Finding primers and complement primer strands that begin with G/C pairs 
########################################

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
