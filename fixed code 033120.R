#####################################################################################################################
##Emma Code for primer precursor, GC
#We will first compute the GC content of the primer precursors generated from 18bp-24bp
#We will use the R package,  seqinr  which is frequently used by biologists for sequence manipulation
install.packages(seqinr)
library(seqinr)

#calculate the GC content of each row in pimer_precursors df
#change NA to N, as the functions are already coded to ignore "N" functions in the default settings
######primer_precursors
#change NA to N
primer_precursors[is.na(primer_precursors)] <- "N"

#Calculate GC 
pp_GC <-apply(as.matrix(primer_precursors), MARGIN = 1, FUN = GC)

#add list to the end of the complement_primer_precursors df called "GC content"
primer_precursors["GC Content"] <- pp_GC

#filter out GC content that is not between 0.4-0.6, using dplyr
#if the GC content >= 0.4 and <=0.6, it will create a new file called primer_precursors_pass_GC_threshold
##output are primers that contain 0.4-0.6 (or, 40%-60%) GC content

primer_precursors_pass_GC_threshold <-filter(primer_precursors, primer_precursors$`GC Content` >= 0.4,primer_precursors$'GC Content' <= 0.6)


################################################################################################################################################################################################################################
#We will now compute the GC content of the complementary_primer_precursors df
#workflow is similar to primer precursors
#change NA to N
complement_primer_precursors[is.na(complement_primer_precursors)] <- "N"

#Calculate GC
cpp_GC <-apply(as.matrix(complement_primer_precursors), MARGIN = 1, FUN = GC)

#add list to the end of the complement_primer_precursors DF called "GC content"
complement_primer_precursors["GC Content"] <- cpp_GC

#filter out GC content that is not between 0.4-0.6
complement_primer_precursors_pass_GC_threshold <-filter(complement_primer_precursors, complement_primer_precursors$`GC Content` >= 0.4,complement_primer_precursors$`GC Content`<= 0.6)

################################################################################################################################################################################################################################
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









