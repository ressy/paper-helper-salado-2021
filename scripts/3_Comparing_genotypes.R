########################################################################
# COMPARING_GENOTYPES.R #
########################################################################

# R script to compare the genotypes generated in Megasat of subsampled fastq
# files with the reference genotypes in order to calculate the number of needed
# reads to obtain a high proportion of correct genotypes.

## Loading required packages:
library(arsenal)
library(ggplot2)

##### 1. Comparison of genotypes tables. ####

## Creating a variable with input directory (with resulting files from Megasat)
## and a variable with output directory (with final csv and graph):

path_input <- 'Your path to directory /Megasat_output with results of Megasat'
dir.create('Your path to the directory /Final_results to save final results (csv and graph)')
path_output <- 'Your path to the directory to save final results (csv and graph) created with dir.create()'

## Creating a df for file with reference genotypes:
reference <- read.delim("Reference_5867genotypes_coverage.csv",
                        header = TRUE, sep = "\t")

## Creating a variable with Megasat output filenames:
files <- list.files(path = path_input, pattern = ".txt",
                    full.names = FALSE, recursive = FALSE)

## Function to compare Megasat ouput with reference:
setwd(path_input)

cmp.lst <- sapply(files, simplify = FALSE, function(x){
  # Read table:
  output_Meg <- read.delim(x, header = TRUE, sep = "\t", stringsAsFactors = F)
  # Order samples in table by PCR ascending:
  output_Meg <- output_Meg[order(output_Meg$Sample_idx1_idx2), ]
  # Replace "X" and "Unscored" by O:
  output_Meg[output_Meg=="X" | output_Meg=="Unscored"] <- 0
  # Converting characters variables in numeric:
  output_Meg[, 2:23] <- sapply(output_Meg[, 2:23], as.integer)
  # Comparing with reference:
  cmp <- comparedf(reference, output_Meg)
  # Calculate proportion of correct genotypes:
  prop <- 1-((n.diffs(cmp)/132))
  # Prop of correct genotypes = 1 - ((N alelos diferentes/N total alelos (132))
  # (6 PCR x 11 usat x 2 alelos)
})

## Transform list into a data frame:
cmp.df <- data.frame(id=names(cmp.lst), prop=unlist(cmp.lst), row.names = NULL)
#str(cmp.df)

## Obtaining part of filename through matching with regular expressions: # Create 2 empty vectors:
Nreads <- c()
Rep <- c()
# Looping through values of list with the filenames:
for (i in files){
  # Strsplit divide the string according to a separator, [[1]][2] select the
  # position of the list, regexpr("\\d+") select only digital numbers and
  # as.numeric convert str in num.
  value <- as.numeric(regmatches(strsplit(i, '-')[[1]][2],
  regexpr("\\d+", strsplit(i, '-')[[1]][2])))
  # uptdate the vector with new data
  Nreads <- c(Nreads, value)
  value2 <- as.numeric(regmatches(strsplit(i, '-')[[1]][3],
  regexpr("\\d+", strsplit(i, '-')[[1]][3])))
  # uptdate the vector with new data
  Rep <- c(Rep, value2)
}

## Joining new variables (Nreads and Rep) to df:
cmp.df <- cbind(cmp.df, Nreads)
cmp.df.final <- cbind(cmp.df, Rep)

## Writing final table:
setwd(path_output)
write.table(cmp.df.final, file = "depth_prop_correct_genotypes.csv",
            sep = "\t", col.names = TRUE, row.names = FALSE)

#### 2. Graph representation. ####

## Creating a summary table with variables: Nreads, Mean, Se, Min, Max:
# Calculating the descriptive variables:
N <- aggregate(prop~Nreads, FUN = length, cmp.df.final)
MEAN <- aggregate(prop~Nreads, FUN = mean, cmp.df.final)
#Se function
se <- function(x) {
  se <- sd(x, na.rm = T)/sqrt(length(x))
  return (se)
}
SE <- aggregate(prop~ Nreads, FUN = se, cmp.df.final)
MIN <- aggregate(prop~ Nreads, FUN = min, cmp.df.final)
MAX <- aggregate(prop~ Nreads, FUN = max, cmp.df.final)
# Joining variables into a table:
TOTAL = merge(N, MEAN, by = 'Nreads')
TOTAL = merge(TOTAL, SE, by = 'Nreads')
TOTAL = merge(TOTAL, MIN, by = 'Nreads')
TOTAL = merge(TOTAL, MAX, by = 'Nreads')
names(TOTAL)[2]<-paste("n")
names(TOTAL)[3]<-paste("mean")
names(TOTAL)[4]<-paste("se")
names(TOTAL)[5]<-paste("min")
names(TOTAL)[6]<-paste("max")
TOTAL
# Including observations of 0 values to start the curve in (0,0) in the graph.
TOTAL <- rbind(c(0,0,0,0,0,0), TOTAL)
TOTAL

## Creating the graph:

# Open pdf file
setwd(path_output)
pdf("depth_rplot.pdf")

# Create the graph
graph <- ggplot(TOTAL, aes(x=Nreads, y=mean)) +
  # setting limits of y axes:
  ylim(0, 1) +
  # setting intervals of x axes
  scale_x_continuous(limits=c(0,200), breaks = c(0,25,50,75,100,150,200)) +
  # marking stablished number of reads in the simulations
  geom_point(size = 3.5, shape = 16) +
  # creating a line to connect points
  geom_line(size = 1) +
  # adding a reference horizontal line for y = 0.95:
  geom_hline(yintercept = 0.9, color = "azure4") +
  # establishing axes titles and main title
  xlab("Number of reads") +
  ylab("Proportion of correct genotypes") +
  # setting size, face and distance for axes text and title:
  theme_bw() +
  theme(axis.text = element_text(color = "grey20", size = 12, face = "plain")) +
  theme(axis.title.y = element_text(color = "black", size = 13, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(color = "black", size = 13, face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0)))

graph

# Close the pdf file
dev.off()
