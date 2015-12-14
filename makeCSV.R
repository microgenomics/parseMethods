library(xlsx)
library(gtools)

args <-commandArgs()

directory<-args[6]
pattr<-args[7]
outtable<-args[8]

setwd(directory)

#filename<-args[6]
#filedata<-args[7]


list_of_files <- list.files(path=directory, pattern = pattr)
# make a list of just the patient names

patient_names<-NULL

for( i in 1:length(list_of_files)){
#  patient_names[i]<-substring(list_of_files[i], 1, 6)
  patient_names[i]<-substr(list_of_files[i],8,nchar(list_of_files[i])-15)
}

# read in each table

#read_counts <- lapply(list_of_files, read.table, sep="\t", header = FALSE, skip =2)
read_counts <- lapply(list_of_files, read.table, sep=" ", header = FALSE)
#read_counts <- lapply(read_counts, function(x) x[, c(1,2)])
#read_counts <- lapply(read_counts, function(x) x[complete.cases(x),])

# for each table make the first col name OTU and the second the patient name

for( i in 1:length(list_of_files)){
  colnames(read_counts[[i]])<- c("ti", patient_names[i])
}
#print(read_counts[1])

# list of lists called otu which stores the first column otu names for each dataframe
otu<-NULL
for( i in 1:length(list_of_files)){
#	print(read_counts[[1]][1])
	name<-paste("ti",as.character(read_counts[[i]][,1]))
  otu[i]<- list(name)
}

# for each dataframe in read_counts transpose and then 

read_counts <- lapply(read_counts, function(x) t(x[,2]))

# add the otus back as the column name

for( i in 1:length(list_of_files)){
  read_counts[[i]]<-data.frame(read_counts[[i]])
#  print(read_counts[[i]])
	#print(otu[i])
  colnames(read_counts[[i]])<-otu[[i]]
  #print(read_counts[i])
  read_counts[[i]]<-data.frame(patient = patient_names[i], read_counts[[i]])
#print(paste("reaaad:",read_counts[i]))
}

# combine the different dataframes together
otu_table <- read_counts[[1]]
for( i in 2:length(list_of_files)){
  otu_table <- smartbind(otu_table, read_counts[[i]], fill = 0)
}

# transpose the table back so that the microbes are the rows and the patients are the col

otu_table<-t(data.matrix(otu_table))
colnames(otu_table)<-patient_names
otu_table<-otu_table[2:nrow(otu_table),]

# remove zeroes
otu_table_noZeroes<-otu_table[apply(otu_table, 1, function(x){ !isTRUE(all.equal(sum(x),0))}),]
write.csv(otu_table_noZeroes,outtable)