export LANG=en_US.UTF-8
set -e
#####################
#Dependences:
#bash version 4
#Internet Connection
#Curl
#R (with gtools and xlsx)

#NOTES
#pathoscope: pathoscope and metamix use the tsv extension, don't merge this files in same folder
#sigma: we assume when you worked with sigma, the names of the each fasta folder were the same that fasta gi (consult the script prepare sigmaDB if you don't have this format).
#metaphlan:we working on this
#constrains:metaphlan dependent
#metamix:same pathoscope note
#PERDONAZO METHOD: this method consist in change the mayor reads assigned in specific tax id to a defined permament tax id, getting by consequence the correct analysis when its compare in real data.
#to apply this change the family of tax id must be the same.

#####################################################################################################################
#####################					PARSE PARAMETERS SECTION			#########################################

workband=0
cfileband=0
statusband=0
titogiband=0
tifamilyband=0

for i in "$@"
do
	case $i in
	"--workpath")
		workband=1
	;;
	"--cfile")
		cfileband=1
	;;
	"--tifamily")
		tifamilyband=1
	;;
	"--help")
		echo -e "Options aviable:\n --workpath path where your files are (included requirement files).\n --cfile configuration file."
		echo "requeriment files: makeCSV.R"
		echo -e "\n to apply Perdonazo method, you must especify in the config file the parameter ABSENT=yes, the script automatically calculate corresponding data"
		echo "if ABUNDANCE is missing in the configuration file, its equals to generate results in reads number instead percents"
		exit
	;;
	*)
		if [ $((workband)) -eq 1 ];then
			RUTAINICIAL=$i
			EXECUTIONPATH=`pwd`
			statusband=$((statusband+1))
			workband=0
		fi
		
		if [ $((cfileband)) -eq 1 ];then
			for parameter in `awk '{print}' $i`
			do
				Pname=`echo "$parameter" |awk 'BEGIN{FS="="}{print $1}'`		
				case $Pname in
					"GENOMESIZEBALANCE")
						GENOMESIZEBALANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`
						#echo "${parameters[$i]}"								
					;;
					"COMMUNITYCOMPLEX")
						COMMUNITYCOMPLEX=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"SPECIES")
						SPECIES=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"ABUNDANCE")
						ABUNDANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"DOMINANCE")
					DOMINANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"READSIZE")
						READSIZE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"ABSENT")
						ABSENT=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"METHOD")
						METHOD=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"permanent")
						permanent=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"tipermanent")
						tipermanent=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"CORES")
						CORES=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"THREADS")
						THREADS=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"METASIMFOLDER")
						METASIMFOLDER=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"MUMERPATH")
						MUMERPATH=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"PATHOSCOPEHOME")
						PATHOSCOPEHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"SIGMAHOME")
						SIGMAHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"METAMIXHOME")
						METAMIXHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
				esac
			done
			statusband=$((statusband+1))
			cfileband=0
		fi
		
		if [ $((tifamilyband)) -eq 1 ]; then
			TIFAMILYFILE=$i
			tifamilyband=0
		fi
	;;
	esac
done

if [[ "$ABSENT" =~ "yes" ]] ; then
	if [ "$tipermanent" == ""  ]; then
		echo "ABSENT=yes, but you don't especify the tax id of your permament genome, it's a requisite to apply perdonazo method"
		exit
	else
		curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$tipermanent" > tmp.xml
		FAMILYPERMANENT=`awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="family"){printf "%s,",prev}}' tmp.xml` #family corresponding to fasta permament
		rm tmp.xml
	fi
fi
#####################################################################################################################

###############################					FUNCTION DECLARATION				#################################
function TakeLineageFunction {
	#################################################################################################################################################
	#this function take a tax id from the result files and fetch the lineage to save it in the same file that has already parsed by other functions.
	#################################################################################################################################################

	Matrix=$1
	if [ -f $Matrix ]; then
		#print the headers for the csv
		echo '"Superkingdom","Phylum","Class","Order","Family","Genus","Specie","Name",' >> TaxonomyPredictionMatrix.csv 
		for ti in `awk 'BEGIN{FS="\""}{if(NR>1){print $2}}' $Matrix`
		do
			echo "fetching lineage from ti: $ti"
			 #warning, no error tolerance (I never get the error for cover the case)
			 #fetch the ti by ncbi api
			curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$ti" > tmp.xml
			name=`awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml` #be careful with \n
			lineage=`awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="superkingdom"){printf "%s,",prev}if($3=="phylum"){printf "%s,",prev}if($3=="class"){printf "%s,",prev}if($3=="order"){printf "%s,", prev}if($3=="family"){printf "%s,",prev}if($3=="genus"){printf "%s,",prev}if($3=="species"){printf "%s,",prev}}' tmp.xml`

			echo "$lineage$name," >> TaxonomyPredictionMatrix.csv
			rm tmp.xml
		done
		paste -d '\0' TaxonomyPredictionMatrix.csv $Matrix > tmp.csv 
		rm TaxonomyPredictionMatrix.csv $Matrix
		mv tmp.csv $Matrix
	else
		echo "$Matrix doesn't exist, nothing changes"
	fi

}

function pathoscopeFunction {
	###################################################################################################################################################################################
	#this function take the pathoscope results (tsv file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for tsvfile in `ls -1 *.tsv`
	do
		#if to recognize if the files for post analysis are in reads number or percent (abundance required)
		if [ "$ABUNDANCE" == "" ]; then			
			awk 'BEGIN{FS="|"}{print $2}' $tsvfile |awk -v abu=$ABUNDANCE '{if(NR>2)print $1, $4}' > pathoids.dat
		else
			awk 'BEGIN{FS="|"}{print $2}' $tsvfile |awk -v abu=$ABUNDANCE '{if(NR>2)print $1, ($4/(abu*2))*100}' > pathoids.dat
		fi					
		##########PERDONAZO METHOD FOR ABSENTS##############
		if [ "$ABSENT" == "yes" ]; then
			timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' pathoids.dat`
			#make sure you have tifamily.dat
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
								
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed -i '' "s/[[:<:]]$timayor[[:>:]]/$tipermament/g" pathoids.dat
				echo "--------------------perdonazo in $tsvfile: yes"
			else
				echo "--------------------perdonazo in $tsvfile: no"
			fi
		fi
		########################################
												
		mv pathoids.dat parsed_$tsvfile.dat
		echo "$tsvfile file formated"
	done	
		total=`ls -1 *.tsv.dat |wc -l`
	if [ $((total)) -le 1 ]; then
		echo "need at least 2 files to make a table"
	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat pathoscope_table.csv
		rm parsed* makeCSV.R
		sed -i '' "s/ti.//g" pathoscope_table.csv
		sed -i '' "s/\"\"/\"ti\"/g" pathoscope_table.csv
		TakeLineageFunction pathoscope_table.csv

	fi
}

function metamixFunction {
	###################################################################################################################################################################################
	#this function take the metamix results (tsv file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for tsvfile in `ls -1 *.tsv`
	do
		awk 'BEGIN{FS="\""}{if(NR>1)print $4}' $tsvfile > taxidasigned
				
		#if to recognize if the files for post analysis are in reads number or percent (abundance required)
		if [ "$ABUNDANCE" == "" ]; then
			awk 'BEGIN{FS="\""}{if(NR>1)print $7}' $tsvfile |awk -v abu=$ABUNDANCE  '{print $1}' > readsasigned
		else
			awk 'BEGIN{FS="\""}{if(NR>1)print $7}' $tsvfile |awk -v abu=$ABUNDANCE  '{print ($1/abu)*100}' > readsasigned
		fi		
		paste taxidasigned readsasigned > metamixids.dat
		
		rm taxidasigned readsasigned

		##########PERDONAZO METHOD FOR ABSENTS##############
		if [ "$ABSENT" == "yes" ]; then			
			timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' metamixids.dat`
			#make sure you have tifamily.dat
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
									
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed -i "s/[[:<:]]$timayor[[:>:]]/$tipermament/g" metamixids.dat
				echo "--------------------perdonazo in $tsvfile: yes"
			else
				echo "--------------------perdonazo in $tsvfile: no"
			fi
		fi
		########################################
							
		mv metamixids.dat parsed_$tsvfile.dat
	done		   			
		total=`ls -1 *.tsv.dat |wc -l`
	if [ $((total)) -le 1 ]; then
		echo "need at least 2 files to make a table"
	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat metamix_table.csv
		rm parsed* makeCSV.R
		sed -i '' "s/ti.//g" metamix_table.csv
		sed -i '' "s/\"\"/\"ti\"/g" metamix_table.csv
		TakeLineageFunction metamix_table.csv
	fi
}

function sigmaFunction {
	#####################################################################################################################################################
	#this function take the pathoscope results (gvector.txt file), and parse it to leave only the tax id (fetch ti by gi is necessary in this function)
	#####################################################################################################################################################

	for gvector in `ls -1 *gvector.txt`
	do
		#to recognize if the files for post analysis are in reads number or percent (abundance required)
		if [ "$ABUNDANCE" == "" ]; then
			mappedread=`awk '{if($1=="+"){print $4-2}}' $gvector`
			awk -v map=$mappedread '{if($1=="*"){printf "%d %d\n",$2, ($3*map)/100}}' $gvector > tmp.dat
		else
			awk '{if($1=="*"){print $2, $3}}' $gvector > tmp.dat
		fi
		cp tmp.dat sigmaids.dat
		awk '{print $1}' tmp.dat > tmp2.dat
		for n in `awk '{print $1}' tmp.dat`
		do
			sigmagi=`awk -v target=$n '{if($1=="@"){if($2==target){print $3}}}' $gvector`
			sed -i '' "s/[[:<:]]$n[[:>:]]/$sigmagi/g" tmp2.dat #\< \> for whole word search
		done
		paste tmp.dat tmp2.dat |awk '{print $3, $2}' > sigmaids.dat
		rm tmp*

		##########PERDONAZO METHOD##############

		if [ "$ABSENT" == "yes" ]; then
			gimayor=`awk 'BEGIN{mayor=-1;gi=1}{if($2>mayor){gi=$1;mayor=$2}}END{print gi}' sigmaids.dat` #sigma col 1 have gi 
			timayor=`grep -w "$gimayor" ${RUTAINICIAL}/$TITOGIFILE | awk '{print $1}'`
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
		
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed -i "s/[[:<:]]$gimayor[[:>:]]/$gipermament/g" sigmaids.dat
				echo "--------------------perdonazo in $gvector: yes"
			else
				echo "--------------------perdonazo in $gvector: no"
			fi
		fi
		########################################
									
		#####trade gi x ti#########
		cp sigmaids.dat tmp.dat
		
		for gi in `awk '{print $1}' tmp.dat`	
		do
			ti=""
			while [ "$ti" == "" ]
			do
				ti=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=$gi" |grep "<Id>"|tail -n1 |awk '{print $1}' |cut -d '>' -f 2 |cut -d '<' -f 1`
				#echo "ti: $ti"
			done
			
			sed -i '' "s/[[:<:]]$gi[[:>:]]/$ti/g" sigmaids.dat

		done
					
		echo "$gvector file formated"
		rm tmp.dat
		mv sigmaids.dat parsed_$gvector.dat

	done
	###########################################################
	##############FETCHING TI BY GI############################
	total=`ls -1 *.gvector.txt.dat |wc -l`
	if [ $((total)) -le 1 ]; then
		echo "need at least 2 files to make a table"
	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . gvector.txt.dat sigma_table.csv
		rm parsed* makeCSV.R
		sed -i '' "s/ti.//g" sigma_table.csv
		sed -i '' "s/\"\"/\"ti\"/g" sigma_table.csv
		TakeLineageFunction sigma_table.csv
	fi
	
}

function metaphlanFunction {
	echo "Metaphlan not yet :D"
}

function constrainsFunction {
	echo "Constrains not yet :D"

}

function makeCSV {
echo 'library(xlsx)
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
write.csv(otu_table_noZeroes,outtable)'
}

################################################################################################################
#begin the code
if [ $((statusband)) -ge 2 ]; then
	cd ${RUTAINICIAL}
	for g in $METHOD
	do
		case $g in
			"PATHOSCOPE")
				pathoscopeFunction
			;;
			"METAPHLAN")
				metaphlanFunction
	   		;;
	   		"METAMIX")
	   			metamixFunction
	   		;;
	   		"SIGMA")
				sigmaFunction
	   		;;
	   		"CONSTRAINS")
	   			constrainsFunction
	   		;;
	   		*)
	   			echo "no method aviable for $METHOD"
	   			exit
	   		;;
		esac
	done

else
	echo "Invalid or Missing Parameters, print --help to see the options"
	echo "Usage: bash parseMethods.bash --workpath [files directory] --cfile [config file]"
	exit
fi
