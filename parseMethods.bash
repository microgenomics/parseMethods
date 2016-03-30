if [[ "$@" =~ "--debug" ]]; then
	set -ex
else
	set -e
fi
#####################
#Dependences:
#bash version 4
#Internet Connection
#Curl
#R (with gtools and xlsx)

#NOTES
#pathoscope: pathoscope and metamix use the tsv extension, so, it will be recognize in form pathoscope<some name>.tsv
#sigma: we assume when you worked with sigma, the names of the each fasta folder were the same that fasta gi (consult the script prepare sigmaDB if you don't have this format).
#metaphlan: the results must have .dat exetension, you can change the actual extension for .dat and the script works anyway
#metamix: pathoscope and metamix use the tsv extension, so, it will be recognize in form metamix<some name>.tsv
#PERDONAZO METHOD: this method consist in change the mayor reads assigned in specific tax id to a defined permament tax id (while they belong the same family), getting by consequence the correct analysis when its compare in real data. This is triggered with ABSENT=YES in config file

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
		echo "Options aviable:"
		echo "--workpath path where your files are"
		echo "--cfile configuration file"
		echo "make sure you have R (with gtools and xlsx)"
		echo -e "\n to apply Perdonazo method, you must especify in the config file the parameter ABSENT=YES, the script automatically calculate corresponding data"
		echo "if ABUNDANCE is missing in the configuration file, metaphlan results will write in percent (default)"
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
					"METHOD")
						METHOD=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"CORES")
						CORES=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"THREADS")
						THREADS=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"PATHOSCOPEHOME")
						PATHOSCOPEHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"SIGMAHOME")
						SIGMAHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"BLASTHOME")
						BLASTHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"METAPHLAN2HOME")
						METAPHLAN2HOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"CONSTRAINSHOME")
						CONSTRAINSHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"SAMTOOLSHOME")
						SAMTOOLSHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`
					;;
					"BOWTIE2HOME")
						BOWTIE2HOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`
					;;
					"KRAKENHOME")
						KRAKENHOME=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`
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

if [[ "$ABSENT" =~ "YES" ]] ; then
	if [ "$tipermanent" == ""  ]; then
		echo "ABSENT=YES, but you don't especify the tax id of your permament genome, it's a requisite to apply perdonazo method"
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
		echo 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Name,' >> TaxonomyPredictionMatrix.csv 
		for ti in `awk 'BEGIN{FS=","}{if(NR>1){print $1}}' $Matrix`
		do
			echo "fetching lineage from ti: $ti"
			 #warning, no error tolerance (I never get the error for cover the case)
			 #fetch the ti by ncbi api
			nofetch=""
			while [ "$nofetch" == "" ]
			do
				curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$ti" > tmp.xml
				nofetch=`cat tmp.xml`
			done
			name=`awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml` #be careful with \n
			lineage=`awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="superkingdom"){printf "%s,",prev}if($3=="phylum"){printf "%s,",prev}if($3=="class"){printf "%s,",prev}if($3=="order"){printf "%s,", prev}if($3=="family"){printf "%s,",prev}if($3=="genus"){printf "%s,",prev}if($3=="species"){printf "%s,",prev}}' tmp.xml`
			cand=`echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "no"}}'`
			if [ "$cand" == "YES" ]; then
				echo "unknow,unknow,unknow,unknow,unknow,unknow,unknow,$name," >> TaxonomyPredictionMatrix.csv
			else
				echo "$lineage$name," >> TaxonomyPredictionMatrix.csv
			fi
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

	for tsvfile in `ls -1 pathoscope*.tsv`
	do
		#if to recognize if the files for post analysis are in reads number or percent (abundance required)
		awk 'BEGIN{FS="|"}{print $2}' $tsvfile |awk '{if(NR>2)print $1, $4}' > pathoids.dat
				
		##########PERDONAZO METHOD FOR ABSENTS##############
		if [ "$ABSENT" == "YES" ]; then
			timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' pathoids.dat`
			#make sure you have tifamily.dat
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
								
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed "s/[[:<:]]$timayor[[:>:]]/$tipermament/g" pathoids.dat > tmp
				rm pathoids.dat
				mv tmp pathoids.dat
				echo "--------------------perdonazo in $tsvfile: YES"
			else
				echo "--------------------perdonazo in $tsvfile: no"
			fi
		fi
		########################################
		tsvfile=`echo "$tsvfile" |sed "s/,/./g"`
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
		sed "s/ti.//g" pathoscope_table.csv > tmp
		sed "s/\"\"/\"ti\"/g" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		rm pathoscope_table.csv tmp tmp2
		mv tmp3 pathoscope_table.csv

		TakeLineageFunction pathoscope_table.csv

	fi
}
function metaphlanFunction {

for datfile in `ls -1 *.dat`
do
		#if to recognize if the files for post analysis are in reads number or percent (abundance required)
		if [ "$ABUNDANCE" == "" ]; then			
			sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print $2}' | sed '1!G;h;$!d' > quantities
		else
			#we assuming that abundance is by paired end reads
			sed '1!G;h;$!d' $datfile |awk -v abu=$ABUNDANCE 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print ($2*(abu*2))/100}' | sed '1!G;h;$!d' > quantities
		fi

		sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $3, $6, $9, $12, $15, $18, $22, $18"_"$22}' > sname
		#sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $22}' > sname
		paste -d '_' sname quantities > metaphlanid.dat
		rm sname quantities
		
		########################################
		datfile=`echo "$datfile" |sed "s/,/./g"`
		mv metaphlanid.dat parsed_$datfile.dat
		echo "$datfile file formated"


		if [ "$ABSENT" == "YES" ]; then			
			family=`awk 'BEGIN{mayor=-1;fam=""}{gsub("_"," ");if($8>mayor){fam=$5;mayor=$8}}END{print fam}' parsed_$datfile.dat`
									
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed "s/[[:<:]]$timayor[[:>:]]/$tipermament/g" parsed_$datfile.dat > tmp
				rm parsed_$datfile.dat
				mv tmp parsed_$datfile.dat
				echo "--------------------perdonazo in parsed_$datfile.dat: YES"
			else
				echo "--------------------perdonazo in $parsed_$datfile.dat: no"
			fi
		fi
done

total=`ls -1 *.dat.dat |wc -l`
if [ $((total)) -le 1 ]; then
	echo "need at least 2 files to make a table"
else
	#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
	#parameters: work_directory pattern_file name_out_table
	makeCSV > makeCSV.R
	Rscript makeCSV.R . .dat.dat metaphlan_table.csv
	rm makeCSV.R parsed*
	sed "s/\"\"/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name/g" metaphlan_table.csv > tmp
	sed "s/\"//g" tmp > metaphlan_table.csv
	rm tmp
	awk 'BEGIN{FS=","}{gsub(/\./,",",$1);gsub(" ",",",$0);print $0}' metaphlan_table.csv > tmpcsv
	rm metaphlan_table.csv
	mv tmpcsv metaphlan_table.csv
fi
}
function metamixFunction {
	###################################################################################################################################################################################
	#this function take the metamix results (tsv file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for tsvfile in `ls -1 metamix*.tsv`
	do
		awk 'BEGIN{FS="\""}{if(NR>1)print $4}' $tsvfile > taxidasigned
				
			awk 'BEGIN{FS="\""}{if(NR>1)print $7}' $tsvfile |awk '{print $1}' > readsasigned

		paste taxidasigned readsasigned > metamixids.dat
		
		rm taxidasigned readsasigned

		##########PERDONAZO METHOD FOR ABSENTS##############
		if [ "$ABSENT" == "YES" ]; then			
			timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' metamixids.dat`
			#make sure you have tifamily.dat
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
									
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed "s/[[:<:]]$timayor[[:>:]]/$tipermament/g" metamixids.dat > tmp
				rm metamixids.dat
				mv tmp metamixids.dat
				echo "--------------------perdonazo in $tsvfile: YES"
			else
				echo "--------------------perdonazo in $tsvfile: no"
			fi
		fi
		########################################
		tsvfile=`echo "$tsvfile" |sed "s/,/./g"`
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
		sed "s/ti.//g" metamix_table.csv > tmp
		sed "s/\"\"/\"ti\"/g" tmp > tmp2
		sed "s/\"//g" tmp2 > metamix_table.csv
		rm tmp tmp2
		TakeLineageFunction metamix_table.csv
	fi
}

function sigmaFunction {
	#####################################################################################################################################################
	#this function take the sigma results (gvector.txt file), and parse it to leave only the tax id (fetch ti by gi is necessary in this function)
	#####################################################################################################################################################

	for gvector in `ls -1 *gvector.txt`
	do
		#to recognize if the files for post analysis are in reads number or percent (abundance required)
			mappedread=`awk '{if($1=="+"){print $4-2}}' $gvector`
			awk -v map=$mappedread '{if($1=="*"){printf "%d %d\n",$2, ($3*map)/100}}' $gvector > tmp.dat

		##########PERDONAZO METHOD##############

		if [ "$ABSENT" == "YES" ]; then
			gimayor=`awk 'BEGIN{mayor=-1;gi=1}{if($2>mayor){gi=$1;mayor=$2}}END{print gi}' sigmaids.dat` #sigma col 1 have gi 
			timayor=`grep -w "$gimayor" ${RUTAINICIAL}/$TITOGIFILE | awk '{print $1}'`
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
		
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed "s/[[:<:]]$gimayor[[:>:]]/$gipermament/g" sigmaids.dat > tmp
				rm sigmaids.dat
				mv tmp sigmaids.dat
				echo "--------------------perdonazo in $gvector: YES"
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
			
			sed "s/[[:<:]]$gi[[:>:]]/$ti/g" sigmaids.dat > tmp
			rm sigmaids.dat
			mv tmp sigmaids.dat

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
		sed "s/ti.//g" sigma_table.csv > tmp
		sed "s/\"\"/\"ti\"/g" tmp > tmp2
		sed "s/\"//g" tmp2 > sigma_table.csv
		rm tmp tmp2
		TakeLineageFunction sigma_table.csv
	fi
	
}


function constrainsFunction {
	
	for profile in `ls -1 *.profiles`
	do
		awk '{if(NR>1){print $1, $4}}' $profile > parsed_$profile
		echo 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Name,' > constrains_table.csv
		while read line
		do
			tosearch=`echo $line |awk '{gsub("_"," ");print $1"%20"$2}'`
			percentreads=`echo $line |awk '{print $2}'`
			lineage=`curl -s "http://www.ebi.ac.uk/ena/data/view/Taxon:$tosearch&display=xml" |awk 'BEGIN{band=0}{if($0~"<lineage>"){band=1;next};if($0~"</lineage>"){band=0};if(band==1){gsub("\"","");gsub("="," ");print $7, $3}}' |tail -r|awk '{if(NR<=2){next};lineage=lineage$2" "}END{print lineage}'`							

		done < <(grep "" parsed_$profile)
		##GETTING LINEAGE
		
		
		##########PERDONAZO METHOD##############

		if [ "$ABSENT" == "YES" ]; then
			gimayor=`awk 'BEGIN{mayor=-1;gi=1}{if($2>mayor){gi=$1;mayor=$2}}END{print gi}' sigmaids.dat` #sigma col 1 have gi 
			timayor=`grep -w "$gimayor" ${RUTAINICIAL}/$TITOGIFILE | awk '{print $1}'`
			family=`grep "$timayor" ${RUTAINICIAL}/$TIFAMILYFILE | awk '{print $2}'`
		
			if [ "$family" == "$FAMILYPERMANENT" ]; then
				sed "s/[[:<:]]$gimayor[[:>:]]/$gipermament/g" sigmaids.dat
				echo "--------------------perdonazo in $gvector: YES"
			else
				echo "--------------------perdonazo in $gvector: no"
			fi
		fi
		########################################

					
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
		sed -i '' "s/\"//g" sigma_table.csv
		TakeLineageFunction sigma_table.csv
	fi

}

function krakenFunction {
	echo "Kraken not yet :D"

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
  patient_names[i]<-substr(list_of_files[i],1,nchar(list_of_files[i])-4)
}

# read in each table

#read_counts <- lapply(list_of_files, read.table, sep="\t", header = FALSE, skip =2)
if(pattr == ".dat.dat"){
	read_counts <- lapply(list_of_files, read.table, sep="_", header = FALSE)
}else{
	read_counts <- lapply(list_of_files, read.table, sep=" ", header = FALSE)
#read_counts <- lapply(read_counts, function(x) x[, c(1,2)])
#read_counts <- lapply(read_counts, function(x) x[complete.cases(x),])
}


# for each table make the first col name OTU and the second the patient name

if(pattr == ".dat.dat"){
	for( i in 1:length(list_of_files)){
  		colnames(read_counts[[i]])<- c(patient_names[i])
  		#print(read_counts[i])
	}
}else{
	for( i in 1:length(list_of_files)){
  		colnames(read_counts[[i]])<- c("ti", patient_names[i])
	}
}

# list of lists called otu which stores the first column otu names for each dataframe
otu<-NULL

if (pattr == ".dat.dat"){
	for( i in 1:length(list_of_files)){
	name<-paste(as.character(read_counts[[i]][,1]))
  	otu[i]<- list(name)
  	#print(otu[i])
}

}else{
	for( i in 1:length(list_of_files)){
	name<-paste("ti",as.character(read_counts[[i]][,1]))
    otu[i]<- list(name)
	}
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
#print(otu_table_noZeroes[,1])
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
