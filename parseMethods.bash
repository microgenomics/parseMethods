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
	"--dbKR")
		dbkrband=1
		invalidband=0
	;;
	"--help")
		echo "Usage: bash parseMethods.bash --workpath . --cfile config"
		echo "Options aviable:"
		echo "--workpath path where your files are"
		echo "--cfile configuration file"
		echo -e "\n Notes:"
		echo "a) If KRAKEN is in your METHODS, Use --dbKR to provide kraken database folder, this is to translate resuts into a tax id and then homologate outputs. "
		echo "b) Make sure you have R (with gtools and xlsx)"
		echo "c) To apply Perdonazo method, you must especify in the config file the parameter ABSENT=YES, the script automatically calculate corresponding data"
		echo "d) If ABUNDANCE is missing in the configuration file, metaphlan and constrains results will write in percent and others in reads number (default), provide a ABUNDANCE will set all results in reads number. This value can be obtained from total lines of your reads files previously used"
		echo "finally, don't use \",\" to name your files, this script will generate a csv, so don't put the character coma in names"
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
		
		if [ $((dbkrband)) -eq 1 ]; then
			dbkrband=0
			if [ -d $i ]; then
				INITIALPATH=`pwd`
				cd $i
				DBKR=`pwd`
				cd $INITIALPATH
			else
				echo "$i file no exist"
				exit
			fi
		fi
	;;
	esac
done


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
			while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]]
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
		awk 'BEGIN{FS="|"}{print $2}' $tsvfile |awk '{if(NR>2)print $1, $4}' > pathoids.dat
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
		sed "s/\"//g" tmp2 > pathoscope_table.csv
		rm  tmp tmp2

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
			sed '1!G;h;$!d' $datfile |awk -v abu=$ABUNDANCE 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print ($2*abu)/100}' | sed '1!G;h;$!d' > quantities
		fi

		sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $3, $6, $9, $12, $15, $18, $22, $18"..."$22}' > sname
		#sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $22}' > sname
		paste -d '_' sname quantities > metaphlanid.dat
		rm sname quantities
		
		########################################
		datfile=`echo "$datfile" |sed "s/,/./g"`
		mv metaphlanid.dat parsed_$datfile.dat
		echo "$datfile file formated"
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
	sed "s/\.\.\./ /g" metaphlan_table.csv > tmp
	sed "s/\"\"/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name/g" tmp > tmp2
	sed "s/\"//g" tmp2 > tmp3
	awk 'BEGIN{FS=","}{gsub("\\.",",",$1);FS=" ";print $0}' tmp3 > metaphlan_table.csv
	rm tmp*
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
		awk '{if(NR>1){print $1, $4}}' $profile > parsed_$profile.dat
		while read line
		do
			genus=`echo "$line" |awk 'BEGIN{FS="_"}{print $1}'`
			species=`echo "$line" |awk 'BEGIN{FS="_| "}{print $2}'`
		
			if [ "$ABUNDANCE" == "" ]; then			
				reads=`echo "$line" |awk '{print $2}'`
			else
				reads=`echo "$line" |awk -v abu=$ABUNDANCE '{print ($2*abu)/100}'`
			fi
		
			lineage=""
			echo "fetching $genus $species lineage"
			while [ "$lineage" == "" ]
			do
				lineage=`curl -s "http://www.ebi.ac.uk/ena/data/view/Taxon:$genus%20$species&display=xml" |awk 'BEGIN{band=0}{if($0~"<lineage>"){band=1;next}if($0~"</lineage>"){band=0};if(band==1){gsub("="," ");print}}' |awk '{toprint="";for(i=1;i<=NF;i++){if($i=="scientificName"){toprint=$(i+1)};if($i=="rank"){toprint=toprint" "$(i+1)}}print toprint}' |tail -r |awk '{gsub("\"","");if($1!="" && $2!=""){if(toprint==""){toprint=$1}else{toprint=toprint" "$1}}}END{print toprint}'`
			done
			echo "$lineage $species $genus...$species""_$reads" >> tmp
		done < <(grep "" parsed_$profile.dat)
		rm parsed_$profile.dat
		mv tmp parsed_$profile.dat

		echo "$profile file formated"

	done

	###########################################################
	##############FETCHING TI BY GI############################
	total=`ls -1 parsed_*.profiles.dat|wc -l`
	if [ $((total)) -le 1 ]; then
		echo "need at least 2 files to make a table"
	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .profiles.dat constrains_table.csv
		rm makeCSV.R parsed*
		sed "s/\.\.\./ /g" constrains_table.csv > tmp
		sed "s/\"\"/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name/g" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		awk 'BEGIN{FS=","}{gsub("\\.",",",$1);FS=" ";print $0}' tmp3 > constrains_table.csv
		rm tmp*
	fi

}

function krakenFunction {
	if [ "$DBKR" == "" ];then
		echo "you must provide the kraken database folder for convertion names to tax id"
		exit
	fi

	for kraken in `ls -1 *.kraken`
	do
		#for % reads
		#totalreads=`wc -l $kraken |awk '{print $1}'`
		#awk -v $totalreads=$totalreads 'BEGIN{FS=";"}{if(NR==FNR && $2!=""){n[$2]+=1}}END{for(key in n){print key";"n[key]/totalreads}}' $kraken > parsed_$kraken
		awk 'BEGIN{FS=";"}{if(NR==FNR && $2!=""){n[$2]+=1}}END{for(key in n){print key";"n[key]}}' $kraken > tmp
		while read line
		do
			name=`echo $line |awk 'BEGIN{FS=";"}{print $1}'`
			reads=`echo $line |awk 'BEGIN{FS=";"}{print $2}'`
			ti=`grep "$name" ${DBKR}/taxonomy/names.dmp |awk '{print $1}'` #this line will get the tax id
			echo "$ti $reads"
		done < <(grep "" tmp) > parsed_$kraken.dat
		rm tmp
		echo "$kraken file formated"	
	done

	#####################################################
	total=`ls -1 *.kraken.dat |wc -l`
	if [ $((total)) -le 1 ]; then
		echo "need at least 2 files to make a table"
	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .kraken.dat kraken_table.csv
		rm parsed* makeCSV.R
		sed "s/ti.//g" kraken_table.csv > tmp
		sed "s/\"\"/ti/g" tmp > tmp2
		sed "s/\"//g" tmp2 > kraken_table.csv
		rm  tmp tmp2

		TakeLineageFunction kraken_table.csv
	fi

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
if(pattr == ".dat.dat" || pattr == ".profiles.dat"){
	read_counts <- lapply(list_of_files, read.table, sep="_", header = FALSE)
}else{
	read_counts <- lapply(list_of_files, read.table, sep=" ", header = FALSE)
#read_counts <- lapply(read_counts, function(x) x[, c(1,2)])
#read_counts <- lapply(read_counts, function(x) x[complete.cases(x),])
}


# for each table make the first col name OTU and the second the patient name

if(pattr == ".dat.dat" || pattr == ".profiles.dat"){
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

if(pattr == ".dat.dat" || pattr == ".profiles.dat"){
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
	   		"KRAKEN")
				krakenFunction
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
