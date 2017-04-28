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
#all files fetch for species name

#####################################################################################################################
#####################					PARSE PARAMETERS SECTION			#########################################

workband=0
cfileband=0
statusband=0
deepsband=0
methodband=0

for i in "$@"
do
	case $i in
	"--workpath")
		workband=1
	;;
	"--cfile")
		cfileband=1
	;;
	"--method")
		methodband=1
	;;
	"--help")
		echo "Usage: bash parseMethods.bash --workpath . --cfile config"
		echo "Options aviable:"
		echo "--workpath path where your files are"
		echo "--cfile configuration file"
		echo "--deeps to convert the % of sigma, metaphlan2 and constrains results into raw abundance reads"
		echo -e "\n Notes:"
		echo "* Make sure you have R (with gtools and xlsx)"
		echo "* don't use ',' to name your files, this script will generate a csv and ',' may cause problems"
		echo "DEPENDENCES: stable internet connection"
		exit
	;;
	*)
		if [ $((workband)) -eq 1 ];then
			RUTAINICIAL=$i
			EXECUTIONPATH=$(pwd)
			statusband=$((statusband+1))
			workband=0
		fi
		
		if [ $((cfileband)) -eq 1 ];then
			for parameter in $(awk '{print}' $i)
			do
				Pname=$(echo "$parameter" |awk -F"=" '{print $1}')	
				case $Pname in
					"DEPTH")
						DEEPS=$(echo "$parameter" | awk -F"=" '{print $2}')					
					;;
					"METHOD")
						METHOD=$(echo "$parameter" | awk -F"=" '{print $2}' | sed "s/,/ /g")				
					;;
			
				esac
			done
			statusband=$((statusband+1))
			cfileband=0
		fi

		if [ $((methodband)) -eq 1 ];then
			METHOD=$(echo "$i" |sed "s/,/ /g")
			statusband=$((statusband+1))
			methodband=0

		fi

	;;
	esac
done


###############################					FUNCTION DECLARATION				#################################
function groupingFunction {
	echo 'args<-commandArgs()
	file<-c(args[6])
	headr<-c(args[7])
	if(headr=="T"){
			df<-read.csv(file, header = T, check.names = F)
			newdf<-aggregate(. ~ Name + Species + Genus + Family + Order + Class + Phylum + Kingdom, df, FUN = sum)
      		newdf<-newdf[,c(8,7,6,5,4,3,2,1,9:length(colnames(newdf)))]

	}else{
			df<-read.csv(file, header = F)
			colnames(df)<-c("COL1","COL2")
			newdf<-aggregate(. ~ COL1, df, FUN = sum)
	}
	write.csv(newdf,file,row.names = F,quote = F)' > grp.R
}
function TakeLineageFunction {
	#################################################################################################################################################
	#this function take a tax id from the result files and fetch the lineage to save it in the same file that has already parsed by other functions.
	#################################################################################################################################################

	Matrix=$1
	if [ -f $Matrix ]; then
		#print the headers for the csv
		echo 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Name,' >> TaxonomyPredictionMatrix.csv 
		for ti in $(awk -F"," '{if(NR>1){gsub("ti\.","");print $1}}' $Matrix)
		do
			echo "fetching lineage from ti: $ti"
			 #warning, no error tolerance (I never get the error for cover the case)
			 #fetch the ti by ncbi api
			nofetch=""
			while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]] || [[ "$nofetch" =~ "Bad Gateway!" ]]
			do
				if curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$ti" > tmp.xml ;then
					touch tmp.xml
					nofetch=$(cat tmp.xml)
				else
					echo "curl error fetch, internet connection?"
					nofetch=""
				fi
			done
			name=$(awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml) #be careful with \n
			spctoawk=$(awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml |awk '{print $2}')
			lineage=$(awk -v emergencyname=$spctoawk 'BEGIN{FS="[<|>]";prev="";superk="";phy="";class="";order="";fam="";gen="";spc=""}{if($2=="ScientificName"){prev=$3}if($3=="superkingdom"){superk=prev}if($3=="phylum"){phy=prev}if($3=="class"){class=prev}if($3=="order"){order=prev}if($3=="family"){fam=prev}if($3=="genus"){gen=prev}if($3=="species"){spc=prev}}
			END{if(superk==""){printf "unknown,"}else{printf "%s,",superk};if(phy==""){printf "unknow,"}else{printf "%s,",phy}; if(class==""){printf "unknow,"}else{printf "%s,",class}; if(order==""){printf "unknow,"}else{printf "%s,",order}; if(fam==""){printf "unknow,"}else{printf "%s,",fam}; if(gen==""){printf "unknow,"}else{printf "%s,",gen}; if(spc==""){if(emergencyname==""){print "unknow,"}else{printf "%s,",emergencyname}}else{printf "%s,",spc}}' tmp.xml)
			cand=$(echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}')
			lineage=$(echo $lineage |awk '{gsub("\\[|\\]","");print $0}')
			if [ "$cand" == "YES" ]; then
				newname=$(echo $name |awk '{print $2, $3}') #be careful with $name, maybe is not made by 3 cols (Candidatus somegenus somespecies)
				echo "unknow,unknow,unknow,unknow,unknow,unknow,$newname,$name," >> TaxonomyPredictionMatrix.csv
			else
				name=$(echo "$name" |awk '{print $1, $2}')
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
	awk '{gsub("\\[|\\]","");print}' $Matrix > tmp
	rm $Matrix
	#only get by species
	awk -F"," '{for(i=1;i<NF;i++){if(i!=9){printf "%s,",$i}}printf "%s\n",$NF}' tmp > $Matrix
	rm tmp

	echo "Grouping otus simulated otus"
	groupingFunction
	Rscript grp.R $Matrix T
	rm grp.R
}
function pathoscopeFunction {
	###################################################################################################################################################################################
	#this function take the pathoscope results (tsv file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for tsvfile in $(ls -1 *report.tsv)
	do
		awk 'BEGIN{FS="|"}{print $2}' $tsvfile |awk '{if(NR>2)print $1, $4/2}' > pathoids.dat
		tsvfile=$(echo "$tsvfile" |sed "s/,/./g")
		mv pathoids.dat parsed_$tsvfile.dat
		echo "$tsvfile file formated"
	done
		total=$(ls -1 *.tsv.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$tsvfile.dat parsed__$tsvfile.dat
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat pathoscope_table.csv
		rm parsed* makeCSV.R

		TakeLineageFunction pathoscope_table.csv
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' pathoscope_table.csv > tmp && mv tmp pathoscope_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat pathoscope_table.csv
		rm parsed* makeCSV.R

		TakeLineageFunction pathoscope_table.csv

	fi
}
function metaphlanFunction {

	for datfile in $(ls -1 metaphlan*.dat)
	do
			#if to recognize if the files for post analysis are in reads number or percent (abundance required)

			#if [ "$READTYPE" == "PAIRED" ]; then
			#	sed '1!G;h;$!d' $datfile |awk -v abu=$DEEPS 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print (($2*abu)*2)/100}' | sed '1!G;h;$!d' > quantities
			#else
				sed '1!G;h;$!d' $datfile |awk -v abu=$DEEPS 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print ($2*abu)/100}' | sed '1!G;h;$!d' > quantities
			#fi

			sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk -F"|" '{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $3, $6, $9, $12, $15, $18, $22, $18"..."$22}' > sname
			#sed '1!G;h;$!d' $datfile |awk 'BEGIN{sum=0}{sum+=$2;if(sum<=100)print}' | sed '1!G;h;$!d' |awk '{print $1}' |awk 'BEGIN{FS="|"}{print $1, $2, $3, $4, $5, $6, $7}' |awk 'BEGIN{FS="_| "}{print $22}' > sname
			paste -d '_' sname quantities > metaphlanid.dat
			rm sname quantities
			
			########################################
			datfile=$(echo "$datfile" |sed "s/,/./g")
			mv metaphlanid.dat parsed_$datfile.dat
			echo "$datfile file formated"
	done

	total=$(ls -1 *.dat.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$datfile.dat parsed__$datfile.dat

		makeCSV > makeCSV.R
		Rscript makeCSV.R . .dat.dat metaphlan_table.csv
		rm makeCSV.R parsed*
		sed "s/\.\.\./ /g" metaphlan_table.csv > tmp
		sed "s/,parsed/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name,parsed/" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		awk 'BEGIN{FS=","}{if(NR==1){gsub("\\.",",",$1);FS=" ";gsub(" ",",",$0);print $0}else{gsub("\\.",",",$1);FS=" ";print $0}}' tmp3 > tmp4
		awk -F"," '{if(NR==1){print $0;next}for(i=1;i<7;i++){printf "%s,",$i};printf "%s,%s,",$8,$8; for(i=9;i<NF;i++){printf "%s,",$i}; printf "%s\n",$NF}' tmp4 > metaphlan_table.csv

		rm tmp*
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' metaphlan_table.csv > tmp && mv tmp metaphlan_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .dat.dat metaphlan_table.csv
		rm makeCSV.R parsed*
		sed "s/\.\.\./ /g" metaphlan_table.csv > tmp
		sed "s/,parsed/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name,parsed/" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		awk 'BEGIN{FS=","}{if(NR==1){gsub("\\.",",",$1);FS=" ";gsub(" ",",",$0);print $0}else{gsub("\\.",",",$1);FS=" ";print $0}}' tmp3 > tmp4
		awk -F"," '{if(NR==1){print $0;next}for(i=1;i<7;i++){printf "%s,",$i};printf "%s,%s,",$8,$8; for(i=9;i<NF;i++){printf "%s,",$i}; printf "%s\n",$NF}' tmp4 > metaphlan_table.csv

		rm tmp*
	fi
}
function metamixFunction {
	###################################################################################################################################################################################
	#this function take the metamix results (tsv file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for tsvfile in $(ls -1 metamix*.tsv)
	do
		awk 'BEGIN{FS="\""}{if(NR>1)print $4}' $tsvfile > taxidasigned
				
		awk 'BEGIN{FS="\""}{if(NR>1)print $7/2}' $tsvfile |awk '{print $1}' > readsasigned

		paste -d " " taxidasigned readsasigned > metamixids.dat
		delete=$(grep "unknown" -n metamixids.dat |awk -F ":" '{print $1}')

		if [ "$delete" != "" ]; then
			sed "${delete}d" metamixids.dat > tmp
			rm metamixids.dat
			mv tmp metamixids.dat
		fi

		rm taxidasigned readsasigned

		tsvfile=$(echo "$tsvfile" |sed "s/,/./g")
		mv metamixids.dat parsed_$tsvfile.dat
	done		   

	total=$(ls -1 *.tsv.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$tsvfile.dat parsed__$tsvfile.dat
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat metamix_table.csv
		rm parsed* makeCSV.R
		sed "s/ti.//g" metamix_table.csv > tmp
		rm metamix_table.csv && mv tmp metamix_table.csv
		TakeLineageFunction metamix_table.csv
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' metamix_table.csv > tmp && mv tmp metamix_table.csv


	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tsv.dat metamix_table.csv
		rm parsed* makeCSV.R
		sed "s/ti.//g" metamix_table.csv > tmp
		rm metamix_table.csv && mv tmp metamix_table.csv
		TakeLineageFunction metamix_table.csv
	fi
}
function sigmaFunction {
	#####################################################################################################################################################
	# This function take the sigma results (gvector.txt file), and parse it to leave only the tax id (fetch ti by gi is necessary in this function)
	# SIGMA NEEDS THE DEEPS TO CONVERT % INTO RAW DATA, PROVIDE IT BY --SigmaDeepSeq
	#####################################################################################################################################################

	for gvector in $(ls -1 *gvector.txt)
	do
		#to recognize if the files for post analysis are in reads number or percent (abundance required)
		#mappedread=`awk '{if($1=="+"){print $4-2}}' $gvector`
		#awk -v map=$mappedread '{if($1=="*"){printf "%d %d\n",$2, ($3*map)/100}}' $gvector > tmp.dat
		awk '{if($1=="@"){print $2, $3}}' $gvector > index.dat
		awk '{if($1=="*"){print $2, $3}}' $gvector > ids.dat

		awk '{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $2, n[$1]}}}' ids.dat index.dat > sigmaids.dat

		cp sigmaids.dat tmp.dat
		#####trade gi x ti#########
		cat tmp.dat |while read line
		do
			gi=$(echo $line |awk '{print $1}')
			abu=$(echo $line |awk -v deep=$DEEPS '{print $2*deep/100}')	

			ti=""
			while [ "$ti" == "" ]
			do
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=$gi" |grep "<Id>"|tail -n1 |awk '{print $1}' |cut -d '>' -f 2 |cut -d '<' -f 1)
				#echo "ti: $ti"
			done
			if [ "$ti" == "$gi" ];then
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$gi" |head -n20 |grep "id" |awk '{print $2}' |head -n 1)				
			fi
		
			echo "$ti $abu"
		
		done > sigmaids.dat

		echo "$gvector file formated"
		rm tmp.dat index.dat ids.dat
		mv sigmaids.dat parsed_$gvector.dat

	done
	###########################################################
	##############FETCHING TI BY GI############################
	total=$(ls -1 *.gvector.txt.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$gvector.dat parsed__$gvector.dat
		makeCSV > makeCSV.R
		Rscript makeCSV.R . gvector.txt.dat sigma_table.csv
		rm parsed* makeCSV.R
		sed "s/ti.//g" sigma_table.csv > tmp
		rm sigma_table.csv && mv tmp sigma_table.csv
		TakeLineageFunction sigma_table.csv
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' sigma_table.csv > tmp && mv tmp sigma_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . gvector.txt.dat sigma_table.csv
		rm parsed* makeCSV.R
		sed "s/ti.//g" sigma_table.csv > tmp
		rm sigma_table.csv && mv tmp sigma_table.csv
		TakeLineageFunction sigma_table.csv
	fi	
}
function constrainsFunction {
	
	for profile in $(ls -1 *.profiles)
	do
		awk '{if(NR>1){print $1, $4}}' $profile > parsed_$profile.dat
		rm -f tmp
		while read line
		do
			genus=$(echo "$line" |awk 'BEGIN{FS="_"}{print $1}')
			species=$(echo "$line" |awk 'BEGIN{FS="_| "}{print $2}')
		
			#if [ "$DEEPS" == "" ]; then			
			#	reads=$(echo "$line" |awk '{print $2}')
			#else
			#	if [ "$READTYPE" == "PAIRED" ]; then			
			#		reads=$(echo "$line" |awk -v abu=$DEEPS '{print ($2*abu*2)/100}')
			#	else
					reads=$(echo "$line" |awk -v abu=$DEEPS '{print ($2*abu)/100}')
			#	fi
			#fi
		
			lineage=""
			echo "fetching $genus $species lineage"
			while [ "$lineage" == "" ]
			do
				lineage=$(curl -s "https://www.ebi.ac.uk/ena/data/view/Taxon:$genus%20$species&display=xml" |awk 'BEGIN{band=0}{if($0~"<lineage>"){band=1;next}if($0~"</lineage>"){band=0};if(band==1){gsub("="," ");print}}' |awk '{toprint="";for(i=1;i<=NF;i++){if($i=="scientificName"){toprint=$(i+1)};if($i=="rank"){toprint=toprint" "$(i+1)}}print toprint}' |tail -r |awk '{gsub("\"","");if($1!="" && $2!=""){if(toprint==""){toprint=$1}else{toprint=toprint" "$1}}}END{print toprint}' |awk '{gsub("/","");print $1, $2, $3, $4, $5}')
			done

			cand=$(echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}')
			if [ "$cand" == "YES" ]; then
				echo "unknow,unknow,unknow,unknow,unknow,unknow,$genus $species,$genus...$species""_$reads" >> tmp
			else
				echo "$lineage $genus $species $genus...$species""_$reads" >> tmp
			fi
		done < <(grep "" parsed_$profile.dat)
		rm parsed_$profile.dat
		mv tmp parsed_$profile.dat

		echo "$profile file formated"
	done
	###########################################################
	##############FETCHING TI BY GI############################
	total=$(ls -1 parsed_*.profiles.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$profile.dat parsed__$profile.dat

		makeCSV > makeCSV.R
		Rscript makeCSV.R . .profiles.dat constrains_table.csv
		rm makeCSV.R parsed*

		sed "s/\.\.\./ /g" constrains_table.csv > tmp
		awk '{FS=",";if(NR==1){printf "Kingdom,Phylum,Class,Order,Family,Genus,Species,Name%s\n",$0;FS=" "}else{FS=",";gsub("\\.",",",$1);FS=" ";print $0}}' tmp > constrains_table.csv
		rm tmp*
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' constrains_table.csv > tmp && mv tmp constrains_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .profiles.dat constrains_table.csv
		rm makeCSV.R parsed*

		sed "s/\.\.\./ /g" constrains_table.csv > tmp
		awk '{FS=",";if(NR==1){printf "Kingdom,Phylum,Class,Order,Family,Genus,Species,Name%s\n",$0;FS=" "}else{FS=",";gsub("\\.",",",$1);FS=" ";print $0}}' tmp > constrains_table.csv
		rm tmp*
	fi
}
function krakenFunction {

	for kraken in $(ls -1 *.kraken)
	do
		#-1 for the root line
		#NF contain numbers
		#$1 reads number
		totallines=$(wc -l $kraken |awk '{print $1-1}')
		awk -F "d__|p__|c__|o__|f__|g__|s__" -v total=$totallines '{if(NR<total){print $NF, $1}}' $kraken |awk '{print $1, $2}' |awk '{gsub("_"," ");if(NF==3)print}' > parsed_$kraken.dat
		rm -f tmp
		while read line
		do
			genus=$(echo "$line" |awk '{print $1}')
			species=$(echo "$line" |awk '{print $2}')
			#KRAKEN DOESN'T NEED ADJUST DEEPS
			reads=$(echo "$line" |awk '{print $3}')
		
			lineage=""
			echo "fetching $genus $species lineage"
			while [ "$lineage" == "" ]
			do
				lineage=$(curl -s "https://www.ebi.ac.uk/ena/data/view/Taxon:$genus%20$species&display=xml" |awk 'BEGIN{band=0}{if($0~"<lineage>"){band=1;next}if($0~"</lineage>"){band=0};if(band==1){gsub("="," ");print}}' |awk '{toprint="";for(i=1;i<=NF;i++){if($i=="scientificName"){toprint=$(i+1)};if($i=="rank"){toprint=toprint" "$(i+1)}}print toprint}' |tail -r |awk '{gsub("\"","");if($1!="" && $2!=""){if(toprint==""){toprint=$1}else{toprint=toprint" "$1}}}END{print toprint}' |awk '{gsub("/","");gsub("\\[|\\]","");print $1, $2, $3, $4, $5}')
			done

			cand=$(echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}')
			if [ "$cand" == "YES" ]; then
				echo "unknow,unknow,unknow,unknow,unknow,$genus $species,$genus...$species""_$reads" >> tmp
			else
				echo "$lineage $genus $species $genus...$species""_$reads" >> tmp
			fi
		done < <(grep "" parsed_$kraken.dat)
		rm parsed_$kraken.dat
		mv tmp parsed_$kraken.dat
		echo "$kraken file formated"
	done

	#####################################################
	total=$(ls -1 *.kraken.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$kraken.dat parsed__$kraken.dat
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .kraken.dat kraken_table.csv
		rm parsed* makeCSV.R
		sed "s/\.\.\./ /g" kraken_table.csv > tmp
		sed "s/,parsed/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name,parsed/" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		awk 'BEGIN{FS=","}{if(NR==1){gsub("\\.",",",$1);FS=" ";gsub(" ",",",$0);print $0}else{gsub("\\.",",",$1);FS=" ";print $0}}' tmp3 > tmp4
		awk -F"," '{if(NR==1){print $0;next}for(i=1;i<7;i++){printf "%s,",$i};printf "%s,%s,",$8,$8; for(i=9;i<NF;i++){printf "%s,",$i}; printf "%s\n",$NF}' tmp4 > kraken_table.csv

		rm tmp*
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' kraken_table.csv > tmp && mv tmp kraken_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .kraken.dat kraken_table.csv
		rm parsed* makeCSV.R
		sed "s/\.\.\./ /g" kraken_table.csv > tmp
		sed "s/,parsed/Kingdom.Phylum.Class.Order.Family.Genus.Species.Name,parsed/" tmp > tmp2
		sed "s/\"//g" tmp2 > tmp3
		awk 'BEGIN{FS=","}{if(NR==1){gsub("\\.",",",$1);FS=" ";gsub(" ",",",$0);print $0}else{gsub("\\.",",",$1);FS=" ";print $0}}' tmp3 > tmp4
		awk -F"," '{if(NR==1){print $0;next}for(i=1;i<7;i++){printf "%s,",$i};printf "%s,%s,",$8,$8; for(i=9;i<NF;i++){printf "%s,",$i}; printf "%s\n",$NF}' tmp4 > kraken_table.csv

		rm tmp*
 
	fi
}
function taxatorFunction {
	
	###################################################################################################################################################################################
	#this function take the taxator results (tax file), and parse it to leave only the tax id and number of mapped reads (% mapped reads if you specify an abundance in config file)
	###################################################################################################################################################################################

	for taxfile in $(ls -1 taxator*.tax)
	do
		grep -v "#" $taxfile |grep -v "@" |awk '{if($1!="")print $2}' |sort |uniq -c |sort -nr |awk '{print $2, $1/2}' > taxid.dat
		tsvfile=$(echo "$tsvfile" |sed "s/,/./g")
		mv taxid.dat parsed_$taxfile.dat
		echo "$taxfile file formated"
	done	
		total=$(ls -1 parsed_*.tax.dat |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "need at least 2 files to make a table"
		cp parsed_$taxfile.dat parsed__$taxfile.dat
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tax.dat taxator_table.csv
		rm parsed* makeCSV.R
		sed "s/ti\.//g" taxator_table.csv > tmp
		rm -f taxator_table && mv tmp taxator_table.csv

		TakeLineageFunction taxator_table.csv
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' taxator_table.csv > tmp && mv tmp taxator_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . tax.dat taxator_table.csv
		rm parsed* makeCSV.R
		sed "s/ti\.//g" taxator_table.csv > tmp
		rm -f taxator_table && mv tmp taxator_table.csv

		TakeLineageFunction taxator_table.csv
	fi
}
function centrifugeFunction {

	for tsvfile in $(ls -1 centrifuge*.tsv)
	do
		awk -F"\t" '{if(NR>1 && $3=="species")print $1"_"$5/2}' $tsvfile > centrifugeid.dat
		mv centrifugeid.dat parsed_$tsvfile.cf
		echo "$tsvfile file formated"
	done	
	total=$(ls -1 parsed*.tsv.cf |wc -l)
	if [ $((total)) -le 1 ]; then
		#echo "* need at least 2 files to make a table."
		cp parsed_$tsvfile.cf parsed__$tsvfile.cf
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .tsv.cf centrifuge_table.csv
		rm parsed*.cf makeCSV.R

		#take lineage
		echo 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Name' > tmp.csv 
		awk '{if(NR>1)print}' centrifuge_table.csv |while read line
		do
			genus=$(echo "$line" |awk -F"," '{print $1}' |awk -F"." '{print $1}')
			species=$(echo "$line" |awk -F"," '{print $1}' |awk -F"." '{print $2}')

			if [ "$species" == "" ];then
				echo "* no species was detected, aborting, check centrifuge table, genus: $genus, species: $species"
				exit
			fi

			nofetch=""
			while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]] || [[ "$nofetch" =~ "Bad Gateway!" ]]
			do
				if curl -s "http://www.ebi.ac.uk/ena/data/view/Taxon:$genus%20$species&display=xml" |awk 'BEGIN{band=0}{if($1=="<lineage>"){band=1;next}if(band==1){print}if($1=="</lineage>"){exit}}' > tmplin ;then
					touch tmplin
					nofetch=$(cat tmplin)
				else
					echo "curl error fetch, internet connection?, retrying after 10 seconds"
					nofetch=""
					wait 10
				fi
			done

			family=$(grep "rank=\"family\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$family" == "" ];then 
				family="unknow"
			fi
			order=$(grep "rank=\"order\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$order" == "" ];then 
				order="unknow"
			fi
			class=$(grep "rank=\"class\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$class" == "" ];then 
				class="unknow"
			fi
			phylum=$(grep "rank=\"phylum\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$phylum" == "" ];then 
				phylum="unknow"
			fi
			superk=$(grep "rank=\"superkingdom\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$superk" == "" ];then 
				superk="unknow"
			fi
			echo "$superk,$phylum,$class,$order,$family,$genus,$species,$genus $species"  >> tmp.csv
		done
		awk -F"," '{$1="";gsub(" ",",");print $0}' centrifuge_table.csv > tmpvalues.csv
		paste -d '\0' tmp.csv tmpvalues.csv > centrifuge_table.csv
		rm -f tmplin tmp.csv tmpvalues.csv
		awk -F"," '{$(NF-1)=$NF;$NF="";for(i=1;i<(NF-1);i++){printf "%s,",$i}printf "%s\n",$(NF-1)}' centrifuge_table.csv > tmp && mv tmp centrifuge_table.csv

	else
		#we call makeCSV.R to merge the results in a single file that contain the "raw data" for several analysis
		#parameters: work_directory pattern_file name_out_table
		makeCSV > makeCSV.R
		Rscript makeCSV.R . .tsv.cf centrifuge_table.csv
		rm parsed*.cf makeCSV.R

		#take lineage
		echo 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Name' > tmp.csv 
		awk '{if(NR>1)print}' centrifuge_table.csv |while read line
		do
			genus=$(echo "$line" |awk -F"," '{print $1}' |awk -F"." '{print $1}')
			species=$(echo "$line" |awk -F"," '{print $1}' |awk -F"." '{print $2}')

			if [ "$species" == "" ];then
				echo "* no species was detected, aborting, check centrifuge table, genus: $genus, species: $species"
				exit
			fi

			nofetch=""
			while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]] || [[ "$nofetch" =~ "Bad Gateway!" ]]
			do
				if curl -s "http://www.ebi.ac.uk/ena/data/view/Taxon:$genus%20$species&display=xml" |awk 'BEGIN{band=0}{if($1=="<lineage>"){band=1;next}if(band==1){print}if($1=="</lineage>"){exit}}' > tmplin ;then
					touch tmplin
					nofetch=$(cat tmplin)
				else
					echo "curl error fetch, internet connection?, retrying after 10 seconds"
					nofetch=""
					wait 10
				fi
			done

			family=$(grep "rank=\"family\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$family" == "" ];then 
				family="unknow"
			fi
			order=$(grep "rank=\"order\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$order" == "" ];then 
				order="unknow"
			fi
			class=$(grep "rank=\"class\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$class" == "" ];then 
				class="unknow"
			fi
			phylum=$(grep "rank=\"phylum\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$phylum" == "" ];then 
				phylum="unknow"
			fi
			superk=$(grep "rank=\"superkingdom\"" tmplin |awk '{print $2}' |awk -F"=" '{gsub("\"","");print $2}')
			if [ "$superk" == "" ];then 
				superk="unknow"
			fi
			echo "$superk,$phylum,$class,$order,$family,$genus,$species,$genus $species"  >> tmp.csv
		done
		awk -F"," '{$1="";gsub(" ",",");print $0}' centrifuge_table.csv > tmpvalues.csv
		paste -d '\0' tmp.csv tmpvalues.csv > centrifuge_table.csv
		rm -f tmplin tmp.csv tmpvalues.csv

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
	  patient_names[i]<-substr(list_of_files[i],1,nchar(list_of_files[i])-(nchar(pattr)))
	}
	
	# read in each table
	
	#read_counts <- lapply(list_of_files, read.table, sep="\t", header = FALSE, skip =2)
	if(pattr == ".dat.dat" || pattr == ".profiles.dat" || pattr == ".kraken.dat" || pattr == ".tsv.cf"){
		read_counts <- lapply(list_of_files, read.table, sep="_", header = FALSE)
	}else{
		read_counts <- lapply(list_of_files, read.table, sep=" ", header = FALSE)
	#read_counts <- lapply(read_counts, function(x) x[, c(1,2)])
	#read_counts <- lapply(read_counts, function(x) x[complete.cases(x),])
	}
	
	
	# for each table make the first col name OTU and the second the patient name
	
	if(pattr == ".dat.dat" || pattr == ".profiles.dat" || pattr == ".kraken.dat" || pattr == ".tsv.cf"){
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
	
	if(pattr == ".dat.dat" || pattr == ".profiles.dat" || pattr == ".kraken.dat" || pattr == ".tsv.cf"){
		for( i in 1:length(list_of_files)){
			name<-paste(as.character(read_counts[[i]][,1]))
		  	otu[i]<- list(name)
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
	write.csv(otu_table_noZeroes,outtable,quote = F)'
}


################################################################################################################
#begin the code
if [[ "$METHOD" =~ "METAPHLAN" ]] || [[ "$METHOD" =~ "CONSTRAINS" ]] || [[ "$METHOD" =~ "SIGMA" ]]; then
	if [ "$DEEPS" == "" ]; then
		echo "Metaphlan2, Constrains and Sigma needs depth of sequencing as parameter, put DEPTH=[SOME INTEGER] in the config file (see help)"
		exit
	fi
fi

if [ $((statusband)) -ge 2 ]; then
	cd ${RUTAINICIAL}
	for g in $METHOD
	do
		case $g in
			"PATHOSCOPE")
				echo "Parsing pathoscope files"
				pathoscopeFunction
			;;
			"METAPHLAN")
				echo "Parsing metaphlan files"
				metaphlanFunction
	   		;;
	   		"METAMIX")
				echo "Parsing metamix files"
	   			metamixFunction
	   		;;
	   		"SIGMA")
				echo "Parsing sigma files"
				sigmaFunction
	   		;;
	   		"CONSTRAINS")
				echo "Parsing constrains files"
	   			constrainsFunction
	   		;;
	   		"KRAKEN")
				echo "Parsing kraken files"
				krakenFunction
			;;
	   		"TAXATOR")
				echo "Parsing taxator files"
				taxatorFunction
			;;
			"CENTRIFUGE")
				echo "Parsing centrifuge files"
				centrifugeFunction
			;;
	   		*)
	   			echo "no method aviable for $METHOD"
	   		;;
		esac
	done

else
	echo "Invalid or Missing Parameters, print --help to see the options"
	echo "Usage: bash parseMethods.bash --workpath [files directory] --cfile [config file]"
	exit
fi
