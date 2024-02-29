#!/bin/bash


#----
# Initialise
#----

# Set environment variables

script=snp2gene.sh
script_dir=$(dirname $(readlink -f $( which $0) ) )
init_dir=`pwd`
shopt -s extglob

# Set default parameters

snp_list=`mktemp`

case `hostname` in
	muckleflugga)
	database=/opt/working/wilson/app_full_stuff/locuszoom_v1.4/data/database/locuszoom_hg19.db
	sqlite3=/usr/local/share/anaconda/bin/sqlite3	
	;;

	*ecdf.ed.ac.uk)
	database=/exports/igmm/eddie/wilson-lab/apps_by_us_full_stuff/locuszoom_v1.4/data/database/locuszoom_hg19.db
	. /etc/profile.d/modules.sh
	which R &> /dev/null || module load R/3.3.2
	sqlite3=/usr/bin/sqlite3
	;;
esac

cytobase=${script_dir}/cytobands.db
window=250000
export=""
verbose=0
summarise=0

# Set error trap

set -eE
trap "rm ${snp_list}*" ERR 

# Define functions

print_help() {

	echo "==========="
	echo "$script"
	echo "==========="
	git --git-dir $script_dir/.git log -1 | awk 'BEGIN{print "Author: Paul Timmers <paul.timmers@ed.ac.uk>"}; $1~/commit/||$1~/Date/ {print}' | sed 's/commit/Version:/; s/:  /:/'
	echo " "


	case "$1" in
		*(-)f|*(-)file)
			echo " "
			echo "Showing help for --file"
			echo " "
    		echo "usage:"
    		echo "$script -f (filename) ..."
    		echo " "
    		echo "Find gene names for all SNPs listed in the file"
    		echo "These should be newline-separated. For example:"
    		echo " "
			echo "cat 'snp_list.txt'"
			echo " > rs6857"
			echo " > rs7976168"
			echo " > 3:52742346_TGA_T"
    		echo " "
    		echo "$script -f snp_list.txt rs429358"
    		echo " "
			echo "Header optional."
    		echo " "
    		;;

		*(-)w|*(-)window)
			echo " "
			echo "Showing help for --window"
			echo " "
    		echo "usage:"
    		echo "$script -w (250000) ..."
    		echo " "
    		echo "Set the number of base pairs around the sentinel SNP that is checked for genes"
    		echo "You can also attach Kb and Mb (without spaces) for 1,000 and 100,000 bases, respectively"
    		echo "For example:"
    		echo ""
			echo "$script -w 0.5Mb rs3764814 3:52742346_TGA_T"
    		echo " "
    		;;

		*(-)s*|*(-)summar*)
			echo " "
			echo "Showing help for --summary"
			echo " "
    		echo "usage:"
    		echo "$script -s ..."
    		echo " "
    		echo "For each SNP,  summarise the results of the search into a single name and append to the output."
    		echo "If multiple gene names are available, a compound is created from the closest upstream gene name"
			echo "and the closest downstream gene name (e.g. CELSR2/PSRC1)"
			echo " "
			echo "NOTE: Pseudogenes and RNAs are not included in compound gene names to keep names short"
    		echo " "
    		;;

		*(-)b|*(-)build)
			echo " "
			echo "Showing help for --build"
			echo " "
    		echo "usage:"
    		echo "$script -b (hg19) ..."
    		echo " "
    		echo "Specify which genome build to check for SNP and gene positions. The default is 'hg19' but you can also"
    		echo "select 'hg38' or 'hg18'. All data is contained within the locuszoom database, so please update that if"
    		echo "you would like to update the genome build."
    		echo " "
    		;;

		*(-)e|*(-)export)
			echo " "
			echo "Showing help for --export"
			echo " "
    		echo "usage:"
    		echo "$script -e (filename) ..."
    		echo " "
    		echo "Save the results of the search into the specified file. For example:"
    		echo ""
			echo "$script -e ../p01_results/snp_info.tsv rs3764814 rs7976168"
    		echo " "
    		;;

		*(-)v|*(-)verbose)
			echo " "
			echo "Showing help for --verbose"
			echo " "
    		echo "usage:"
    		echo "$script -v ..."
    		echo " "
    		echo "Make the script more talkative. Output will contain script version, the number of SNPs,"
    		echo "as well as any other arguments entered. While searching for SNPs, % progress is reported"
    		echo "and the final output is displayed in an easy-to-read format"
    		echo " "
    		;;

		*)
    		echo "This  script  will  search  around  each  specified  SNP  for  nearby  genes."
    		echo "The top 3 closest genes are reported by name and distance in basepairs, where"
    		echo "negative  distance implies upstream of coding sequence and positive  distance"
    		echo "implies downstream of coding sequence."
    		echo " "
    		echo "usage:"
    		echo "$script [options] [snp1] [snp2] ..."
    		echo " "                   
    		echo "options:"
    		echo "-h   --help [option]           show brief help, or detailed help for specific option"
    		echo "-f   --file (filename)         read multiple SNP names from single file"
    		echo "-w   --window (250000)         number of base pairs flanking SNP checked for genes"
    		echo "-b   --build (hgXX)            select genome build to use for SNP and gene positions"
    		echo "-e   --export (filename)       export results to a file with the name of your choice"
    		echo "-s   --summary                 summarise results into a single name per SNP"
		echo "-s2  --summary2                like --summary but shortens overlapping gene names"
    		echo "-v   --verbose                 output more info, such as version, arguments, and progress"
    		echo " "
    		;;
	esac
}


find_gene() {
	snp=$1
	chr=$2
	pos=$3
	window=$4
	database=$5
	cytobase=$6
	sqlite3=$7
	min=`echo "$pos - $window" | bc`
	max=`echo "$pos + $window" | bc`
	build=`echo $database | grep -q "hg38" && echo hg38 || echo hg19`

	cyto=`${sqlite3} -separator "	" $cytobase "SELECT cyto FROM cyto_pos WHERE build = '$build' AND chr=$chr AND start <= $pos AND end >= $pos LIMIT 1"`

	chrs=({1..22} X Y M MT)
	${sqlite3} -separator "	" $database \
	"SELECT DISTINCT geneName,cdsStart,cdsEnd FROM refFlat WHERE chrom=\"chr${chrs[$chr-1]}\" AND cdsEnd > $min AND cdsStart < $max" \
	| awk -v OFS="\t" -v snp=$snp -v chr=$chr -v pos=$pos -v window=$window -v cyto=$cyto '
	function abs(v) {return v < 0 ? -v : v} 
	function dist(v,a,b) { if(a<=v && b>=v) {return 0} else {return abs(a-v) < abs(b-v) ? a-v : b-v}} 
	gene[$1] == 0 {print dist(pos,$2,$3),$1,dist($2,$3); gene[$1]++}' \
	| sort -k1Vd,1 -k3n,3 | cut -f1-2 \
	| awk -v snp=$snp -v chr=$chr -v pos=$pos -v window=$window -v cyto=$cyto \
	'{distance[NR]=$1; gene[NR]=$2} 
	END{printf "\"%s\"\t%i\t%i\t%i\t%i\t",snp,chr,pos,NR,window; for (i=1; i<=3; i++) { if(i <= NR){printf "\"%s\"\t%i\t",gene[i],distance[i]} else {printf "%s\t%s\t","NA","NA"}} {printf "\"%s\"\n",cyto}}'
	
	

}

export -f find_gene


#----
# Get options
#----

# Show help if no arguments specified

if [[ -z $1 ]]; then
    print_help 
    exit 1
fi


# Otherwise, parse through arguments

while test $# -gt 0; do
        case "$1" in
                +(-)h|+(-)help)
                        shift
                        print_help $1 
                        exit 0
                        ;;

                +(-)w|+(-)window)
						shift
						
						if `echo $1 | grep -iqE "[ac-jlno-z]"`; then
							echo "$script: invalid window size '$1'"
                            exit 1
						fi

						if [[ -z $1 ]]; then
							echo "$script: no window size specified"
							exit 1
						fi

						window=`echo $1 | sed 's/[^0-9\.]*//g'`
						echo $1 | grep -qiE "k[b]?" && window=$(bc -l <<< "$window * 1000") 
						echo $1 | grep -qiE "m[b]?" && window=$(bc -l <<< "${window} * 1000000")

						shift
						;;

                +(-)b|+(-)build)
						shift
						if [ $# -gt 0 ] && [ ${1:0:1} != "-" ]; then
							database=`echo $database | sed "s/hg19/$1/"`
							if [[ ! -f $database ]]; then
								echo -en "$script: invalid build specified.\n  Please choose from: "
								ls `dirname $database`/locuszoom*.db | sed 's=.*locuszoom_==g; s/.db//g' | awk -v ORS=" " '1'
								echo ""
								exit 1
							fi
						else
                            echo -en "$script: invalid build specified.\n  Please choose from: "
							ls `dirname $database`/locuszoom*.db | sed 's=.*locuszoom_==g; s/.db//g' | awk -v ORS=" " '1'
							echo ""
							exit 1
                        fi
                        shift
						;;

                +(-)e|+(-)export)
						shift
						if [ $# -gt 0 ] && [ ${1:0:1} != "-" ]; then
							export=$1
						else
                            echo "$script: no output filename specified"
                            exit 1
                        fi
                        shift
						;;

				+(-)s|+(-)summar*[^2])
						shift
						summarise=1
						;;
				+(-)s2|+(-)summar*2)
                                                shift
                                                summarise=2
                                                ;;


                +(-)v|+(-)verbose)
						shift
						verbose=1
						;;

                +(-)f|+(-)file)
                        shift
                        if [ $# -gt 0 ] && [ ${1:0:1} != "-" ]; then
                        	
                        	if [[ ! -f $1 ]]; then
                        		echo "$script: '$1': No such file or directory"
                        		exit 1
                        	fi

                        	cat=cat
                        	file `readlink -f $1` | grep -q "gzip" && cat=zcat 

                        	eval $( $cat $1 | head -1 | awk '{for(i = 1; i <= NF; i++) if(tolower($i)~/rs[^0123456789]|marker|snp|id$|name/) {print "rsid_col="i; break }}' )
                        	
                        	if [[ -z $rsid_col ]]; then
                        		$cat $1 | awk '{print $1}' >> $snp_list
                        	else
                        		$cat $1 | awk -v rs=$rsid_col 'NR>1 {print $rs}' >> $snp_list
                        	fi
                        	
                            shift
                        else
                            echo "$script: no file specified"
                            exit 1
                        fi
                        ;;

                *)
					if [ ${1:0:1} == "-" ]; then
							echo "$script: unrecognised argument '$1'"
                            exit 1
					fi

						if [[ -f $1 ]]; then
							cat=cat
                        	file `readlink -f $1` | grep -q "gzip" && cat=zcat 

                        	eval $( $cat $1 | head -1 | awk '{for(i = 1; i <= NF; i++) if(tolower($i)~/rs[^0123456789]|marker|snp|id$|name/) {print "rsid_col="i; break }}' )
                        	
                        	if [[ -z $rsid_col ]]; then
                        		$cat $1 | awk '{print $1}' >> $snp_list
                        	else
                        		$cat $1 | awk -v rs=$rsid_col 'NR>1 {print $rs}' >> $snp_list
                        	fi
                        else
						echo $1 >> $snp_list
						fi

						shift
                    	;;
        esac
done


#----
# Compile SNP list
#----



n_snps=$(timeout 5s bash <<EOT
wc -l < $snp_list
EOT
)

if [[ $n_snps -eq 0 ]]; then
	echo "$script: no SNPs specified"
	exit 1
fi


if [[ $n_snps -le 10000 ]]; then
	awk 'snp[$0] == 0 {snp[$0]++; print}' ${snp_list} > ${snp_list}.t && mv ${snp_list}.t ${snp_list}
else
	n_snps=`echo "${n_snps}"`
fi

awk '$1 ~ /rs/ {print > FILENAME".rsid"; next} {print > FILENAME".other"}' $snp_list


if [[ -f ${snp_list}.rsid ]]; then
	rsids=`awk '{printf "\"%s\",",$1}' ${snp_list}.rsid | sed 's/,$//'`
fi


if [[ -f ${snp_list}.other ]]; then
	sed 's/chr//I' ${snp_list}.other | awk -v list=${snp_list}.other -F":|_" '{getline f < list} {printf "%s\t%i\t%i\n",f,$1,$2}' > ${snp_list}.otherpos
fi




#----
# Print header
#----

if [[ $verbose -gt 0 ]]; then
	echo "==========="
	echo "snp2gene.sh"
	echo "==========="
	git --git-dir $script_dir/.git log -1 | awk 'BEGIN{print "Author: Paul Timmers <paul.timmers@ed.ac.uk>"}; $1~/commit/||$1~/Date/ {print}' | sed 's/commit/Version:/; s/:  /:/'
	echo " "
	echo " "
	echo "ARGUMENTS"
	(echo -e "n_snps:\t$n_snps"
	echo -e "window:\t$window"
	echo -e "database:\t`basename $database`"
	if [[ ! -z $export ]]; then
		echo -e "export:\t`readlink -f $export`"
	fi
	) | column -t
	echo " "
fi





#----
# START
#----

if [[ $verbose -gt 0 ]]; then
	echo -en "Finding chromosome and base pair positions... "
fi

if [[ $n_snps -le 10000 ]]; then
	${sqlite3} -separator "	" $database \
	"CREATE TEMP VIEW snp_pos_trans AS SELECT rs_orig as snp,chr,pos FROM snp_pos p INNER JOIN refsnp_trans t ON (t.rs_current = p.snp); 
	SELECT snp,chr,pos FROM snp_pos_trans WHERE snp IN ($rsids)" > ${snp_list}.snppos
else
(${sqlite3} -separator "	" $database <<EOF
.load '$script_dir/csv'
CREATE VIRTUAL TABLE temp.t1 USING csv(filename='$snp_list');
CREATE TEMP VIEW snp_pos_trans AS SELECT rs_orig as snp,chr,pos FROM snp_pos p INNER JOIN refsnp_trans t ON (t.rs_current = p.snp); 
SELECT snp,chr,pos FROM snp_pos_trans p INNER JOIN t1 t ON (p.snp = t.c0) 
EOF
) > ${snp_list}.snppos
fi


if [[ $verbose -gt 0 ]]; then
	echo -e "done."
	echo -e "Finding closest genes... "
fi

if [[ -f ${snp_list}.otherpos ]]; then
	awk '$2 > 0 && $3 > 0' ${snp_list}.otherpos >> ${snp_list}.snppos
fi


if [[ $verbose -gt 0 ]]; then
	echo -e "rsid\tchr\tpos\tn_genes\twindow\tgene1\tdist1\tgene2\tdist2\tgene3\tdist3\tcyto" > ${snp_list}.genepos
	cat ${snp_list}.snppos | parallel --bar --col-sep="\t" -j 20 --no-notice find_gene {1} {2} {3} $window $database $cytobase $sqlite3 >> ${snp_list}.genepos
	

	if [[ $summarise -gt 0 ]]; then
		echo -en "Summarising results... "
		Rscript ${script_dir}/summarise${summarise}.R ${snp_list}.genepos ${snp_list}
		echo -e "done.\n"
		sum="NA"
	fi

	(awk -v sum=$sum -v OFS="\t" 'NR==1 {print} ARGIND == 1 {gsub(/"/,"",$0); rsid[$1]=$0; next} rsid[$1] > 0 {print rsid[$1]; next} {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",sum}' ${snp_list}.genepos ${snp_list}) | tee ${snp_list}1.genepos | column -t
	mv ${snp_list}1.genepos ${snp_list}.genepos

else
	if [[ $summarise -gt 0 ]]; then
		echo -e "rsid\tchr\tpos\tn_genes\twindow\tgene1\tdist1\tgene2\tdist2\tgene3\tdist3\tcyto" > ${snp_list}.genepos
		cat ${snp_list}.snppos | parallel --col-sep="\t" -j 20 --no-notice find_gene {1} {2} {3} $window $database $cytobase $sqlite3 >> ${snp_list}.genepos
	
	
		Rscript ${script_dir}/summarise${summarise}.R ${snp_list}.genepos ${snp_list}
		awk -v OFS="\t" 'NR==1 {print} ARGIND == 1 {gsub(/"/,"",$0); rsid[$1]=$0; next} rsid[$1] > 0 {print rsid[$1]; next} {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"}' ${snp_list}.genepos ${snp_list} | tee ${snp_list}1.genepos
		mv ${snp_list}1.genepos ${snp_list}.genepos

	else

		echo -e "rsid\tchr\tpos\tn_genes\twindow\tgene1\tdist1\tgene2\tdist2\tgene3\tdist3\tcyto" | tee ${snp_list}.genepos
		cat ${snp_list}.snppos | parallel --col-sep="\t" -j 20 --no-notice find_gene {1} {2} {3} $window $database $cytobase $sqlite3 | sed 's/"//g' | tee -a ${snp_list}.genepos
		awk -v OFS="\t" 'ARGIND == 1 {rsid[$1]=$0; next} rsid[$1] == 0 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"}' ${snp_list}.genepos ${snp_list} | tee -a ${snp_list}.genepos
	fi

	
fi


if [[ ! -z $export ]]; then
	cat ${snp_list}.genepos > $export
fi

rm ${snp_list}*
