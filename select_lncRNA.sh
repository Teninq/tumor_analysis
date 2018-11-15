lnc_file=$1
while read line
do
	ens_file=`grep $line gencode.v19.long_noncoding_RNAs.gtf`
	OLD_IFS="$IFS"
	IFS=$'\n'
	ens_gtf=($ens_file)
	all_exon=""
	for item in ${ens_gtf[@]}
	do
		type=`echo $item | awk '{print $3}'`
		if [[ "$type" == "transcript" ]] && [[ "$all_exon" != "" ]];then
			trans_id=`echo $item | awk '{print substr($12,2,length($12)-3)}'`
			all_exon=`echo $all_exon | tr "atcg" "ATCG"`
			if [[ "$strand" == "+" ]]; then
				echo $all_exon
			else
				all_exon=`echo $all_exon | tr "ATCG" "TAGC"`
				echo $all_exon
			fi
			echo ">"$line"_"$trans_id
			all_exon=""
		elif [[ "$type" == "transcript" ]] && [[ "$all_exon" == "" ]]; then
			trans_id=`echo $item | awk '{print substr($12,2,length($12)-3)}'`
			echo ">"$line"_"$trans_id
		fi

		strand=`echo $item | awk '{print $7}'`
		if [[ "$type" == "exon" ]];then
			chrid=`echo $item | awk '{print $1}'`
			exon_start=`echo $item | awk '{print $4}'`
			exon_end=`echo $item | awk '{print $5}'`
			exon_location=$chrid":"$exon_start"-"$exon_end
			# exon_len=`$exon_end - $exon_start`
			# exon_len=${exon_len#-}
			exon_seq=`samtools faidx /media/_EXTend2016/Lvchenkai2017/reference/ucsc.hg19.fasta $exon_location | awk '!/^>/{a=a""$0}END{print a}'`
			all_exon=${all_exon}${exon_seq}
		fi
	done
	all_exon=`echo $all_exon | tr "atcg" "ATCG"`
	if [[ "$strand" == "+" ]]; then
		echo $all_exon
	else
		all_exon=`echo $all_exon | tr "ATCG" "TAGC"`
		echo $all_exon
	fi
done < $lnc_file
