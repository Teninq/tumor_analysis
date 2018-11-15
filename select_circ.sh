#!/bin/bash
len=$#
if [ $len != "1" ];then
	echo ""
	echo "Usage sh select_circ.sh test.txt"
	echo ""
	exit 0
fi
circ_file=$1
while read line
do
	res=`grep $line "all.output" | awk '{split($10,id,",");print id[1]"\t"$11}'`
	echo "=================================="
	echo $res
	echo $line
	chrid=`echo $line | awk '{split($0,a,":");print a[1]}'`
	ens_id=`echo $res | awk '{print $1}'`
	#echo $ens_id
	circ_start=`echo $line | awk '{split($0,a,":");split(a[2],b,"|");print b[1]}'`
	circ_end=`echo $line | awk '{split($0,a,":");split(a[2],b,"|");print b[2]}'`
	ens_file=`grep $ens_id /media/_EXTend2016/zhangli/HCC-RNAseq/index/gencode.v19.annotation.gtf`
	OLD_IFS="$IFS"
	IFS=$'\n'
	ens_gtf=($ens_file)
	#ens_gtf=`echo $ens_file | awk '{print $0}'`
	for item in ${ens_gtf[@]}
	do
		type=`echo $item | awk '{print $3}'`
		if [ $type == "transcript" ];then
			trans_start=`echo $item | awk '{print $4}'`
			trans_end=`echo $item | awk '{print $5}'`
			if [ $trans_start -le $circ_start ] && [ $trans_end -ge $circ_end ];then
				echo "------------------------------------"
				echo "transcript:"$chrid":"$trans_start"-"$trans_end
				flag="TRUE"
			else
				flag="FALSE"
				#echo "------------------------------------"
				continue
			fi
		fi
		if [ $type == "exon" ] && [ $flag == "TRUE" ];then
			exon_start=`echo $item | awk '{print $4}'`
			exon_end=`echo $item | awk '{print $5}'`
			dis_start=`expr $exon_start - $circ_start`
			dis_end=`expr $exon_end - $circ_end`
			dis_start=${dis_start#-}
			dis_end=${dis_end#-}
			#echo "Dis:"$dis_start $dis_end
			if [ $dis_start -le 2 ] || [ $dis_end -le 2 ];then
				echo "Exon:"$chrid":"$exon_start"-"$exon_end
			fi
		fi
	done
	#IFS="$OLD_IFS"
    #echo $ens_file | awk '{print $3}'

#        if [$type == "transcript"];then
#            trans_start=`echo $line1 | awk '{print $4}'`
#            echo $trans_start
done < $circ_file














