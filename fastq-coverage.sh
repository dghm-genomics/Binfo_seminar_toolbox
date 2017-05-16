#!/bin/bash
## fastq-coverage.sh -g [Option] [-g=GENOME_LENGTH FILES] - calculation of theoretical coverage values given an expected genome length / DSMZ / 170424
## Copyright (c) 2017 Leibniz-Institut DSMZ - Deutsche Sammlung von Mikroorganismen und Zellkulturen GmbH
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
## documentation files (the "Software"), to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
## and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or substantial portions
## of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
## TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
## CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.

VERSION=V1.2
DATE=$(date +"%Y-%m-%d-%H-%M-%S")
CID=$(date +"%N")
set -eu

open_file(){
   first_file=$(echo "$1"| head -1)
   ftype=${first_file#*.}
   case "$ftype" in
      fastq.xz) xzcat $1 ;;
      fsq.xz) xzcat $1 ;;
      fq.xz) xzcat $1 ;;
      fastq.gz)	zcat $1 ;;
      fsq.gz)	zcat $1 ;;
      fq.gz) 	zcat $1 ;;
      fastq)	cat $1 ;;
      fsq)	cat $1 ;;
      fq)	cat $1 ;;
   esac
}

### create help panel
#
if [ $# -eq 0 ]
then
  echo -e "
  \033[1mProgram\033[0m
      fastq-coverage.sh (${VERSION}) - analyse fastq files and infer coverage

  \033[1mSYNOPSIS\033[0m
      fastq-coverage.sh [OPTION]... [-g=GENOME_LENGTH FASTQ-FILE...]

  \033[1mDESCRIPTION\033[0m

    The script figures the amount of reads, bases, and bases with a Q-value >= 20
    from fastq files and calculates:

      coverage = bp / expected genome length

      weighted coverage = weighted bp / expected genome length.

    All values will be stored in a fastq-coverage table file.

    It is possible to filter fastq files by their coverages (-s/-d).

    \033[1m-d=MIN_COVERAGE\033[0m
      specify coverage to deselect fastq file(s) >= MIN_COVERAGE

    \033[1m-g=GENOME_LENGTH\033[0m
      specify expected genome length

    \033[1m-l\033[0m
      just list the file ids but don't copy selection/deselection

    \033[1m-s=MIN_COVERAGE\033[0m
      specify coverage to select fastq file(s) >= MIN_COVERAGE

    \033[1m-t=STAT_TABLE_FILE\033[0m
      specify name of already existing fastq-coverage table file
  "
  exit
fi

### read user options
#
declare -i fc_cov_select=-1 fc_cov_deselect=-1
fc_gl=
fc_list=
fc_stat_table=
fc_selected_ids=
while getopts d:g:ls:t: opt # 2>/dev/null
do
   case $opt in
      d) fc_cov_deselect=$((OPTARG + 0));;
      g) fc_gl=$((OPTARG + 0));;
      l) fc_list=ON;;
      s) fc_cov_select=$((OPTARG + 0));;
      t) fc_stat_table=$OPTARG;;
      *) echo "\n [fastq-coverage] ERROR invalid user option!\n"; exit 1;;
   esac
done
shift $((OPTIND-1))

if [ -z "$fc_gl" ] && [ -z "$fc_stat_table" ]
then
    echo -e "\n [fastq-coverage] ERROR: please specify genome length (-g).\n"
    exit 1
fi





if [ -z "$fc_stat_table" ]
then
    echo -e "file-id\tfile size\tread count\tbases total\tproper weighted bases\taverage read length\ttheoretical coverage\tweighted coverage"| tee "fastq-coverage_${DATE}.txt"
    for fq_file
    do
      file_info=$(du -h "$fq_file"| awk '{print $2"\t"$1}')
      base_count=$(open_file "$fq_file"| sed -n '4~4p'| tr -d '\n'| wc -c)
      w_base_count=$(open_file "$fq_file"| sed -n '4~4p'| tr -d '\n!"#$%&'"'"'()*+,-./01234'| wc -c)
      read_count=$(open_file "$fq_file"| wc -l| awk '{print $1 / 4}')
      calculations=$(echo "$base_count $w_base_count $read_count $fc_gl"| awk '{printf "%4.0f\t%4.1f\t%4.1f\n", ($1 / $3), ($1 / $4), ($2 / $4)}')
      echo -e "$file_info\t$read_count\t$base_count\t$w_base_count\t$calculations"| tee -a "fastq-coverage_${DATE}.txt"
      fc_stat_table="fastq-coverage_${DATE}.txt"
    done
fi

if [ "$fc_cov_deselect" -ne -1  ]
then
    fc_selected_ids=$(awk -v min_cov="$fc_cov_deselect" 'NR != 1 && $8 < min_cov {print $1}' "$fc_stat_table")
elif [ "$fc_cov_select" -ne -1 ]
then
    fc_selected_ids=$(awk -v min_cov="$fc_cov_select" 'NR != 1 && $8 >= min_cov {print $1}' "$fc_stat_table")
fi
if [ "$fc_selected_ids"  ]
then
    file_count=$(echo -e "$fc_selected_ids"| wc -l)
    echo -e "$fc_selected_ids\n$file_count"
    if [ -z "$fc_list" ]
    then
	fc_copy_ids=$(
	  echo -e "$fc_selected_ids"|
	    awk '
	      {
		read_file=$1
		sub(/_+[Rr][12]/,"",$1)
		if (read_file ~ /[Ff][Aa][Ss][Tt][Qq]/ && read_file ~ /[Rr]0*1/){
		  read[$1] = read_file
		}
		if (read_file ~ /[Ff][Aa][Ss][Tt][Qq]/ && read_file ~ /[Rr]0*2/){
		  read2[$1] = read_file
		}
	      }
	      END{
		for(i in read2){
		  if (i in read){
		    printf read[i]" "read2[i]" "
		  }
		}
	      }'
	)
	mkdir "fastq_file_selection_${DATE}" || exit 1
	cp ${fc_copy_ids} "fastq_file_selection_${DATE}"
	echo -e "\n [fastq-coverage] selection\n\n$(echo $fc_copy_ids| tr ' ' '\n'| nl)\n"
	echo -e "\n copied to fastq_file_selection_${DATE}\n"
    fi
fi
