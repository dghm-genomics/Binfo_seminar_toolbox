#!/bin/bash
## bwamapping.sh [OPTIONS] [-r=REFERENCE-FASTA FASTQ-FILES]
## ! READ1-FASTQ READ2-FASTQ (paired-end reads)
## arrange bwa-mem mapping / DSMZ / 170421
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
##
##
## It is expected to call "bwa" (versions 0.7.5 >= X <= 0.7.12) by "bwa" command! Otherwise you may adapt the input of the toolcheck function.
## It is expected to call samtools (version 0.1.19) by "samtools" command! Otherwise you may adapt the input of the toolcheck function.

VERSION=V1.1
DATE=$(date +"%Y-%m-%d-%H-%M-%S")
CID=$(date "+%N")
set -eu

toolcheck(){
      tool=$(if [ -x "$1" ]; then echo "$PWD/$1"; fi)
      tool=$(if [ -z "$tool" ]; then which "$1"; fi)
      if [ -x "$tool" ]
      then
          echo "$tool"
      else
          echo -e "\n [vs-varcall] ERROR: could not locate $1\n"
          sleep 2s
          exit
      fi
 }

### create help panel
#
if [ $# == 0 ]
then
    echo -e "\n\033[1mProgram\033[0m
      \tbwamapping.sh (${VERSION}) - perform mapping with bwa-mem [1]
\033[1mSYNOPSIS\033[0m
      \tbwamapping.sh [OPTIONS] [-r=RERENCE-DIRECTORY INPUT-FASTQ(s)]
\033[1mDESCRIPTION\033[0m

    \t This script does not support interleaved read files!

    \t\033[1m-d\033[0m\n\t  specify target dir
    
    \t\033[1m-r\033[0m\n\t  reference (fasta) or reference(s)-DIRECTORY
        
    \t\033[1m-t\033[0m\n\t  number of threads [1]
    
    
    \t\033[1m[1]\033[0m\t BWA-MEM: H. Li, ‘Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM’, ArXiv13033997 Q-Bio, Mar. 2013.
    "
    exit
fi

### read user options
#
target_dir= 
bm_ref=
bm_threads=1
while getopts :d:r:t: opt
do
   case $opt in
      d) target_dir="$OPTARG";;
      r) bm_ref="$OPTARG";;
      t) bm_threads="$OPTARG";;
      *) echo "\n [bwamapping] ERROR invalid user option!\n"; exit 1;;
   esac
done
shift $((OPTIND-1))

if [ -z "$bm_ref" ]
then
    echo -e "\n [bwamapping] ERROR: Reference sequence was not specified!\n"
    sleep 2s
    exit
fi

### check for dependencies and provide target dir
#
bm_samtools=$(toolcheck samtools)
bm_bwa=$(toolcheck bwa)
$bm_samtools 2> .toolcheck${CID}~ || true
SAMTOOLS_VERSION=$(grep "Version" .toolcheck${CID}~)
rm -f .toolcheck${CID}~
$bm_bwa 2> .toolcheck${CID}~ || true
BWA_VERSION=$(grep "Version" .toolcheck${CID}~)
rm -f .toolcheck${CID}~





### main section
### providing reference
#
ref_files=$(ls "$bm_ref"| grep "\.[Ff][Aa][Ss][Tt][Aa]$\|\.[Ff][AaSsRrNn][AaNnRrSs]$")
if [ -d "$bm_ref" ]
then
    bm_ref=$(echo "$bm_ref"| sed 's%/$%%g')
else
    bm_ref=.
fi

for ref in $ref_files
do
  if [ -z "$target_dir" ]
  then
	target_dir="${ref%.*}_mapping_results_$DATE"
	mkdir "${ref%.*}_mapping_results_$DATE"
  fi
  RDIR="REF$CID"
  mkdir "$RDIR"
  cp "$bm_ref/$ref" "$RDIR" || exit 1
  locref="$RDIR/${ref##*/}"
  input=$(echo "$@"|
    tr ' ' '\n'|
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
	      print read[i]";"read2[i]
	    }
	  }
	}'
  )
  for input_reads in $input
  do
    r1=$(echo "$input_reads"| cut -d ';' -f1)
    r2=$(echo "$input_reads"| cut -d ';' -f2)

    ### generate target id from read file 1 and reference id
    #
    targetid="$(echo "${r1%%.*}"|
    tr '_' '\n'|
     awk '$1 !~ /^[LlRrSs][[:digit:]][[:digit:]]{,1}[[:digit:]]{,1}$/ && $1 !~ /^00[[:digit:]]$/ && $1 !~ /^$/ {if(NR == 1) {i=i$1} else {i=i"-"$1}} END {print i}'
    )$(
        echo "_REF_${ref#*/}"|
         sed 's/\..\+$//g'|
          cut -c-40
     )"

    ### start mapping
    #
    logfile="$target_dir/$targetid.log"
    time(
      ref_len=$(grep -v ">" "$locref"| tr -d '\n'| wc -c)
      ref_header=$(grep ">" "$locref")
      echo -e "BWA mapping V1.0 logfile [$(date +"%Y-%m-%d-%H-%M-%S")]\n" > "$logfile"
      echo -e "### Reference" >> "$logfile"
      echo -e "#" >> "$logfile"
      echo -e "filename:\t${locref##*/}" >> "$logfile"
      echo -e "dir:\t\t${locref%%/*}" >> "$logfile"
      echo -e "fasta-header:\t${ref_header}" >> "$logfile"
      echo -e "length:\t\t${ref_len} bp" >> "$logfile"
      echo -e "size:" >> "$logfile"
      du -h "$r1" "$r2"| sed 's/^/\t/g' >> "$logfile"
      echo -e "\n### Mapping" >> "$logfile"
      echo -e "#" >> "$logfile"
      echo -e "bwa version: $BWA_VERSION" >> "$logfile"
      echo -e "samtools version: $SAMTOOLS_VERSION" >> "$logfile"
      if [ ! -f "$locref.amb" -o ! -f "$locref.ann" -o ! -f "$locref.bwt" -o ! -f "$locref.pac" -a ! -f "$locref.sa" ]
      then
	  $bm_bwa index -a is "${locref}"
      fi
      if [ ! -f "$locref.fai" ]
      then
	  $bm_samtools faidx "${locref}"
      fi
      $bm_bwa mem -t $bm_threads "$locref" "$r1" "$r2" > "$targetid.sam"
      $bm_samtools view -bt "$locref.fai" "$targetid.sam" -o "${targetid}.bam"
      $bm_samtools sort "${targetid}.bam" > "${targetid}_sorted.bam"
    ) 2>> "$logfile"

    ### bam-sort and statistics
    #
    [ -f "${targetid}_sorted.bam" ] && { mv "${targetid}_sorted.bam" "$target_dir/${targetid}.bam"; }
    echo -e "\n### Bam characteristics" >> "$logfile"
    echo -e "#" >> "$logfile"
    $bm_samtools flagstat "${targetid}.bam" >> "$logfile"

    ### environmental protection
    #
    rm -f "$targetid.sam" "${targetid}_sorted.bam" "$targetid.bam"

    echo -e "\n [bwamapping.sh] done with \n
      \t\t$ref
      \t\t$r1
      \t\t$r2\n"
  done
  rm -f "$locref.amb" "$locref.ann" "$locref.bwt" "$locref.pac" "$locref.sa" "$locref.fai" "$locref"
  rmdir "$RDIR"
done
