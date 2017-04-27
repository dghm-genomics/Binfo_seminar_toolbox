#!/bin/bash
## vs-varcall.sh REFERENCE BAM(s) / DSMZ / 170413
## variant calling, using VarScan, and generating consensus fasta
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
## It is expected to call "VarScan" (version 2) by "VarScan.jar" command! Otherwise you may adapt the input of the toolcheck function.
## It is expected to call samtools (version 0.1.19) by "samtools" command! Otherwise you may adapt the input of the toolcheck function.

VERSION=V1.1

DATE=$(date +"%Y-%m-%d-%H-%M-%S")
CID=$(date "+%N")
readonly DATE
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

generate-vcf(){
  echo -e "$4" >> "$2/$3"
  $vsvc_samtools faidx "$1"
  $vsvc_samtools mpileup -B -q 30 -f "$1" "$4"| $vsvc_varscan mpileup2cns --output-vcf 1 --min-coverage $vsvc_mincoverage --min-reads2 $vsvc_minreads2 --min-avg-qual $vsvc_minavgqual --min-var-freq $vsvc_minvarfreq --min-freq-for-hom $vsvc_minfreqforhom --p-value $vsvc_pvalue --strand-filter $vsvc_strandfilter > "$2/${4%.*}_vs_$DATE.vcf" 2>> "$2/$3"
  rm -f "${ref}.fai"
  echo "$2/${4%.*}_vs_${DATE}.vcf"
}

vcf-filter(){
  grep -v "^#" "$1"|
   tr ',' '.'|
   awk -v OFS="\t" '
   BEGIN{i = 1}
   {
      split($10,PVAL,":")
      split($4,c4,"")
      split($5,c5,"")
      if ($2$1 in SNP) {
        if (PVALMAX[$2$1] < PVAL[8]) {
          ID[$2$1] = $1
          REF[$2$1] = c4[1]
          SNP[$2$1] = c5[1]
          PVALMAX[$2$1] = PVAL[8]
        }
      }
      else {
        index_c[i] = $2$1
        i++
        ID[$2$1] = $1
        pos[$2$1] = $2
        REF[$2$1] = c4[1]
        SNP[$2$1] = c5[1]
         PVALMAX[$2$1] = PVAL[8]
      }
   }
   END {
      for (j = 1; j <= i-1; j++)
	print pos[index_c[j]],ID[index_c[j]],REF[index_c[j]],SNP[index_c[j]]
   }
   ' > "${1%.*}_${DATE}_filtered_vcf.txt"
  echo "${1%.*}_${DATE}_filtered_vcf.txt"
}

generate-fasta(){
  echo ">${4%.*}" > "$2/${4%.*}_cns.fasta"
  awk -v ipf="$1" '
    BEGIN {
      while (getline ipl < ipf) {
            if (ipl ~ /^>/) {
                i++
                split(ipl,seed)
                sub(/^>/,"",seed[1])
                header[i] = seed[1]
                pos = 1
             }else{
                for (j=pos; j<=(pos+length(ipl)-1); j++) {
                        bases_total++
                        index_c[bases_total] = j""header[i]
                        seq[index_c[bases_total]] = "N"
                }
                pos = j
              }
      }
    }
    {
      if ($4 == "."){
              seq[$1$2] = $3
      }
      else {
            seq[$1$2] = $4
      }
    }
    END {
      for(k=1; k<=bases_total; k++){
            print seq[index_c[k]]
      }
    }
  ' "$3"|
    tr -d '\n'|
    sed '$a\'|
    fold -w60 >> "$2/${4%.*}_cns.fasta"
}

### create help panel
#
if [ $# == 0 ]
then
    echo -e "\n\033[1mProgram\033[0m
      \tvs-varcall.sh (${VERSION}) - generate consensus fasta using VarScan2 [1] variant calling
\033[1mSYNOPSIS\033[0m
      \tvs-varcall.sh [OPTION]... [-r=REF INPUT-BAM(s)]...
\033[1mDESCRIPTION\033[0m

  Realise variant calling using VarScan2, and generate consensus fasta.

  \033[1m-c=MINCOVERAGE\033[0m
    VarScan parameter \"mincoverage\" [10]

  \033[1m-d=TARGET_DIR\033[0m
    specify target directory
      
  \033[1m-D\033[0m
    delete vcf file after run
      
  \033[1m-f=MINVARFREQ\033[0m
    VarScan parameter \"minvarfreq\" [0.8]
      
  \033[1m-h=MINFREQFORHOM\033[0m
    VarScan parameter \"minfreqforhom\" [0.75]
      
  \033[1m-p=PVALUE\033[0m
    VarScan parameter \"pvalue\" [0.01]
      
  \033[1m-q=MINAVGQUAL\033[0m
    VarScan parameter \"minavgqual\" [20]
      
  \033[1m-r=REF\033[0m
    REFERENCE-SEQUENCE-FASTA
      
  \033[1m-R=MINREADS\033[0m
    VarScan parameter \"minreads2\" [6]
      
  \033[1m-s=STRANDFILTER\033[0m
    VarScan parameter \"strandfilter\" [1]
	  
  \033[1m[1]\033[0m	VarScan 2: Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin,
	L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2:
	Somatic mutation and copy number alteration discovery in cancer by exome
	sequencing Genome Research DOI: 10.1101/gr.129684.111
	    
  \033[1mURL\033[0m	http://varscan.sourceforge.net
  "
  exit
fi


### read user options
#
declare -i vsvc_mincoverage=10 vsvc_minreads2=6 vsvc_minavgqual=20 vsvc_strandfilter=1
vsvc_CLEAR=0
target_dir=varcall-consensus-vs_run_$DATE 
vsvc_minvarfreq=0.8 
vsvc_minfreqforhom=0.75
vsvc_pvalue=0.01
vsvc_samtools=
vsvc_varscan=
ref=
while getopts c:d:Df:h:P:q:r:R:s: opt # 2>/dev/null
do
   case $opt in
      c) vsvc_mincoverage=$((OPTARG + 0));;
      d) target_dir=$OPTARG;;
      D) vsvc_CLEAR=1;;
      f) vsvc_minvarfreq=$OPTARG;;
      h) vsvc_minfreqforhom=$OPTARG;;
      P) vsvc_pvalue=$OPTARG;;
      q) vsvc_minavgqual=$((OPTARG + 0));;
      r) ref=$OPTARG;;
      R) vsvc_minreads2=$((OPTARG + 0));;
      s) vsvc_strandfilter=$((OPTARG + 0));;
      *) echo "\n [call-vcvs.sh] ERROR invalid user option!\n"
   esac
done
shift $((OPTIND-1))

if [ -z "$ref" ]
then
    echo -e "\n [vs-varcall] ERROR: Reference sequence was not specified!\n"
    sleep 2s
    exit 1
fi

### check for dependencies and provide target dir
#
vsvc_samtools=$(toolcheck samtools)
$vsvc_samtools 2> .toolcheck${CID}~ || true
SAMTOOLS_VERSION=$(grep "Version" .toolcheck${CID}~)
rm -f .toolcheck${CID}~
vsvc_varscan=$(toolcheck varscan)
[ ! -d "$target_dir" ] && { mkdir "$target_dir"; }





### main section
#
echo -e "\ncalling samtools:\n\n\t $vsvc_samtools\n\t $SAMTOOLS_VERSION\n" >> "$target_dir/varcall-vs_${ref%.*}_${DATE}.log"

for inpf
do
  ### call consensus with VarScan, generating vcf files
  #
  vcf_file=$(generate-vcf "$ref" "$target_dir" "varcall-vs_${ref%.*}_${DATE}.log" "$inpf")
  ### generate reduced table from vcf file including consensus
  ### bases and snps only
  #
  filtered_vcf=$(vcf-filter "$vcf_file")
  ### finally build consensus fasta
  #
  generate-fasta "$ref" "$target_dir" "$filtered_vcf" "$inpf"
  rm -f "$target_dir/$filtered_vcf" "$filtered_vcf"
  if [ $vsvc_CLEAR -eq 1 ]
  then
      rm -f "$vcf_file"
  fi
done
