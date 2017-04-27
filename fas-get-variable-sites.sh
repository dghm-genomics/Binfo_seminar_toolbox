#!/bin/bash
## fas-get-variable-sites.sh VCF-FILE / DSMZ /170426
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

VERSION=V1.1

CID=$(date "+%N")
readonly CID
set -eu

linebreaker(){
  awk 'NR == 1 {
        fheader = $0
     }
     $1 !~ ">" {
       fseq = fseq $1
     }
     $1 ~ ">" && NR != 1 {
        print fheader
        print fseq
        fseq = ""
        fheader = $0
      }
      END{
        print fheader
        print fseq
      }
    ' $1
}

printcolumns(){
  mkdir .fragments_$CID
  cd .fragments_$CID
  cp ../$2 ./
  i=1000
  awk '{
    printf $0","
    frac = int(NR/15000)
    if (NR - frac * 15000 == 0){print ""}}' "$2"|
    sed 's/^,\|,$//g'| sed '$a\'|
    while read inpl
    do
        i=$((i+1))
        cut -c"$inpl" "../$1" > p_$i
    done
    rm "$2"
    paste -d '' $(ls| tr '\n' ' ')
    cd ..
    rm -r .fragments_$CID
}

### create help panel
#
if [ $# == 0 ]
then
  echo -e "
  \033[1mProgram\033[0m
      fas-get-variable-sites.sh (${VERSION}) - extract variable sites from an alignment
      
  \033[1mSYNOPSIS\033[0m
      fas-get-variable-sites.sh FASTA-FILE
      
  \033[1mDESCRIPTION\033[0m
  
    \033[1m-n\033[0m
      filter \"N symbols\"
  "
  exit
fi

### read user options
#
while getopts n opt # 2>/dev/null
do
   case $opt in
       n) fasgtvs_no_n=1;;
       :) echo echo "\n [fas-get-variable-sides] ERROR invalid user option!\n"
   esac
done
shift $((OPTIND-1))

: ${fasgtvs_no_n:=NDEF}





grep ">" "$1" > header_$CID
linebreaker "$1"| grep -v ">" > "${1}_sl_${CID}"

awk 'NR == 1{
  ref_c = $0
}{
  print ref_c
}' "${1}_sl_${CID}" >> "ref_sl_${CID}"

glen=$(awk -v FS="" 'END{print NF}' "ref_sl_${CID}")

if [ "$fasgtvs_no_n" == NDEF ]
then
    cmp -l "${1}_sl_${CID}" "ref_sl_${CID}"|
    sed 's/\s\+/\t/g;s/^\t//g'|
    cut -f1|
    awk -v glen="$glen" '{
      a = int( $1 / glen )
      print $1 - ( a * glen ) - a
    }'| sort -ug > "vsites_${1%.*}_$CID.txt"
else
 cmp -bl "${1}_sl_${CID}" "ref_sl_${CID}"|
     sed 's/^\s\+//g;s/\s\+/\t/g'| 
     cut -f1,3,5|
      awk -v glen="$glen" '{
	a = int( $1 / glen )
	indx = $1 - ( a * glen ) - a
	b[indx]
	if ($2 == "N" || $3 == "N")
	n[indx]
      }END{
	for (i in b){
	  if (i in n){
	  }else{
	    print i
	  }
	}
      }
      '| 
      sort -ug > "vsites_${1%.*}_$CID.txt"
fi

printcolumns "${1}_sl_${CID}" "vsites_${1%.*}_$CID.txt" > alg_segm_$CID

paste header_$CID alg_segm_$CID| tr '\t' '\n' > "snp_algnm_${1%.*}.fasta"

rm -f "${1}_sl_${CID}" "ref_sl_${CID}" alg_segm_$CID header_$CID

echo -e "\n [fas-get-variable-sides] variable sites: $(wc -l "vsites_${1%.*}_$CID.txt")\n"
