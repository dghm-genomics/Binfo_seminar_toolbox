#!/bin/bash
## fas-cut.sh [FASTA-FILE INTERVAL-FILE] / DSMZ / 170421
## cut specified positions from an alignment (fasta format)
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

DATE=$(date +"%Y-%m-%d-%H-%M-%S")
CID=$(date "+%N")
readonly CID
set -eu

### generate a "single-line per sequence" version
### of input (fasta file)
#
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
  ' "$1"
}

### remove overlaps in the given interval set
### and generate a continuous position vector
#
interval_expand(){
  awk '
    {
      if ($1 > $2){
        sp = $2
        ep = $1
      }else{
        sp = $1
        ep = $2
      }
      for (i=sp; i<=ep; i++){
        c[i] = i
      }
    }END{
      for (j in c){
        print c[j]
      }
    }
  ' "$1"
}

### create help panel
#
if [ $# == 0 ]
then
  echo -e "\n\033[1mProgram\033[0m
    \tfas-cut.sh (${VERSION}) - cut columns from an alignment
\033[1mSYNOPSIS\033[0m
    \tfas-cut.sh [FASTA-FILE INTERVAL-TABLE]
\033[1mDESCRIPTION\033[0m

    \t Intervals must be given as tab-separated list (1  5)
    "
	    exit
fi

### main section
#
interval_expand "$2" > cut_pos_$CID
linebreaker "$1" > "${1}_sl_${CID}"

awk -v cut_pos_f="cut_pos_$CID" -v header_f="header_$CID" -v FS="" '
  BEGIN{
    while (getline inpl < cut_pos_f){
        cutpos[inpl] = inpl
    }
  }
  $1 ~ /^>/ {
    print $0
  }
  $1 !~ /^>/ {
    for (k=1; k<=NF; k++){
      if (k in cutpos){
      }else{
        printf $k
      }
    }
    print ""
  }
' "${1}_sl_${CID}"

### environmental protection
#
rm -f "${1}_sl_${CID}" cut_pos_$CID
