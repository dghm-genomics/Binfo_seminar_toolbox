#!/bin/bash
## fas-info.sh [OPTIONS FASTA-FILE] / DSMZ / 170421
## check for fasta characteristics
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
set -eu

analyse-fas(){
  awk -v fi_seq="$1" -v fi_verbose="$2" -v OFS="\t" '
    BEGIN{
      if (fi_seq != 1 && fi_verbose != 1){
        print "sequence\tlength\tA\tC\tG\tT\tN (%)\tgap\tbases total (%)\tGC content (%)"
      }
      if (fi_seq != 1 && fi_verbose == 1){
        print "sequence\tlength\tA\tC\tG\tT\ta\tc\tg\tt\tN (%)\tgap\tbases total (%)\tGC content (%)\tR\tY\tW\tS\tM\tK\tH\tB\tV\tD"
      }
      basecount = 0
      A = 0; a = 0; C = 0; c = 0; G = 0; g = 0; T = 0; t = 0; N = 0
      R = 0; Y = 0; W = 0; S = 0; M = 0; K = 0; H = 0; B = 0; V = 0; D = 0
      gap = 0
    }
    $1 !~ ">" {
          if (seql == 0){
            seql =  $1
          }else{
            seql = seql $1
          }
    }
    $1 ~ ">" && NR == 1 {
      sub(/>/,"",$0)
      sheader = $0
      seql = 0
    }
    $1 ~ ">" && NR != 1 {
      if(fi_seq == 0){
        split (seql,sequencev,"")
        for (i in sequencev){
          if (sequencev[i] == "A"){A++}
          if (sequencev[i] == "a"){a++}
          if (sequencev[i] == "C"){C++}
          if (sequencev[i] == "c"){c++}
          if (sequencev[i] == "G"){G++}
          if (sequencev[i] == "g"){g++}
          if (sequencev[i] == "T"){T++}
          if (sequencev[i] == "t"){t++}
          if (sequencev[i] == "N"){N++}
          if (fi_verbose == 1){
            if (sequencev[i] == "R"){R++}
            if (sequencev[i] == "Y"){Y++}
            if (sequencev[i] == "W"){W++}
            if (sequencev[i] == "S"){S++}
            if (sequencev[i] == "M"){M++}
            if (sequencev[i] == "K"){K++}
            if (sequencev[i] == "H"){H++}
            if (sequencev[i] == "B"){B++}
            if (sequencev[i] == "V"){V++}
            if (sequencev[i] == "D"){D++}
          }
          if (sequencev[i] == "-"){gap++}
          basecount++
        }
        if (fi_seq == 0 && fi_verbose == 0){
          print sheader,basecount,A+a,C+c,G+g,T+t,N" ("N/basecount*100")",gap,A+C+G+T+a+c+g+t" ("(A+C+G+T+a+c+g+t)/basecount*100")",C+c+G+g" ("(C+c+G+g)/basecount*100")"
        }
        if (fi_seq == 0 && fi_verbose == 1){
          print sheader,basecount,A,C,G,T,a,c,g,t,N" ("N/basecount*100")",gap,A+C+G+T+a+c+g+t" ("(A+C+G+T+a+c+g+t)/basecount*100")",C+c+G+g" ("(C+c+G+g)/basecount*100")",R,Y,W,S,M,K,H,B,V,D
        }
        sub(/>/,"",$0)
        sheader = $0
        basecount = 0
        A = 0; a = 0; C = 0; c = 0; G = 0; g = 0; T = 0; t = 0; N = 0
        R = 0; Y = 0; W = 0; S = 0; M = 0; K = 0; H = 0; B = 0; V = 0; D = 0
        gap = 0
      }else{
        print seql
      }
      seql = 0
    }END{
      if(fi_seq == 0){
        split (seql,sequencev,"")
        for (i in sequencev){
          if (sequencev[i] == "A"){A++}
          if (sequencev[i] == "a"){a++}
          if (sequencev[i] == "C"){C++}
          if (sequencev[i] == "c"){c++}
          if (sequencev[i] == "G"){G++}
          if (sequencev[i] == "g"){g++}
          if (sequencev[i] == "T"){T++}
          if (sequencev[i] == "t"){t++}
          if (sequencev[i] == "N"){N++}
          if (fi_verbose == 1){
            if (sequencev[i] == "R"){R++}
            if (sequencev[i] == "Y"){Y++}
            if (sequencev[i] == "W"){W++}
            if (sequencev[i] == "S"){S++}
            if (sequencev[i] == "M"){M++}
            if (sequencev[i] == "K"){K++}
            if (sequencev[i] == "H"){H++}
            if (sequencev[i] == "B"){B++}
            if (sequencev[i] == "V"){V++}
            if (sequencev[i] == "D"){D++}
          }
          if (sequencev[i] == "-"){gap++}
          basecount++
        }
        if (fi_seq == 0 && fi_verbose == 0){
          print sheader,basecount,A+a,C+c,G+g,T+t,N" ("N/basecount*100")",gap,A+C+G+T+a+c+g+t" ("(A+C+G+T+a+c+g+t)/basecount*100")",C+c+G+g" ("(C+c+G+g)/basecount*100")"
        }
        if (fi_seq == 0 && fi_verbose == 1){
          print sheader,basecount,A,C,G,T,a,c,g,t,N" ("N/basecount*100")",gap,A+C+G+T+a+c+g+t" ("(A+C+G+T+a+c+g+t)/basecount*100")",C+c+G+g" ("(C+c+G+g)/basecount*100")",R,Y,W,S,M,K,H,B,V,D
        }
        sub(/>/,"",$0)
        sheader = $0
        basecount = 0
        A = 0; a = 0; C = 0; c = 0; G = 0; g = 0; T = 0; t = 0; N = 0
        R = 0; Y = 0; W = 0; S = 0; M = 0; K = 0; H = 0; B = 0; V = 0; D = 0
        gap = 0
      }else{
        print seql
      }
      seql = 0
    }
  ' "$3"
}

### create help panel
#
if [ $# -eq 0 ]
then
    echo -e "
\033[1mProgram\033[0m
    \tfas-info.sh (${VERSION}) - print fasta attributes
\033[1mSYNOPSIS\033[0m
    \tfas-info.sh [OPTION]... [FASTA-FILE]
\033[1mDESCRIPTION\033[0m\

  \tprint fasta characteristics
  
  \t\033[1m-c\033[0m\n\t  count sequences
  
  \t\033[1m-h\033[0m\n\t  print header only

  \t\033[1m-s\033[0m\n\t  print sequence only

  \t\033[1m-v\033[0m\n\t  verbose

  \t Results will be promted to std out
  "
  exit
fi

### read user options
#
fi_verbose=0
fi_header=0
fi_sequence=0
fi_count=0
while getopts chsv opt # 2>/dev/null
do
   case $opt in
      c) fi_count=1;;
      h) fi_header=1;;
      s) fi_sequence=1;;
      v) fi_verbose=1;;
      *) echo "\n [fas-info] ERROR invalid user option!\n"; exit 1;;
   esac
done
shift $((OPTIND-1))





### main section
#
if [ $fi_header == 1 ]
then
    grep ">" "$1"
    exit;
fi

if [ $fi_count == 1 ]
then
    grep -c ">" "$1"
    exit
fi

analyse-fas "$fi_sequence" "$fi_verbose" "$1"
