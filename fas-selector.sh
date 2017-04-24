#!/bin/bash
## fas-selector.sh FASTA-FILE / DSMZ / 170421
## selection and deletion of sequences in a multi-fasta file
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
readonly DATE
set -eu

### realising sequence selection specifying sequence number
#
num_select(){
  awk -v inpl="$1" -v sdel="$2" -v seq_count="$3" 'BEGIN{
          split (inpl,num_inpv)
          for (j in num_inpv){
            numv[num_inpv[j]]
          }
    }
    $1 !~ ">" {
        i++
        sel_seq[i] =  $1
    }
    $1 ~ ">" && NR == 1 {
        n++
        sel_header = $0
        if(n in numv){
           sel_flag = 1
        }else{
           sel_flag = 0
        }
    }
    $1 ~ ">" && NR != 1 {
        if (sel_flag == 1){
            if(sdel == 0){
              print sel_header
              for (j=1; j<=i; j++){
                  print sel_seq[j]
              }
            }
        }else{
            if(sdel == 1){
              print sel_header
              for (j=1; j<=i; j++){
                print sel_seq[j]
              }
            }
        }
        i = 0
        delete sel_seq
        n++
        sel_header = $0
        if(n in numv){
            sel_flag = 1
        }else{
            sel_flag = 0
        }
    }
    END{
        if (sel_flag == 1){
            if(sdel == 0){
              print ">"sel_header
              for (j=1; j<=i; j++){
                  print sel_seq[j]
              }
            }
        }else{
            if(sdel == 1){
              print sel_header
              for (j=1; j<=i; j++){
                 print sel_seq[j]
              }
            }
        }
        i = 0
        delete sel_seq
    }
  ' "$4"
}

### realising sequence selection specifying fasta header(s)
#
id_select(){
  awk -v header_inp="$1" -v sdel="$2" 'BEGIN{
          split (header_inp,header_set)
       }
       $1 !~ ">" {
         if (sel_flag == 1){
           i++
           sel_seq[i] =  $1
         }else{
           delete sel_seq
         }
       }
       $1 ~ ">" && NR == 1 {
         if (sdel == 1){
           sel_flag = 1
           sel_header = $0
         }else{
           sel_flag = 0
         }
         for (h in header_set){
           if($0 ~ header_set[h]){
             if (sdel == 1){
               sel_flag = 0
             }else{
               sel_flag = 1
               sel_header = $0
             }
           }
         }
       }
       $1 ~ ">" && NR != 1 {
          if (sel_flag == 1){
            print sel_header
            for (j=1; j<=i; j++){
              print sel_seq[j]
            }
            i=0
          }
          if (sdel == 1){
            sel_flag = 1
            sel_header = $0
          }else{
            sel_flag = 0
          }
          for (h in header_set){
            if($0 ~ header_set[h]){
              if (sdel == 1){
                sel_flag = 0
              }else{
                sel_flag = 1
                sel_header = $0
              }
            }
          }
         }
        END{
          if (sel_flag == 1){
            print sel_header
            for (j=1; j<=i; j++){
              print sel_seq[j]
            }
            i=0
            delete sel_seq
          }
        }
  ' "$3"
}

### create help panel
#
if [ $# == 0 ]
then
  echo -e "
\033[1mProgram\033[0m
    \tfas-selector.sh (${VERSION}) - select or delete sequences from a multi-fasta file
\033[1mSYNOPSIS\033[0m
    \tfas-selector.sh [OPTIONS] [FASTA-FILE]
\033[1mDESCRIPTION\033[0m
    
    \t\033[1m-d\033[0m\n\t  delete selection
    
    \t\033[1m-f\033[0m\n\t  selection has to be specified by id(s) provided in a file
    
    \t\033[1m-h\033[0m\n\t  selection has to be specified by id(s)

    \t\033[1m-n\033[0m\n\t  selection has to be specified by line numbers
    
    \t\033[1m-N\033[0m\n\t  selection has to be specified by line numbers provided in a file
    
    
    \t Results will be promted to std out
    "
    exit
fi

### read user options
#
sel_del=0
sel_header=
sel_file=
sel_num=
sel_num_file=
while getopts df:h:n:N: opt # 2>/dev/null
do
   case $opt in
      d) sel_del=1;;
      f) sel_file="$OPTARG";;
      h) sel_header="$OPTARG";;
      n) sel_num="$OPTARG";;
      N) sel_num_file="$OPTARG";;
      *) echo "\n [call-vcvs.sh] ERROR invalid user option!\n"; exit 1;;
   esac
done
shift $((OPTIND-1))



### main section
#
if [ ! -z "$sel_header" ]
then
    id_select "$sel_header" "$sel_del" "$1"
fi

if [ ! -z "$sel_file" ]
then
    sheader=$(tr '\n' ' ' < "$sel_file")
    id_select "$sheader" "$sel_del" "$1"
fi

if [ ! -z "$sel_num_file" ]
then
    seq_count=$(grep -c ">" "$1")
    snum=$(tr '\n' ' ' < "$sel_num_file")
    num_select "$snum" "$sel_del" "$seq_count" "$1"
fi

if [ ! -z "$sel_num" ]
then
    seq_count=$(grep -c ">" "$1")
    num_select "$sel_num" "$sel_del" "$seq_count" "$1"
fi
