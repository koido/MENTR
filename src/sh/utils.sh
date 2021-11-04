#!/bin/bash

# == help ==
function help(){
  inf=$1
  awk '{
    if(NR>2){
      if(/^#/){
        sub("^# ?", "");
        print
      }else{
        exit
      }
    }
  }' $inf
  exit 1
}

# == INT check ==
function int_chk(){
  x=$1
  name=$2
  expr ${x} + 1 >/dev/null 2>&1

  if [ $? -ne 1 -a $? -ne 0 ]; then
    echo -e "Error: ${name} must be INT.\n"
    help
  elif [ ${x} -lt 1 ]; then
    echo -e "Error: ${name} must be >=1\n"
    help
  fi
}
