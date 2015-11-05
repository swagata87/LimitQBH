#!/bin/bash

function get_number(){
root -l $1 << END

limit->GetEntries()

END
}

cd output_masses

dirlist=`ls | grep Mass | grep output | grep -v txt`

startdir=`pwd`

#dirlist=`echo "Mass_700_input_histos "`

for direct in $dirlist
do

  bobserved=0

  echo ""
  echo "Directory: $direct"
  echo ""

  cd $direct/condor

#   if [ -e observed_limit.root ]
#       then
#       echo "observed_limit.root already exists in $direct "
#       cd $startdir
#       continue
#   fi

#   filelist=`ls | grep higgsCombine | grep root`

#  for file in $filelist
#  do
#    number_entries=`get_number $file | grep "(const Long64_t)" | cut -d ")" -f2`
#    echo $file

#    if [ ${number_entries} == 1 ]
# 	then
# 	mv $file observed_limit.root	
# 	bobserved=1
#    fi

#    if [ $bobserved == 1 ]
# 	then
# 	break
#    fi

#  done


#    Nobserved=`ls | grep higgsCombine | grep root | grep mH50 | head -n 1 | wc -l`

#   if [ $Nobserved -ne 1 ]
#        then
#        echo "$direct: number of observed limit(s) is $Nobserved, should be 1"
#        cd $startdir
#        continue      
#    else
#        file_observed=`ls | grep higgsCombine | grep root | grep mH50 | head -n 1`
#        mv ${file_observed} observed_limit.root
#    fi
      

  if [ -e expected.root ]
      then
      rm -f expected.root
  fi  

  number_expected_files=`ls | grep higgsCombineTest | wc -l`
  if [ ${number_expected_files} -eq 0 ]
      then
      echo "no expected limit files, omitting ..."
      cd $startdir
      continue
  fi

  hadd expected.root higgsCombineTest*.root

  echo ""
  number_entries_tot=`get_number expected.root | grep "(const Long64_t)" | cut -d ")" -f2`
  echo "Number toys total: ${number_entries_tot}"
  echo ""

  cd $startdir


done