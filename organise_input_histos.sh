#!/bin/bash

#input txt for limits

if [ $# -ne 2 ]
then
    echo "You have to provide the approach to limit setting as input:
          shape
          shape_two_bkg
          cut-and-count
          single-bin-counting

          AND

          name of the directory used for the limit calculation
          "
    exit 0
fi

echo "MAKE sure that you edited EmuSpectrum_datacard.txt which is used as input for the limit calculation!!!!!"
sleep 5

mass_min=500
mass_max=1500

mass_binning=500

./create_input_dirs.sh $2 ${mass_min} ${mass_max} ${mass_binning}

if [ "$1" == "shape" ]
then

    ./prepare_run.sh $2 ${mass_min} ${mass_max} ${mass_binning} 

    if ! [ -f EmuSpectrum_datacard.txt ]
	then
	echo "file EmuSpectrum_datacard.txt needed as input is missing. Exiting ..."
	exit 0
    #rm -f EmuSpectrum_datacard.txt
    fi
    
    mass=${mass_min}
    
    while [ $mass -le ${mass_max} ] 
      do

      echo $mass

      if [ -d $2/input_masses/Mass_${mass}_input ]
	  then
	  rm -f $2/input_masses/Mass_${mass}_input/*
      else
	  mkdir $2/input_masses/Mass_${mass}_input
      fi

      
      if ! [ -f root_out/out_mass_${mass}.root ]
	  then
	  echo "Did not find root file with histos for mass point: ${mass}"
	  continue
      else
	  cp root_out/out_mass_${mass}.root $2/input_masses/Mass_${mass}_input
	  cp root_out/out_mass_${mass}.root $2/output_masses/Mass_${mass}_output
	  
      fi
      
      if ! [ -f txt_out/normalization_Mass_${mass}_input_histos.txt ]
	  then
	  echo "Did not find normalization txt file with histos for mass point: ${mass}"
	  continue
      else
	  cp txt_out/normalization_Mass_${mass}_input_histos.txt $2/input_masses/Mass_${mass}_input/normalization.txt
      fi  
      
      if [ -f $2/input_masses/Mass_${mass}_input/EmuSpectrum_datacard.txt ]
	  then
	  rm -f $2/input_masses/Mass_${mass}_input/EmuSpectrum_datacard.txt
      fi
      
      cp EmuSpectrum_datacard.txt $2/input_masses/Mass_${mass}_input/tmp.txt
      
      number_signal_events=$(cat $2/input_masses/Mass_${mass}_input/normalization.txt | grep "signal $mass" | awk '{print $3}')
      num_observed_events=$(cat $2/input_masses/Mass_${mass}_input/normalization.txt | grep "data data_obs" | awk '{print $3}')
      
      echo "number_signal_events = " $number_signal_events
      echo "num_observed_events  = " $num_observed_events

      cat $2/input_masses/Mass_${mass}_input/normalization.txt | grep 'bkg ' | grep -v Up | grep -v Down > tmp.txt
      
  #rates=`echo @rate_signal@`
      rates=`echo ${number_signal_events}`
      echo "rates = " $rates
      while read line
	do
	
    #echo "name bkg: $line" 
	
	rates_tmp=$(echo $line | awk '{print $3}')
	rates=`echo "${rates} ${rates_tmp}"`
	
      done < "tmp.txt"
      
      rm -f tmp.txt
      echo "Now do sed"
      sed -e s/@root_input_file@/out_mass_${mass}.root/g -e  s/@rate_signal@/"${number_signal_events}"/g -e  s/@num_observed_events@/"${num_observed_events}"/g -e  s/@rates@/"${rates}"/g  $2/input_masses/Mass_${mass}_input/tmp.txt >> $2/input_masses/Mass_${mass}_input/EmuSpectrum_datacard.txt
      
      cp $2/input_masses/Mass_${mass}_input/EmuSpectrum_datacard.txt $2/output_masses/Mass_${mass}_output/EmuSpectrum_datacard.txt

      rm -f $2/input_masses/Mass_${mass}_input/tmp.txt
      
#      if [ $mass -lt $mass_int_1 ]
#	  then
#	  let mass=$mass+${mass_binning}
#      elif [ $mass -ge $mass_int_1 ] && [ $mass -lt $mass_int_2 ]
#	  then
#	  let mass=$mass+${mass_binning_2}
#      elif [ $mass -ge $mass_int_2 ]
#	  then
#	  let mass=$mass+${mass_binning_3}
#      fi      


      if [ $mass -le $mass_max ]
         then
         let mass=$mass+${mass_binning}           
      fi

    done
    
    if [ -f nohup.out ]
	then
	rm -f nohup.out
    fi

fi

exit