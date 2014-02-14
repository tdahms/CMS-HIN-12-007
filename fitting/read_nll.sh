#!/bin/bash

ptbins=(65-30 3-30 3-65)
rapbins=(0-24 0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)

DIRECTORY="20140210_M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError_bkgStudy/"
echo "pp"
for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${rap} != "16-24" && ${pt} != "65-30" ]];
	then
	    continue;
	fi;
	for cent in "${centbins[0]}";
	do
	echo "rap${rap}_pT${pt}_cent${cent}"
	grep NLL ${DIRECTORY}ppFracLogCBG_*_rap${rap}_pT${pt}_cent${cent}.txt | awk '{print $2}'
	echo $'\n\n'
	done;
    done;
done;

echo "-------------------------------"
echo $'\n\n'
echo "PbPb"

for cent in "${centbins[@]}";
do
    for pt in "${ptbins[@]}";
    do
	for rap in "${rapbins[@]}";
	do
	    if [[ ${rap} != "16-24" && ${pt} != "65-30" ]];
	    then
		continue;
	    fi;

	    if [[ ${cent} != "0-100" && ${pt} != "3-30" && ${rap} == "16-24" ]];
	    then
		continue;
	    fi;
	    echo "rap${rap}_pT${pt}_cent${cent}"
	    grep NLL ${DIRECTORY}fracLogCBG_*_rap${rap}_pT${pt}_cent${cent}.txt | awk '{print $2}'
	    echo $'\n\n'
	done;
    done;
done;

