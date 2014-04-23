#! /bin/bash

DATEIN="20140423"
SUFFIXIN="SimFits_M2242_DblMu0_AllCent_WithSyst_final"
DIRECTORYIN="${DATEIN}_${SUFFIXIN}"

DATEOUT="20140423"
SUFFIXOUT="CFscan_WithSyst"
DIRECTORYOUT="${DATEOUT}_${SUFFIXOUT}"

if [ ! -d  ${DIRECTORYIN} ];
then
    echo "${DIRECTORYIN} does not exist."
    exit 1
else
    echo "Reading results in ${DIRECTORYIN}"
fi

if [ ! -d  ${DIRECTORYOUT} ];
then
    echo "Creating ${DIRECTORYOUT}:"
    mkdir ${DIRECTORYOUT}
else
    echo "${DIRECTORYOUT} already exist, change name of output directory!"
    exit 1
fi

ptbins=(3-30 65-30)
#rapbins=(0-16 16-24)
rapbins=(0-16)
#centbins=(HI020 HI2040 HI40100)
centbins=(HI020)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${pt} == "65-30" && ${rap} == "0-16" ]];
	then
	    bkg="PbPbpol1_HI020_pol1_HI2040_pol0_HI40100_pppol3";
	elif [[ ${pt} == "3-30" && ${rap} == "16-24" ]];
	then
	    bkg="PbPbpol3_HI020_pol2_HI2040_pol1_HI40100_pppol1";
	else
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	nsteps=10
	rmax=0.1
	if [[ ${rap} == "16-24" && ${pt} == "3-30" ]];
	then
	    ntoys=1000
	    nsteps=11
	    rmin=0.0
	    rmax=2.0
	elif [[ ${rap} == "0-16" && ${pt} == "65-30" ]];
	then
	    ntoys=1000
	    nsteps=11
	    rmin=0.0
	    rmax=1.0
	fi;
	for cent in "${centbins[@]}";
	do
	    if [[ ${rap} == "16-24" && ${pt} == "3-30" ]];
	    then
		if [[ ${cent} == "HI020" ]];
		then
		    ntoys=1000
		    nsteps=15
		    rmin=0.8
		    rmax=3.6
		elif [[ ${cent} == "HI2040" ]];
		then
		    ntoys=1000
		    nsteps=12
		    rmin=0.0
		    rmax=2.2
		elif [[ ${cent} == "HI40100" ]];
		then
		    ntoys=1000
		    nsteps=11
		    rmin=0.0
		    rmax=2.0
		fi
	    elif [[ ${rap} == "0-16" && ${pt} == "65-30" ]];
	    then
		if [[  ${cent} == "HI020" ]];
		then
		    ntoys=1000
		    nsteps=21
		    rmin=0.0
		    rmax=1.0
		elif [[  ${cent} == "HI2040" ]];
		then
		    ntoys=1000
		    nsteps=15
		    rmin=0.0
		    rmax=1.4
		elif [[ ${cent} == "HI40100" ]];
		then
		    ntoys=1000
		    nsteps=13
		    rmin=0.0
		    rmax=0.6
		fi;
	    fi;
	    if [ ! -d  ${DIRECTORYOUT}/${cent} ];
	    then	    
		mkdir "${DIRECTORYOUT}/${cent}"
	    fi;
	    echo "Running F-C scan for rap${rap} pT${pt} ${bkg} ${cent}."
	    echo "Scanning with: ntoys = ${ntoys} in ${nsteps} steps from ${rmin} to ${rmax}:"
	    root -b -q StandardHypoTestInvDemo.C+\(\"${DIRECTORYIN}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_centMult_Workspace.root\",\"workspace\",\"model_${cent}\",\"B_only_model_${cent}\",\"redDataSim\",0,2,false,${nsteps},${rmin},${rmax},${ntoys},false\) &> ${DIRECTORYOUT}/rap${rap}_pt${pt}_cent${cent}.log 
	    gzip "${DIRECTORYOUT}/rap${rap}_pt${pt}_cent${cent}.log"
	    mv ${DIRECTORYOUT}/*.* ${DIRECTORYOUT}/${cent}/
	    echo "Done with F-C scan for rap${rap} pT${pt} ${bkg} ${cent}."
	    sleep 30
	done;
    done;
done;
