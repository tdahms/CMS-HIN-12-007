#! /bin/bash

TODAY="20140820"
SUFFIX="M2242_DblMu0_Split"

DIRECTORY="${TODAY}_${SUFFIX}"


#bkgfunctions=(pol0 pol1 pol2 pol3 pol4 pol5 pol6 expPol1 expPol2 expPol3 expPol4 expPol5 expPol6 expPol7)
ptbins=(65-30 3-30)
rapbins=(0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)

ratio_tmp=0;
ratio=0;
relRatio=0;
min=0;
max=0;

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ (${rap} == "16-24" && ${pt} == "65-30")|| 
	      (${rap} != "16-24" && ${pt} == "3-30")]];
	then
	    continue;
	fi;
	if [[ ${rap} == "0-24" ]];
	then
	    rapAlt="0.0-2.4"
	elif [[ ${rap} == "0-16" ]];
	then
	    rapAlt="0.0-1.6"
	elif [[ ${rap} == "16-24" ]];
	then
	    rapAlt="1.6-2.4"
	fi;	    
	if [[ ${pt} == "65-30" ]]
	then
	    ptAlt="6.5-30.0"
	elif [[ ${pt} == "3-30" ]]
	then
	    ptAlt="3.0-30.0"
	fi;
        echo "pp:"
	    echo ${DIRECTORY}/pp_fitCBG_pol*_rap${rapAlt}_pt${ptAlt}_cent0-100.log.gz
	    zgrep 'R_psi(2S) ' ${DIRECTORY}/pp_fitCBG_pol*_rap${rapAlt}_pt${ptAlt}_cent0-100.log.gz | awk '{print $3; print $5}'
        for cent in "${centbins[@]}";
        do
	    echo ${DIRECTORY}/pbpb_fitCBG_pol*_rap${rapAlt}_pt${ptAlt}_cent${cent}.log.gz
	    zgrep 'R_psi(2S) ' ${DIRECTORY}/pbpb_fitCBG_pol*_rap${rapAlt}_pt${ptAlt}_cent${cent}.log.gz | awk '{print $3; print $5}'
        done;
    done;
done;
