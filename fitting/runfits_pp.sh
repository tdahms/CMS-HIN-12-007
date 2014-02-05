#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

bkgfunctions=(pol1 pol2 pol3 pol4 pol5 pol6)
ptbins=(6.5-30.0 3.0-30.0 3.0-6.5)
rapbins=(0.0-2.4 0.0-1.6 1.6-2.4)
centbins=(0-100)

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${rap} != "1.6-2.4" && ${pt} != "6.5-30.0" ]];
	then
	    echo "skipping " ${rap} " " ${pt};
	    continue;
	fi;
	echo "fitting " ${rap} " " ${pt};
	for bkg in "${bkgfunctions[@]}";
	do
	    for cent in "${centbins[@]}";
	    do
		echo "fitting " ${cent} " " ${bkg};
		./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v signalCB1 signalCB1P ${bkg} -d ${DIRECTORY}/ppFracLogCB -p ${pt} -y ${rap} -t ${cent} -x 4 -r 1 -l 1 -s 0 -n 1 &> ${DIRECTORY}/pp_fitCB_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
		./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -x 4 -r 1 -l 1 -s 0 -n 1 &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}.log
	    done;
	done;
    done;
done;

gzip -f ${DIRECTORY}/pp_fitCB*.log
