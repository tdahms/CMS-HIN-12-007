#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="M2242_FixCBtail_NoSumw2_DblMu0_NoMinos_finalCBtailFits_fixWideFactorPbPbMC_corrPsi2SError_systStudies"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi


pt="6.5-30.0"
rap="0.0-2.4"
bkg="pol5"
cent="0-100"
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log

pt="6.5-30.0"
rap="0.0-1.6"
bkg="pol3"
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log

pt="6.5-30.0"
rap="1.6-2.4"
bkg="pol5"
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log

pt="3.0-30.0"
rap="1.6-2.4"
bkg="pol5"
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log

pt="3.0-6.5"
rap="1.6-2.4"
bkg="pol1"
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -a &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeAlpha.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -n &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeN.log
./Fit1DDataPbPb -f ../root_files/ppData2013_DblMu0_cent${cent}_M22-42.root -v sigCB1G2 sigCB1G2P ${bkg} -d ${DIRECTORY}/ppFracLogCBG -p ${pt} -y ${rap} -t ${cent} -r 1 -l 1 -s 0 -b 1 -g &> ${DIRECTORY}/pp_fitCBG_${bkg}_rap${rap}_pt${pt}_cent${cent}_freeGwidth.log

gzip -f ${DIRECTORY}/pp_fitCB*.log
