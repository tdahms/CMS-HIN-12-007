#! /bin/bash
TODAY=`date +%Y%m%d`
echo "Today is ${TODAY}"

SUFFIX="DblMu0"

DIRECTORY="${TODAY}_${SUFFIX}"

if [ ! -d  ${DIRECTORY} ];
then
    echo "Creating new directory with name ${DIRECTORY}."
    mkdir ${DIRECTORY}
else
    echo "Directory ${DIRECTORY} already exists."
    echo "Continuing anyways, will overwrite files!"
fi

./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4
./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 0-100 -x 4

#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 0-20 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent0-20.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 0-20 -x 4
#
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 20-40 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent20-40.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 20-40 -x 4
#
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 0.0-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 0.0-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 6.5-30.0 -y 1.6-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 6.5-30.0 -y 1.6-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-30.0 -y 1.6-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-30.0 -y 1.6-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v signalCB2WN signalCB2WNP expFunct_HI expFunct_pp  -d simLogCB  -u -p 3.0-6.5 -y 1.6-2.4 -t 40-100 -x 4
#./Fit1DDataPbPbSim -f ../root_files/Data2011_cent40-100.root ../root_files/ppData2011_cent0-100.root -v sigCB2WNG1 sigCB2WNG1P expFunct_HI expFunct_pp -d simLogCBG -u -p 3.0-6.5 -y 1.6-2.4 -t 40-100 -x 4
