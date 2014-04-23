#! /bin/bash
DIRECTORY="20140320_M2242_DblMu0/"

ptbins=(65-30 3-30 3-65)
rapbins=(0-24 0-16 16-24)
centbins=(0-100 0-20 20-40 40-100)
bkg=""
expBkg=""

chi2=(0 0 0 0 0)
dof=(0 0 0 0 0)

i=0

echo "pp"

for pt in "${ptbins[@]}";
do
    for rap in "${rapbins[@]}";
    do
	if [[ ${pt} == "65-30" && ${rap}  == "0-24" ]];
	then 
	    bkg="pol7";
	    expBkg="expPol1";
	elif [[ ${pt} == "65-30" && ${rap}  == "0-16" ]]
	then
	    bkg="pol3";
	    expBkg="expPol3";
	elif [[ ${pt} == "65-30" && ${rap}  == "16-24" ]];
	then
	    bkg="pol1";
	    expBkg="expPol1";
	elif [[ ${pt} == "3-30" && ${rap}  == "16-24" ]];
	then
	    bkg="pol1";
	    expBkg="expPol3";
	elif [[ ${pt} == "3-65" && ${rap}  == "16-24" ]];
	then
	    bkg="pol1";
	    expBkg="expPol1";
	else
	    continue;
	fi;
	chi2[i]=`grep chi2 ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent0-100_fitResult.txt | awk '{print $2}'`
	dof[i]=`grep DOF ${DIRECTORY}/ppFracLogCBG_${bkg}_rap${rap}_pT${pt}_cent0-100_fitResult.txt | awk '{print $2}'`
	echi2[i]=`grep chi2 ${DIRECTORY}/ppFracLogCBG_${expBkg}_rap${rap}_pT${pt}_cent0-100_fitResult.txt | awk '{print $2}'`
	edof[i]=`grep DOF ${DIRECTORY}/ppFracLogCBG_${expBkg}_rap${rap}_pT${pt}_cent0-100_fitResult.txt | awk '{print $2}'`
	i=$((${i}+1))
    done;
done;

#echo -n "chi2: "

for i in $(seq 0 $((${#chi2[@]} - 1)))
do
#    echo "${i}";
    echo -n "${chi2[${i}]} "
done;

echo ""
#echo -n "dof: "

for i in $(seq 0 $((${#dof[@]} - 1)))
do
    echo -n "${dof[${i}]} "
done;

echo ""
echo ""
echo ""
#echo -n "chi2: "

for i in $(seq 0 $((${#echi2[@]} - 1)))
do
    echo -ne "${echi2[${i}]} "
done;

echo ""
#echo -n "dof: "

for i in $(seq 0 $((${#edof[@]} - 1)))
do
    echo -ne "${edof[${i}]} "
done;
echo ""


echo "----------"
echo "PbPb"
i=0

for cent in "${centbins[@]}"
do
    for rap in "${rapbins[@]}";
    do
	for pt in "${ptbins[@]}";
	do
	if [[ ${pt} == "65-30" && ${rap}  == "0-24" ]];
	then 
	    if [[ ${cent} == "0-100" ]];
	    then
		bkg="pol3";
		expBkg="expPol3";
	    elif [[ ${cent} == "0-20" ]];
	    then
		bkg="pol1";
		expBkg="expPol1";
	    elif [[ ${cent} == "20-40" ]];
	    then
		bkg="pol2";
		expBkg="expPol2";
	    elif [[ ${cent} == "40-100" ]];
	    then
		bkg="pol0";
		expBkg="expPol1";
	    fi;
	elif [[ ${pt} == "65-30" && ${rap}  == "0-16" ]]
	then
	    if [[ ${cent} == "0-100" ]];
	    then
		bkg="pol1";
		expBkg="expPol1";
	    elif [[ ${cent} == "0-20" ]];
	    then
		bkg="pol1";
		expBkg="expPol1";
	    elif [[ ${cent} == "20-40" ]];
	    then
		bkg="pol1";
		expBkg="expPol1";
	    elif [[ ${cent} == "40-100" ]];
	    then
		bkg="pol0";
		expBkg="expPol1";
	    fi;
	elif [[ ${pt} == "65-30" && ${rap}  == "16-24" && ${cent} == "0-100" ]];
	then
	    bkg="pol3";
	    expBkg="expPol2";
	elif [[ ${pt} == "3-30" && ${rap}  == "16-24" ]];
	then
	    if [[ ${cent} == "0-100" ]];
	    then
		bkg="pol3";
		expBkg="expPol2";
	    elif [[ ${cent} == "0-20" ]];
	    then
		bkg="pol3";
		expBkg="expPol3";
	    elif [[ ${cent} == "20-40" ]];
	    then
		bkg="pol2";
		expBkg="expPol2";
	    elif [[ ${cent} == "40-100" ]];
	    then
		bkg="pol1";
		expBkg="expPol1";
	    fi;
	elif [[ ${pt} == "3-65" && ${rap}  == "16-24" && ${cent} == "0-100" ]];
	then
	    bkg="pol2";
	    expBkg="expPol2";
	else
	    continue;
	fi;
	echo ${rap} " " ${pt} " " ${cent}
	chi2[i]=`grep chi2 ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_fitResult.txt | awk '{print $2}'`
	dof[i]=`grep DOF ${DIRECTORY}/fracLogCBG_${bkg}_rap${rap}_pT${pt}_cent${cent}_fitResult.txt | awk '{print $2}'`
	echi2[i]=`grep chi2 ${DIRECTORY}/fracLogCBG_${expBkg}_rap${rap}_pT${pt}_cent${cent}_fitResult.txt | awk '{print $2}'`
	edof[i]=`grep DOF ${DIRECTORY}/fracLogCBG_${expBkg}_rap${rap}_pT${pt}_cent${cent}_fitResult.txt | awk '{print $2}'`
	i=$((${i}+1))
	done;
    done;
done;

#echo -n "chi2: "

for i in $(seq 0 $((${#chi2[@]} - 1)))
do
#    echo "${i}";
    echo -n "${chi2[${i}]} "
done;

echo ""
#echo -n "dof: "

for i in $(seq 0 $((${#dof[@]} - 1)))
do
    echo -n "${dof[${i}]} "
done;

echo ""
echo ""
echo ""
#echo -n "chi2: "

for i in $(seq 0 $((${#echi2[@]} - 1)))
do
    echo -ne "${echi2[${i}]} "
done;

echo ""
#echo -n "dof: "

for i in $(seq 0 $((${#edof[@]} - 1)))
do
    echo -ne "${edof[${i}]} "
done;
echo ""
