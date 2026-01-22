#!/bin/bash
echo What model would you like to analyze?
read model
echo $model selected...
echo Enter submodel start index
read si
echo Enter submodel end index
read ei
echo Bind heatmap true/false
read heat
echo Time average cxn
read timeavg
loc=$(pwd)

echo "NMO,IMO" >> ${model}.si
for ((i=si; i<= ei; i++))
do
python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nux --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nue --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor anue --submodel=${i} --time_average=${timeavg} &&
python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nux --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nue --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor anue --submodel=${i} --time_average=${timeavg} &&
python sigma_average.py $model -p NoTransformation -flavor nux --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p NoTransformation -flavor nue --submodel=${i} --time_average=${timeavg} && python sigma_average.py $model -p NoTransformation -flavor anue --submodel=${i} --time_average=${timeavg} &&
echo Now binding data for error analysis &&
python binder.py -nux ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_AdiabaticMSW_IMO --heatmap=${heat} &&
python binder.py -nux ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_AdiabaticMSW_NMO --heatmap=${heat} &&
python binder.py -nux ./sigmas/${model}_s${i}_NoTransformation_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_NoTransformation_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_NoTransformation_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_NoTransformation --heatmap=${heat} &&
echo -n ${loc}/binders/${model}_s${i}_AdiabaticMSW_NMO_cumsum_unfolded.csv >> ${model}.si
echo ,${loc}/binders/${model}_s${i}_AdiabaticMSW_IMO_cumsum_unfolded.csv >> ${model}.si
done
