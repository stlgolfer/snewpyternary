#!/bin/bash
echo What model would you like to analyze?
read model
echo $model selected...
echo Enter submodel start index
read si
echo Enter submodel end index
read ei
echo Bind heatma t/f
read heat

for ((i=si; i<= ei; i++))
do
python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nux --submodel=${i} && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nue --submodel=${i} && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor anue --submodel=${i} &&
python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nux --submodel=${i} && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nue --submodel=${i} && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor anue --submodel=${i} &&
echo Now binding data for error analysis &&
python binder.py -nux ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_AdiabaticMSW_IMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_AdiabaticMSW_IMO --heatmap=${heat} &&
python binder.py -nux ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s${i}_AdiabaticMSW_NMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s${i}_AdiabaticMSW_NMO --heatmap=${heat}
done
