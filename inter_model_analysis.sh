#!/bin/bash
#echo What model would you like to analyze?
#read model
#echo $model selected...
echo Bind heatma t/f
read heat

for model in Kuroda_2020 Nakazato_2013 Tamborra_2014 Walk_2018 Walk_2019 Warren_2020 Zha_2021
do
python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nux && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor nue && python sigma_average.py $model -p AdiabaticMSW_IMO -flavor anue &&
python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nux && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor nue && python sigma_average.py $model -p AdiabaticMSW_NMO -flavor anue &&
echo Now binding data for error analysis &&
python binder.py -nux ./sigmas/${model}_s0_AdiabaticMSW_IMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s0_AdiabaticMSW_IMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s0_AdiabaticMSW_IMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s0_AdiabaticMSW_IMO --heatmap=${heat} &&
python binder.py -nux ./sigmas/${model}_s0_AdiabaticMSW_NMO_BstChnl_nu_mu_sigma_average.csv -nue ./sigmas/${model}_s0_AdiabaticMSW_NMO_BstChnl_nu_e_sigma_average.csv -anue ./sigmas/${model}_s0_AdiabaticMSW_NMO_BstChnl_nu_e_bar_sigma_average.csv --title=${model}_s0_AdiabaticMSW_NMO --heatmap=${heat}
done
