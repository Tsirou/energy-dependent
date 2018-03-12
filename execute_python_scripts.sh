#!/bin/bash

total_nb_analyses=9

analyses_folders=("HiRes_4-30TeV" "HiRes_4-10TeV" "HiRes_10-20TeV" "4-30TeV" "4-10TeV" "10-20TeV" "2-4_4-8TeV" "4-8_9-6TeV" "gt9-6TeV")

save_path="/home/tsirou/Documents/Analyses/Results/Energy-dependent/MultiEnergyBins/2d0_EHE"

no_norm=("no_cmap_norm" "")


normed=0

for (( j=0; j<$total_nb_analyses; j++ ))
do
    python2.7 main.py $j
   
    cp $save_path/*.txt $save_path/${no_norm[$normed]}/${analyses_olders[$j]}/

    cp $save_path/*.png $save_path/${no_norm[$normed]}/${analyses_olders[$j]}/
    
    echo "$(($j+1)) out of $total_nb_analyses done --- $((($j+1)*100./$total_nb_analyses)) % ---"

done


