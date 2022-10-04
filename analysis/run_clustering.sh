module load imp
module load python3/pyrmsd

export cl=3

cp ../A_models_clust${cl}.txt ../scoresA.txt
cp ../B_models_clust${cl}.txt ../scoresB.txt
cp ../../density.txt .
nohup python imp-sampcon/pyext/src/exhaust.py \
       -n Pks13 -p ../ -ra A_models_clust${cl}.rmf3 -rb B_models_clust${cl}.rmf3 -d density.txt \
       -m cpu_omp -c 8 -g 3.0 > clustering.log &
