# simulate scenario where 1% of population splits from another and no migration occurs between them
# calculate allele frequency correlation with python script
# repeat simulation 100 times for different divergence times
# demographic parameters inferred by dadi


#Alex Samano 2024

msms=/usr/local/bin/msms/lib/msms.jar

# 500 gens, no migration
for i in {1..100}
do
java -jar $msms 156 1 -t 382.4826 -I 2 36 120 -ej 0.002217203 2 1 -n 2 0.01  > split_500gens.ms
python3 calculate_shared_afs_from_ms.py split_500gens.ms
done > split_500gens_sims_0.01.txt

# 1000 generations since divergence
for i in {1..100}
do
java -jar $msms 156 1 -t 382.4826 -I 2 36 120 -ej 0.004434405 2 1 -n 2 0.01 > split_1000gens.ms
python3 calculate_shared_afs_from_ms.py split_1000gens.ms
done > split_1000gens_sims_0.01.txt

# 1500 generations since divergence
for i in {1..100}
do
java -jar $msms 156 1 -t 382.4826 -I 2 36 120 -ej 0.006651608 2 1 -n 2 0.01 > split_1500gens.ms
python3 calculate_shared_afs_from_ms.py split_1500gens.ms
done > split_1500gens_sims_0.01.txt

# 2500 generations since divergence
for i in {1..100}
do
java -jar $msms 156 1 -t 382.4826 -I 2 36 120 -ej 0.01108601 2 1 -n 2 0.01 > split_2500gens.ms
python3 calculate_shared_afs_from_ms.py split_2500gens.ms
done > split_2500gens_sims_0.01.txt

# 5000 generations since divergence
for i in {1..100}
do
java -jar $msms 156 1 -t 382.4826 -I 2 36 120 -ej 0.02217202 2 1 -n 2 0.01 > split_5000gens.ms
python3 calculate_shared_afs_from_ms.py split_5000gens.ms
done > split_5000gens_sims_0.01.txt




