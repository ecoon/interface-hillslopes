for i in 0.1 0.2 0.5 1 2 5 10; do 
  sbatch manning_n-$i.sbatch
done
