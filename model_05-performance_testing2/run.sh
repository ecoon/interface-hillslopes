for i in $(cat runlist.txt); do 
  sbatch $i.sbatch
done
