for i in $(cat runlist.txt); do 
  cp TEMPLATE.sbatch.in $i.sbatch
  sed -i "s/MYJOBNAME/$i/g" $i.sbatch
  mkdir $i
done
