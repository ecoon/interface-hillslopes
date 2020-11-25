for i in 0.1 0.2 0.5 1 2 5 10; do 
  cp TEMPLATE.sbatch.in manning_n-$i.sbatch
  sed -i "s/MANN/$i/g" manning_n-$i.sbatch
  cp input/manning_n-TEMPLATE.xml.in input/manning_n-$i.xml
  sed -i "s/MANN/$i/g" input/manning_n-$i.xml
  mkdir manning_n-$i
done
