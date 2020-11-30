for i in $(cat runlist.txt); do 
    rm -f $i.sbatch
    rm -rf $i
done

