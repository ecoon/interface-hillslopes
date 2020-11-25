for i in 0.1 0.2 0.5 1 2 5 10; do 
    rm -f ./manning_n-$i.sbatch
    rm -f ./input/manning_n-$i.xml
    rm -rf ./manning_n-$i/
done

