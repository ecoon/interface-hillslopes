for i in $(cat runlist.txt); do 
  echo $i $(( $(date -d "`tail -n 1 ${i}/out.log`" +%s) - $(date -d "`head -n 1 ${i}/out.log`" +%s) ))
done
