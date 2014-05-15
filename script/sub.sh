#!/bin/bash

checkStatus(){
    while [ '1' == '1' ]; do
        mesg=`squeue -u GGSTkllt | grep GGSTkllt`
        if [ "$mesg" == ''  ]; then
            break
        fi
        sleep 2m
    done

}

for i in `ls [mp][cf]*.sh | sort -t '.' -nk 2 `; do
    # skip finished jobs
    flog=`grep '#SBATCH -o' $i | awk '{print $3}'`
    echo -ne "${i}\t"
    if [ -f $flog ]; then
        echo "   Skipped"
        continue;
    fi

    sbatch $i

    checkStatus
done

