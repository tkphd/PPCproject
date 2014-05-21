#!/bin/bash

# configs  change to your desired configurations:
#N ; number of node
#n ; number of rank
#p ; number of thread

#configs=(
#"32     32     64 "
#"32     64     32 "
#"32     128    16 "
#"32     256     8 "
#"32     512     4 "
#"32     2048    1 "
#"64     64     64 "
#"64     128    32 "
#"64     256    16 "
#"64     512     8 "
#"64     1024    4 "
#"64     4096    1 "
#"128    128    64 "
#"128    256    32 "
#"128    512    16 "
#"128    1024    8 "
#"128    2048    4 "
#"128    8192    1 "
#"256    256    64 "
#"256    512    32 "
#"256    1024   16 "
#"256    2048    8 "
#"256    4096    4 "
#"256    16384   1 "
#)

configs=(
"16     16     16 "
"16     32      8 "
"16     64      4 "
"16     128     2 "
"16     256     1 "
"32     32     16 "
"32     64      8 "
"32     128     4 "
"32     256     2 "
"32     512     1 "
"64     64     16 "
"64     128     8 "
"64     256     4 "
"64     512     2 "
"64     1024    1 "
"128    128    16 "
"128    256     8 "
"128    512     4 "
"128    1024    2 "
"128    2048    1 "
"256    256    16 "
"256    512     8 "
"256    1024    4 "
"256    2048    2 "
"256    4096    1 "
"512    512    16 "
"512    1024    8 "
"512    2048    4 "
"512    4096    2 "
"512    8192    1 "
)

rm [mp][cf].*.*.sh

for conf in "${configs[@]}" ; do
    nnode=`echo $conf | awk '{print $1}'`
    nrank=`echo $conf | awk '{print $2}'`
    nthread=`echo $conf | awk '{print $3}'`

    if [ $nnode -ge 128 ]; then
        partition="medium"
    else
        partition="small"
    fi

    sed -e "s/PARTITION/$partition/g" \
        -e "s/NNODE/$nnode/g" \
        -e "s/NRANK/$nrank/g" \
        -e "s/NTHREAD/$nthread/g" \
        srun_MC.tmp > mc.$nnode.$nrank.sh


    sed -e "s/PARTITION/$partition/g" \
        -e "s/NNODE/$nnode/g" \
        -e "s/NRANK/$nrank/g" \
        -e "s/NTHREAD/$nthread/g" \
        srun.tmp > pf.$nnode.$nrank.sh
done
