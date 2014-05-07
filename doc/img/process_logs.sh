#!/bin/bash
rm pr_*.dat
rm prMC_*.dat

rsync cci:~/barn/project/pr*.log ./

for f in proj_32_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "32\t" >>  pr_N32.dat
		cat $f >> pr_N32.dat
	fi
done
for f in prMC_32_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "32\t" >>  prMC_N32.dat
		cat $f >> prMC_N32.dat
	fi
done
#
for f in proj_64_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "64\t" >>  pr_N64.dat
		cat $f >> pr_N64.dat
	fi
done
for f in prMC_64_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "64\t" >>  prMC_N64.dat
		cat $f >> prMC_N64.dat
	fi
done
#
for f in proj_128_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "128\t" >>  pr_N128.dat
		cat $f >> pr_N128.dat
	fi
done
for f in prMC_128_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "128\t" >>  prMC_N128.dat
		cat $f >> prMC_N128.dat
	fi
done
#
for f in proj_256_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "256\t" >>  pr_N256.dat
		cat $f >> pr_N256.dat
	fi
done
for f in prMC_256_*.log
do
	len=`wc -l $f | awk '{print $1}'`
	if [ len=="1" ]
	then
		echo -ne "256\t" >>  prMC_N256.dat
		cat $f >> prMC_N256.dat
	fi
done

for f in *.plt; do gnuplot $f; done
