#!/bin/bash
#echo $0 $1

#make -f Makefile.linux

rm a.csv
rm a.dat

#for nodeid in 1 2 3 4 5 6 7
for nodeid in 1111
do

	if [ $nodeid -ge 10 ]; then
	  subfolder=100$nodeid
	else
	  subfolder=1000$nodeid
	fi
	./arrtime $HOME/data/MSIMultiple/$subfolder $subfolder arrival.config
done
