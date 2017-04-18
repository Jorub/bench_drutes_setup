#!/bin/bash

for i in $(seq 1 $1); do
	if [ ! -d "$i" ]; then
 		mkdir $i
	fi
	cp -a drutemp/* $i
done

