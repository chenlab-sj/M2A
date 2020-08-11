#!/bin/bash
# Combine files without duplicating headers.


head -q -n 1 $1  | sort | uniq > Features_Curated.txt
tail -q -n +2 $@ >> Features_Curated.txt 
