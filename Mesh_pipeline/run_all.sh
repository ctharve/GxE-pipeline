#!/bin/bash

plate=$1

#python run_analysis.py D1P6_KP39351_T6C1.txt

ls ${plate}_*.txt | while read f; do python run_analysis.py $f; done
