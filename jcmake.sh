#!/bin/bash

rm joaninhasCalorosas ; make -k joaninhasCalorosas ; /usr/bin/time ./joaninhasCalorosas <<<"100 100 100 0 22.765 18 27 0.1 30 0.1 40 10000 2"
# &> jc-out.txt ; less out.txt
