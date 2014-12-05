#!/bin/bash

rm joaninhasCalorosas ; make -k joaninhasCalorosas ; /usr/bin/time ./joaninhasCalorosas 10 10 12 0 1 1 1 0.0 1 0.0 1 1000 2 && cp ./jc-out.txt ./jc2-out.txt && /usr/bin/time ./joaninhasCalorosas 10 10 12 0 1 1 1 0.0 1 0.0 1 1000 1 && diff ./jc-out.txt ./jc2-out.txt
