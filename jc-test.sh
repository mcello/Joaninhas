#!/bin/bash

#rm joaninhasCalorosas ; make -k joaninhasCalorosas ; 
#time ./joaninhasCalorosas 10 10 12 0 1 1 1 0.0 1 0.0 1 1000 2 && cp ./jc-out.txt ./jc2-out.txt && 
#time ./joaninhasCalorosas 10 10 12 0 1 1 1 0.0 1 0.0 1 1000 1 && diff ./jc-out.txt ./jc2-out.txt

rm joaninhasCalorosas ; make -k joaninhasCalorosas ; 
/usr/bin/time ./joaninhasCalorosas 100 100 120 42 30 27 38 0.4 2 0.4 3 1000 1 && cp ./jc-out.txt ./jc1-out.txt && \
/usr/bin/time ./joaninhasCalorosas 100 100 120 42 30 27 38 0.4 2 0.4 3 1000 2 && cp ./jc-out.txt ./jc2-out.txt && \
/usr/bin/time ./joaninhasCalorosas 100 100 120 42 30 27 38 0.4 2 0.4 3 1000 4 && cp ./jc-out.txt ./jc4-out.txt && \
/usr/bin/time ./joaninhasCalorosas 100 100 120 42 30 27 38 0.4 2 0.4 3 1000 8 && cp ./jc-out.txt ./jc8-out.txt && \
diff ./jc1-out.txt ./jc8-out.txt

