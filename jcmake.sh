#!/bin/bash

rm joaninhasCalorosas ; make -k joaninhasCalorosas ; ./joaninhasCalorosas &> out.txt ; less out.txt
