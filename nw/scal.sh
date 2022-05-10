#!/bin/bash

get_hardware_info(){
    max_process=$(LANG=en_US.UTF-8; lscpu | grep "CPU(s)" | head -n 1 | awk '{print $2}' | bc)
    threads_per_core=$(LANG=en_US.UTF-8; lscpu | grep Thread | awk '{print $4}' | bc)
    max_cores=$((max_process/threads_per_core))
}

get_hardware_info

p=1
while (($p < $max_cores))
do
    export OMP_NUM_THREADS=$p; ./diag
    p=$((p*2))
done