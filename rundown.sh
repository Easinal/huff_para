#!/bin/bash
make
numactl -i all ./huff_para > parallel.log
#taskset -c 0 ./huff_para > sequential.log

