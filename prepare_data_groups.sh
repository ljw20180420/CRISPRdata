#!/bin/bash

for group in spcas9 spymac ispymac
do
    case $group in
        spcas9)
            regex="^(A2-|A7-|D2-)"
        ;;
        spymac)
            regex="^(X-|x-|B2-|36t-)"
        ;;
        ispymac)
            regex="^(i10t-|i83-)"
        ;;
    esac
    mkdir -p /home/ljw/sdc1/SX/$group/algs
    for file in $(ls /home/ljw/sdc1/SX/algs | grep -E $regex)
    do
        if ! [ -e /home/ljw/sdc1/SX/$group/algs/$file ]
        then
            ln /home/ljw/sdc1/SX/algs/$file /home/ljw/sdc1/SX/$group/algs/$file
        fi
    done

    if ! [ -e /home/ljw/sdc1/SX/$group/$group.json ]
    then
        ./prepare_data.py --data_dir /home/ljw/sdc1/SX/$group/algs --name $group.json
    fi
done