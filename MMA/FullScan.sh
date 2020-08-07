#!/bin/bash

for typ in 1 2
    do
        for proc in 1 2
            do
                for MHL in 70 80 90 100
                    do
                        printf -v MMM "%03d" ${MHL}
                        screen -S TYPE${typ}PROC${proc}MHL${MMM} -d -m bash -c "./ScanCS.m $typ $proc $MHL; exec bash;"
                        sleep 5
                    done
            done
    done
