#!/bin/bash

for typ in 1 2
    do
        for proc in 1 2
            do
                for MHL in 70 80 90 100
                    do
                        printf -v MMM "%03d" ${MHL}
                        screen -X -S TYPE${typ}PROC${proc}MHL${MMM} quit
                    done
            done
    done
