#!/bin/bash

for a in $(seq 0 1 9); do
    mv waveform_$a.txt waveform_000$a.txt
done

for a in $(seq 10 1 99); do
    mv waveform_$a.txt waveform_00$a.txt
done

for a in $(seq 100 1 999); do
    mv waveform_$a.txt waveform_0$a.txt
done

