#!/usr/bin/env bash

logfile="perf_all_manual_improvements.log"
rm -f "$logfile"

for opt in -O1 -O2 -O3 -Ofast; do
  for arch in "" "-xhost" "-xCORE-AVX512"; do
    for ipo in "" "-ipo"; do
      for alias in "" "-fno-alias"; do
        if [ "$arch" = "-xCORE-AVX512" ]; then
            for avx in "" "-qopt-zmm-usage=high "; do

            flags="$opt $arch $ipo $alias $avx"
            echo "=== Flags: $flags ===" >> "$logfile"

            make clean 
            make CFLAGS="$flags"

            for i in {1..2}; do
            echo "Durchlauf $i" >> "$logfile"
            { time ./heat jacobi.dat; time ./heat gauss_seidel.dat;} 2>>"$logfile" >>"$logfile"
            done

            printf '\n%.0s' {1..6} >> "$logfile"
            done
        else
            flags="$opt $arch $ipo $alias"
            echo "=== Flags: $flags ===" >> "$logfile"

            make clean 
            make CFLAGS="$flags"

            for i in {1..2}; do
            echo "Durchlauf $i" >> "$logfile"
            { time ./heat jacobi.dat; time ./heat gauss_seidel.dat;} 2>>"$logfile" >>"$logfile"
            done

            printf '\n%.0s' {1..6} >> "$logfile"
        fi


      done
    done
  done
done

