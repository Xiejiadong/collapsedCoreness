## collapsedCoreness

The code is for the Collapsed Coreness Problem, published in the paper "Discovering Key Users for Defending Network Structural Stability", Fan Zhang · Jiadong Xie · Kai Wang · Shiyu Yang · Yu Jiang, WWWJ 2021.

## Files

main.cc - source code for the paper

heap.cc / heap.h / defs.h / core.h / glist.cc / glist.h / treap.cc / treap.h - from the paper "A fast order-based approach for core maintenance. In: ICDE, pp. 337–348 (2017). DOI 10.1109/ICDE.2017.93"

data.txt - toy friendship data with 4039 vertices and 88234 edges

Makefile - compile commands

## compile and run

complie with g++ and -O3, use "make" in the file

run like "./core data budge", budge is a number

the program ouputs in data-budge-GCC.txt

## note

If you have any question, please contact me by xiejiadong0623@gmail.com.

If you used this code, please kindly cite the paper.