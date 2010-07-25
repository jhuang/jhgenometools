#!/bin/sh

bin/gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db /home/jiabin/largesequences/EcoliK12.fa -indexname pck -sprank -dna -ssp -des -sds -pl
bin/gt prebwt -pck pck -maxdepth 12

env -i GT_MEM_BOOKKEEPING=off bin/gt dev maxmat4 -maxmatch -p -l 20 -L pck /home/jiabin/largesequences/EcoliO157H7.fa | egrep '# TIME overall' > tmp.prot
gprof -p bin/gt gmon.out | less
