#!/bin/sh

#bin/gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev \
                       #-db testdata/at1MB -indexname at1MBpck \
                       #-sprank -dna -ssp -des -sds -pl
#env -i GT_MEM_BOOKKEEPING=off bin/gt dev maxmat4 -mumreference -b -l 25 \
#-L -s -c -showtime at1MBpck testdata/U89959_genomic.fas | egrep '# TIME overall' > tmp.prot


bin/gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db /home/jiabin/largesequences/EcoliK12.fa -indexname pck -sprank -dna -ssp -des -sds -pl

env -i GT_MEM_BOOKKEEPING=off bin/gt dev maxmat4 -maxmatch -l 20 -L pck /home/jiabin/largesequences/EcoliO157H7.fa | egrep '# TIME overall' > tmp.prot
gprof -p bin/gt gmon.out | less
