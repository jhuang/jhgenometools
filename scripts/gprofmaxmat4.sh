#!/bin/sh
#Yeast=`ls ${GTTESTDATA}/ltrharvest/s_cer/chr[01][0-9].*.gz`
#bin/gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev\ 
                       #-db ${Yeast}\
                       #-indexname Yeastpck\
                       #-sprank -dna -ssp -des -sds -pl

bin/gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev \
                       -db testdata/at1MB -indexname at1MBpck \
                       -sprank -dna -ssp -des -sds -pl
env -i GT_MEM_BOOKKEEPING=off bin/gt dev maxmat4 -mumreference -b -l 25 \
-L -s -c -showtime at1MBpck testdata/U89959_genomic.fas | egrep '# TIME overall' > tmp.prot
gprof bin/gt gmon.out
