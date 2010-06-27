#!/bin/sh

# set -e -x

if test $# -ne 1
then
  echo "Usage: $0 [small|all]"
  exit 1
fi

#repetitivefiles="random trnaglutamine"

case $1 in
  small) allfiles="at1MB atinsert duplicate randomsmall random"
         ;;
  all)   allfiles="at1MB atinsert duplicate randomsmall random random159 random160 randomn tttsmall trnaglutamine"
         ;;
  *)     allfiles=$1
         ;;
esac

code2file()
{
  case $1 in
    at1MB)
      echo "testdata/at1MB";;
    atinsert)
      echo "testdata/Atinsert.fna";;
    duplicate)
      echo "testdata/Duplicate.fna";;
    randomsmall)
      echo "testdata/Random-Small.fna";;
    random)
      echo "testdata/Random.fna";;
    random159)
      echo "testdata/Random159.fna";;
    random160)
      echo "testdata/Random160.fna";;
    randomn)
      echo "testdata/RandomN.fna";;
    tttsmall)
      echo "testdata/TTT-small.fna";;
    trnaglutamine)
      echo "testdata/trna_glutamine.fna";;
    *)
      echo "$0: illegal filecode $1"
      exit 1;;
  esac
}

#checkrepetitive()
#{
  #filename=$1                        # 如果接收的文件 在repetitivefiles表里, 错误.
  #for cfc in ${repetitivefiles}
  #do
    #if test ${cfc} == ${filename}
    #then
      #return 1
    #fi
  #done
  #return 0
#}

mkpackedindex()
{
  fc=$1
  filename=`code2file $1`
  shift
  printf "# RUN $fc $*\n"
  time gt packedindex mkindex -bsize 10 $* -dir rev -db ${filename} -indexname pck-idx -sprank -dna -ssp -des -sds -pl | egrep '# TIME overall|# space peak'
  #${RUNNER} gt suffixerator -v -showtime -indexname sfx-id -tis -suf -db ${filename} $* | egrep '# TIME overall|# space peak'
}

#mkesa()
#{
  #fc=$1
  #printf "# RUN $fc mkesa\n"
  #filename=`code2file $1`
  #runmkesa-sfx.sh ${filename}
#}

maxmat4()
{
  fc=$1
  printf "# RUN maxmat4 comparing with $fc\n"
  filename=`code2file $1`
  #mkesa -p mkesa-idx -b D -g suf -v -d $1
  #echo ${RUNNER}
  time gt dev maxmat4 $* -L pck-idx ${filename} | egrep '# TIME overall|# space peak'
  #time gt suffixerator -indexname pck-idx -dna -v -suf -tis -showtime -pl -dc 64 -db $1 
  #cmp -s sfx-idx.suf mkesa-idx.suf
  #runmkesa-sfx.sh ${filename}
}

for rfc in $allfiles
do
  fn=`code2file ${rfc}`
  if test ! -f ${fn}
  then
    echo "FAILURE: ${fn} does not exist"
    exit 1
  fi
done

# suffixerator ecoli2 -sat uint32 -dc 128
# exit 0

echo "# DATE `date +%Y-%m-%d-%H:%M`"
export GT_MEM_BOOKKEEPING=on
export GT_ENV_OPTIONS=-spacepeak
for rfc in $allfiles
do
  #checkrepetitive ${rfc}
  if test $? -eq 0
  then
    #echo "weg 1\n"
    mkpackedindex ${rfc} -locfreq 8
  fi
  for matchoption in "-mumreference -maxmatch -mum"
  do
    for minlength in 15 19 23
    do
      for rfc2 in $allfiles
      do 
        if test ${rfc2} != ${rfc}   
        then  
          #echo "weg 2\n"
          maxmat4 ${rfc2} ${matchoption} -l ${minlength}
        fi
      done    
    done
  done
  #mkesa ${rfc}
  rm -f pck-idx.* #mkesa-idx.*
done
