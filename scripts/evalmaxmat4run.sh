#!/bin/sh

# set -e -x

USAGE="Usage: $0 [small|all]"

if test $# -ne 1
then
  echo ${USAGE}
  exit 1
fi

case $1 in
  #small) allfiles="at1MB atinsert duplicate randomsmall random"
  small) allfiles="at1MB atinsert duplicate random"
         ;;
  all)   allfiles="at1MB atinsert duplicate randomsmall random random159 random160 randomn tttsmall trnaglutamine"
         ;;
  *)     allfiles=$1
         ;;
esac

export GT_MEM_BOOKKEEPING=on
export GT_ENV_OPTIONS=-spacepeak

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

mkpackedindex()
{
  dbfilename=$1
  dbfilepath=`code2file $dbfilename`
  shift  # so that one gets only the rest parameters from command line with $*
  printf "# RUN packedindex mkindex $dbfilename $*\n"
  ${RUNNER} gt packedindex mkindex -bsize 10 $* -dir rev -db ${dbfilepath} -indexname pck -sprank -dna -ssp -des -sds -pl | egrep '# TIME overall'
}

maxmat4()
{
  queryfilename=$2
  printf "# RUN maxmat4 $3 $4 $5 ($1/$queryfilename)\n"
  queryfilepath=`code2file $queryfilename`

  ${RUNNER} gt dev maxmat4 $3 $4 $5 -L -showtime pck ${queryfilepath} | egrep '# TIME overall|# space peak'
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

echo "# DATE `date +%Y-%m-%d-%H:%M`"

matchoptions=( -mumreference -maxmatch -mum )
#matchoptions=( -mumreference )

for rfc in $allfiles
do
  if test $? -eq 0
  then
    mkpackedindex ${rfc} -locfreq 8
  fi
  for matchoption in ${matchoptions[@]}
  do
    for minlength in 22
    do
      for rfc2 in $allfiles
      do 
        if test ${rfc2} != ${rfc}   
        then  
          maxmat4 ${rfc} ${rfc2} ${matchoption} -l ${minlength}
        fi
      done    
    done
  done
  rm -f pck.*
done
