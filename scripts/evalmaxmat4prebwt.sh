#!/bin/sh

# set -e -x

USAGE="Usage: $0"

if test $# -ne 0
then
  echo ${USAGE}
  exit 1
fi


dbfiles=(EcoliK12 Afumigatus Scerevisiae Dmelanogaster_2L Homosapiens_21)
queryfiles=(EcoliO157H7 Anidulans Spombe Dpseudoobscura_2a3 Musmusculus_16)


export GT_MEM_BOOKKEEPING=on
export GT_ENV_OPTIONS=-spacepeak

code2file()
{
  case $1 in
    EcoliK12)
      echo "/home/jiabin/largesequences/EcoliK12.fa";;
    Afumigatus)
      echo "/home/jiabin/largesequences/Afumigatus.fa";;
    Scerevisiae)
      echo "/home/jiabin/largesequences/Scerevisiae.fa";;
    Dmelanogaster_2L)
      echo "/home/jiabin/largesequences/Dmelanogaster_2L.fa";;
    Homosapiens_21)
      echo "/home/jiabin/largesequences/Homosapiens_21.fa";;
    EcoliO157H7)
      echo "/home/jiabin/largesequences/EcoliO157H7.fa";;
    Anidulans)
      echo "/home/jiabin/largesequences/Anidulans.fa";;
    Spombe)
      echo "/home/jiabin/largesequences/Spombe.fa";;
    Dpseudoobscura_2a3)
      echo "/home/jiabin/largesequences/Dpseudoobscura_2a3.fa";;
    Musmusculus_16)
      echo "/home/jiabin/largesequences/Musmusculus_16.fa";;
    *)
      echo "$0: illegal filecode $1"
      exit 1;;
  esac
}

mkpackedindexandprebwt()
{
  dbfilename=$1
  dbfilepath=`code2file $dbfilename`
  shift  # so that one gets only the rest parameters from command line with $*
  printf "# RUN packedindex mkindex $dbfilename $*\n"
  ${RUNNER} gt packedindex mkindex -bsize 10 $* -dir rev -db ${dbfilepath} -indexname pck -sprank -dna -ssp -des -sds -pl | egrep '# TIME overall'
  printf "# RUN gt prebwt $dbfilename -maxdepth 12\n"
  ${RUNNER} gt prebwt -pck pck -maxdepth 12 | egrep '# TIME overall'
}

maxmat4()
{
  queryfilename=$2
  printf "# RUN maxmat4 $3 -p $4 $5 ($1/$queryfilename)\n"
  queryfilepath=`code2file $queryfilename`
  #printf "# gt dev maxmat4 $3 -p $4 $5 -L -showtime pck ${queryfilepath}\n"
  ${RUNNER} gt dev maxmat4 $3 -p $4 $5 -L -showtime pck ${queryfilepath} | egrep '# TIME overall|# space peak'
}

for rfc in $dbfiles
do
  fn=`code2file ${rfc}`
  if test ! -f ${fn}
  then
    echo "FAILURE: ${fn} does not exist"
    exit 1
  fi
done

for rfc in $queryfiles
do
  fn=`code2file ${rfc}`
  if test ! -f ${fn}
  then
    echo "FAILURE: ${fn} does not exist"
    exit 1
  fi
done

echo "# DATE `date +%Y-%m-%d-%H:%M`"

#matchoptions=( -mumreference -maxmatch -mum )
matchoptions=( -mumreference )


for (( i = 0 ; i < ${#dbfiles[@]} ; i++ ))
do
  mkpackedindexandprebwt ${dbfiles[$i]} -locfreq 8
  for matchoption in ${matchoptions[@]}
  do
    for minlength in 20
    do 
      maxmat4 ${dbfiles[$i]} ${queryfiles[$i]} ${matchoption} -l ${minlength}
    done
  done
  rm -f pck.*
done
