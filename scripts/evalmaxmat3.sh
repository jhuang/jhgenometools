#!/bin/sh

# set -e -x

USAGE="Usage: $0"

if test $# -ne 0
then
  echo ${USAGE}
  exit 1
fi


dbfiles=(EcoliK12 Afumigatus Scerevisiae Dmelanogaster_2L Homosapiens_21)
queryfiles=(EcoliO157H7 Anidulans Spombe Dpseudoobscura Musmusculus_16)


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
    Dpseudoobscura)
      echo "/home/jiabin/largesequences/Dpseudoobscura.fa";;
    Musmusculus_16)
      echo "/home/jiabin/largesequences/Musmusculus_16.fa";;
    *)
      echo "$0: illegal filecode $1"
      exit 1;;
  esac
}

maxmat3()
{
	dbfilename=$1
  dbfilepath=`code2file $dbfilename`
  queryfilename=$2
  printf "# RUN maxmat3 $3 $4 $5 ($dbfilename/$queryfilename)\n"
  queryfilepath=`code2file $queryfilename`
  ~/maxmat3/mm3src/maxmat3.x $3 $4 $5 -L -s -n ${dbfilepath} ${queryfilepath} | egrep '# TIME overall|# space peak'
  #${RUNNER} gt dev maxmat3 $3 $4 $5 -L -showtime pck 
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
  for matchoption in ${matchoptions[@]}
  do
    for minlength in 20
    do 
      maxmat3 ${dbfiles[$i]} ${queryfiles[$i]} ${matchoption} -l ${minlength}
    done
  done
  rm -f pck.*
done
