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

largesequencesfile_path="/home/jiabin"
code2file()
{
  case $1 in
    EcoliK12)
      echo "${largesequencesfile_path}/largesequences/EcoliK12.fa";;
    Afumigatus)
      echo "${largesequencesfile_path}/largesequences/wgs.AAHF.1.fsa_nt";;
    Scerevisiae)
      echo "${largesequencesfile_path}/largesequences/wgs.ACVY.1.fsa_nt";;
    Dmelanogaster_2L)
      echo "${largesequencesfile_path}/largesequences/Dmelanogaster_2L.fa";;
    Homosapiens_21)
      echo "${largesequencesfile_path}/largesequences/Homosapiens_21.fa";;
    EcoliO157H7)
      echo "${largesequencesfile_path}/largesequences/EcoliO157H7.fa";;
    Anidulans)
      echo "${largesequencesfile_path}/largesequences/wgs.AACD.1.fsa_nt";;
    Spombe)
      echo "${largesequencesfile_path}/largesequences/Spombe.fa";;
    Dpseudoobscura)
      echo "${largesequencesfile_path}/largesequences/wgs.AADE.1.fsa_nt";;
    Musmusculus_16)
      echo "${largesequencesfile_path}/largesequences/Musmusculus_16.fa";;
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
matchoptions=( -maxmatch )


for (( i = 0 ; i < ${#dbfiles[@]} ; i++ ))
do
  mkpackedindex ${dbfiles[$i]} -locfreq 8
  for matchoption in ${matchoptions[@]}
  do
    for minlength in 20
    do 
          maxmat4 ${dbfiles[$i]} ${queryfiles[$i]} ${matchoption} -l ${minlength}
    done
  done
  rm -f pck.*
done
