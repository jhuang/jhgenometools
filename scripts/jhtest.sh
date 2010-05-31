#!/bin/sh

USAGE="Usage: $0 [-memcheck]"

rungenerate=0
if test $# -eq 0
then
  MC=""
  SE=""
else
  if test $# -eq 1
  then
    if test "$1" = "-memcheck"
    then
      MC="-memcheck" 
    else 
      if test "$1" = "-g"
      then 
        rungenerate=1
      else
        echo ${USAGE}
        exit 1
      fi
    fi 
  else 
    if test $# -eq 2
    then
      if test "$1" = "-select"
      then
        SE="-select $2"      #TODO: check if $2 is integer        
      else
        echo ${USAGE}
        exit 1
      fi
    fi  
  fi
fi

cerr()
{
  $*
  if [ $? -ne 0 ]
  then
    echo "failure: $*"
    exit 1
  fi
}

runtestsuite=1
if test $runtestsuite -eq 1
then
  cd testsuite
  if test $rungenerate -eq 1
  then
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_mumreference_benchmark' \
				 -gttestdata ${GTTESTDATA}
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_maxmatch_benchmark' \
				 -gttestdata ${GTTESTDATA}
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_mum_benchmark' \
				 -gttestdata ${GTTESTDATA}				 
	else
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_mumreference_with_mummer' \
				 -gttestdata ${GTTESTDATA}
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_maxmatch_with_mummer' \
				 -gttestdata ${GTTESTDATA}
		env -i GT_MEM_BOOKKEEPING=on ./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_mum_with_mummer' \
				 -gttestdata ${GTTESTDATA}				 
		env -i GT_MEM_BOOKKEEPING=on GTTESTDATA=${HOME}/gttestdata ./testsuite.rb \
       ${MC} -keywords 'check_maxmat4_with_repfind' \
       -gttestdata ${GTTESTDATA}
	fi
  cd ..
fi
