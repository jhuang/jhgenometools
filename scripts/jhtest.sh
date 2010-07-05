#!/bin/sh
#
# Copyright (c) 2010 Jiabin Huang <jiabin.huang@studium.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#set -e -x

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
		./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_mumreference_benchmark' \
				 -gttestdata ${GTTESTDATA}
		./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_maxmatch_benchmark' \
				 -gttestdata ${GTTESTDATA}
		./testsuite.rb \
				 ${MC} -keywords 'generate_mummer_mum_benchmark' \
				 -gttestdata ${GTTESTDATA}				 
	else
		./testsuite.rb \
				 ${MC} ${SE} -keywords 'gt_maxmat4_inputerror'
		./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_mumreference_with_mummer' \
				 -gttestdata ${GTTESTDATA}
		./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_maxmatch_with_mummer' \
				 -gttestdata ${GTTESTDATA}
		./testsuite.rb \
				 ${MC} ${SE} -keywords 'check_maxmat4_mum_with_mummer' \
				 -gttestdata ${GTTESTDATA}				 
		./testsuite.rb \
         ${MC} ${SE} -keywords 'check_maxmat4_with_repfind' \
         -gttestdata ${GTTESTDATA}
	fi
  cd ..
fi
