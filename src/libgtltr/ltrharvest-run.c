/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-maxpairs.pr"
#include "libgtmatch/intcode-def.h"
#include "libgtmatch/sfx-mappedstr.pr"
#include "libgtmatch/pos2seqnum.pr"

#include "ltrharvest-opt.h"
#include "repeats.h"
#include "searchforLTRs.h"
#include "duplicates.h"
#include "outputstd.h"
#include "outputfasta.h"
#include "outputgff3.h"

static int runltrharvest(LTRharvestoptions *lo, Env *env)
{
  Sequentialsuffixarrayreader *ssar; /* suffix array */
  Seqpos *markpos = NULL;
  unsigned long numofdbsequences;
  bool had_err = false;

  env_error_check(env);

  ssar = newSequentialsuffixarrayreaderfromfile(lo->str_indexname,
		                  SARR_LCPTAB | SARR_SUFTAB |
				  SARR_ESQTAB | SARR_DESTAB,
				  false,
				  env);
  if (ssar == NULL)
  {
    return -1;
  }

  /* test if motif is valid and encode motif */
  if (testmotifandencodemotif (&lo->motif,
                             alphabetSequentialsuffixarrayreader(ssar),
                             env) != 0)
  {
    had_err = true;
  }

  /* show defined option and values */
  if (!had_err && lo->verbosemode)
  {
    showuserdefinedoptionsandvalues(lo);
  }

  numofdbsequences = numofdbsequencesSequentialsuffixarrayreader(ssar);
  /* calculate markpos array for sequences offset */
  if ( !had_err && numofdbsequences > 1)
  {
    markpos = encseq2markpositions(
	encseqSequentialsuffixarrayreader(ssar),
	numofdbsequencesSequentialsuffixarrayreader(ssar),
	env);
    lo->markpos = markpos;
    if (markpos == NULL)
    {
      had_err = true;
    }
  }

  /* init array for maximal repeats */
  INITARRAY (&lo->repeatinfo.repeats, Repeat);
  lo->repeatinfo.ssarptr = ssar;

  /* search for maximal repeats */
  if (!had_err && enumeratemaxpairs(ssar,
		       getnumofcharsAlphabet(
		         alphabetSequentialsuffixarrayreader(ssar)),
		       encseqSequentialsuffixarrayreader(ssar),
		       readmodeSequentialsuffixarrayreader(ssar),
                       (unsigned int)lo->minseedlength,
		       (void*)simpleexactselfmatchstore,
		       lo,
		       NULL,
		       env) != 0)
  {
    had_err = true;
  }

  /* init array for candidate pairs */
  INITARRAY(&lo->arrayLTRboundaries, LTRboundaries);

  /* apply the filter algorithms */
  if (!had_err && searchforLTRs (ssar, lo, markpos, env) != 0)
  {
    had_err = true;
  }

  /* free array for maximal repeats */
  FREEARRAY(&lo->repeatinfo.repeats, Repeat);

  /* remove exact duplicates */
  if (!had_err)
  {
    removeduplicates(&lo->arrayLTRboundaries);
  }

  /* remove overlapping predictions if desired */
  if (!had_err && (lo->nooverlapallowed || lo->bestofoverlap))
  {
    removeoverlapswithlowersimilarity(&lo->arrayLTRboundaries,
                                      lo->nooverlapallowed);
  }

  /* print multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutput)
  {
    if (showpredictionsmultiplefasta(lo,
          markpos,
	  false,
	  60,
	  ssar,
	  true,
	  env) != 0)
    {
      had_err = true;
    }
  }

  /* print inner region multiple FASTA file of predictions */
  if (!had_err && lo->fastaoutputinnerregion)
  {
    if (showpredictionsmultiplefasta(lo,
          markpos,
	  true,
	  60,
	  ssar,
	  true,
	  env) != 0)
    {
      had_err = true;
    }
  }

  /* print GFF3 format file of predictions */
  if (!had_err && lo->gff3output)
  {
    if (printgff3format(lo,
	  ssar,
	  markpos,
	  env) != 0 )
    {
      had_err = true;
    }
  }

  if (!had_err && numofdbsequences > 1)
  {
    FREESPACE(markpos);
  }

  /* print predictions to stdout */
  if (!had_err)
  {
    showinfoiffoundfullLTRs(lo, ssar, env);
  }

  /* free prediction array */
  FREEARRAY(&lo->arrayLTRboundaries, LTRboundaries);
  /* free suffixarray */
  freeSequentialsuffixarrayreader(&ssar, env);

  return had_err ? -1 : 0;
}

int parseargsandcallltrharvest(int argc,const char *argv[],Env *env)
{
  LTRharvestoptions lo;
  int had_err;

  had_err = ltrharvestoptions(&lo,argc,argv,env);
  if (!had_err) {
    printargsline(argv,argc);
    had_err = runltrharvest(&lo,env);
  }
  wrapltrharvestoptions(&lo,env);

  return had_err;
}
