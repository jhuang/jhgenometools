/*
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <stdbool.h>
#include <stdio.h>

//#include "core/array2dim_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"

#include "match/maxmat4-dfs.h"

/*
  if (currentdepth > 1UL)
  {
    outcol->prefixofsuffixbits
      = incol->prefixofsuffixbits &
        (mti->eqsvector[currentchar] >> (currentdepth-1));
  } else
  {
    outcol->prefixofsuffixbits = mti->eqsvector[currentchar];
  }
*/
int gt_pck_bitparallelism(const FMindex *index,
                            GT_UNUSED const GtEncseq *encseq,
                            unsigned long numofchars,
                            unsigned long totallength,
                            GtProgressTimer *timer,
                            GtLogger *logger,
                            GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtStackMaxmat4Node stack;
  Maxmat4Node root, current;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs;
  unsigned long resize = 64UL; /* XXX make this softcoded */
  unsigned long /*stackdepth,*/ maxdepth;
  
  unsigned long rangesize, idx, num_of_rows;
  unsigned int offset;
  Maxmat4Node child;

  rangeOccs = gt_calloc((size_t) GT_MULT2(numofchars), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (numofchars + 3), sizeof (*tmpmbtab ));
  //numoffiles = gt_encseq_num_of_files(encseq);
  //numofseq = gt_encseq_num_of_sequences(encseq);
  
  GT_STACK_INIT(&stack, resize);
  if (timer != NULL)
  {
    gt_progress_timer_start_new_state(timer,
                                      "obtain special pos",
                                      stdout);
  }
  
//  GT_STACK_NEXT_FREE(&stack,root);
  
//  root = gt_calloc(1, sizeof(Maxmat4Node));
//  GT_STACK_PUSH(&stack, *root);
  root.depth = 0;
  root.lower = 0;
  root.upper = totallength + 1;
  GT_STACK_PUSH(&stack,root);

  if (timer != NULL)
  {
    gt_progress_timer_start_new_state(timer,
                                      "traverse virtual tree",
                                      stdout);
  }
  //stackdepth = 1UL;
  maxdepth = 0;
  while (!GT_STACK_ISEMPTY(&stack))
  {
    //if (maxdepth < stackdepth)
    //  maxdepth = stackdepth;
    //current = stack.space + stack.nextfree - 1;

    //GT_STACK_DECREMENTTOP(&stack);
    //--stackdepth;
    current = GT_STACK_POP(&stack);

    gt_assert(current.lower < current.upper);
    num_of_rows = current.upper - current.lower;

    // 唯一的区别是 current.lower and current.upper 值不同
    rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                                  rangeOccs,
                                                  index,
                                                  current.lower,
                                                  current.upper);
    gt_assert(rangesize <= numofchars);

    printf("range size=%lu,current.lower=%lu,current.upper=%lu\n",rangesize,current.lower, current.upper);
    offset = 0;
    for (idx = 0; idx < rangesize; idx++)
    {
      gt_assert (tmpmbtab[idx].lowerbound <= tmpmbtab[idx].upperbound);
      gt_assert ((tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound) <=
                num_of_rows);
      num_of_rows -= (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound);
      printf("------tmpmbtab[%lu].lowerbound=%lu, tmpmbtab[idx].upperbound=%lu\n",idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound);

      if (tmpmbtab[idx].lowerbound != tmpmbtab[idx].upperbound)
      {
        if (tmpmbtab[idx].lowerbound + 1 ==
            tmpmbtab[idx].upperbound)
        { /* we found a leave on parent */
            printf("---leaf---tmpmbtab[idx].upperbound=%lu,current.depth=%lu\n",tmpmbtab[idx].upperbound,current.depth);
        } else
        {
            /* tmpmbtab[idx] is a branch of parent node */
            //GT_STACK_NEXT_FREE(&stack,child);
            
            child.lower = tmpmbtab[idx].lowerbound;
            child.upper = tmpmbtab[idx].upperbound;
            child.depth = current.depth + 1;
            GT_STACK_PUSH(&stack,child);
            printf("---push---child.lower=%lu, child.upper=%lu, child.depth=%lu\n",child.lower,child.upper,child.depth);
            //stackdepth++;
        }
      }
    }


  }
  gt_logger_log(logger, "max stack depth = %lu", maxdepth);
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  return had_err;
}
