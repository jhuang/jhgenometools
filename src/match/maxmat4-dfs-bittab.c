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

#include <stdio.h>

#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"
#include "match/maxmat4-dfs-bittab.h"
#include "maxmat4-initeqsvec.h"

#ifndef S_SPLINT_S
/** the function check if the mapped sequence is left maximal */
static bool bittab_isleftmaximal(const GtEncseq *encseq,
                         unsigned long subjectpos,
                         const GtUchar *query,
                         const GtUchar *qstart)
{
  bool leftmaximal = false;
  if (subjectpos==0 || qstart==query ) {
    leftmaximal = true;
  } else {
    GtUchar dbleftchar = gt_encseq_get_encoded_char(encseq,
                              subjectpos-1,
                              GT_READMODE_FORWARD);
    if (ISSPECIAL(dbleftchar) || dbleftchar != *(qstart-1) ) {
      leftmaximal = true;
    }  else {
      leftmaximal = false;
    }
  }
  return leftmaximal;
}

/**
 * the function calculate the long common prefix
 * between reference sequence and query sequence
 * from the qnewstart position.
 * qnewstart position refer to the first position >= leastlength
 * after general map-process with bwt
 */
static unsigned long bittab_lcp(const GtEncseq *encseq,
                         unsigned long dbrightbound,
                         unsigned long totallength,
                         const GtUchar *qnewstart,
                         const GtUchar *qend)
{
  const GtUchar *qptr = qnewstart;
  gt_assert(dbrightbound < totallength);
  GtUchar dbrightboundchar = gt_encseq_get_encoded_char(encseq,
                                           dbrightbound,
                                           GT_READMODE_FORWARD);

  while (qptr < qend &&
        (dbrightbound < totallength) && !ISSPECIAL(dbrightboundchar) &&
        *qptr == dbrightboundchar )
  {
    qptr++;
    dbrightbound++;
    if (dbrightbound < totallength)
    {
      dbrightboundchar = gt_encseq_get_encoded_char(encseq,
                                     dbrightbound,
                                     GT_READMODE_FORWARD);
    }
  }
  return (unsigned long) (qptr-qnewstart);
}

static int initialise_bittab_node(void *node)
{
  int had_err = 0;
  Maxmat4NodeBittab *tobeinitialised;
  tobeinitialised = (Maxmat4NodeBittab *) node;
  tobeinitialised->prefixofsuffixbits = NULL;
  return had_err;
}

int gt_pck_bitparallelism_bittab(const GtUchar *query,
                                  unsigned long querylen,
                                  const BWTSeq *bwtSeq,
                                  const GtEncseq *encseq,
                                  const Mbtab **mbtab,
                                  unsigned int maxdepth,
                                  unsigned long totallength,
                                  unsigned long leastlength,
                                  const GtMatchmode matchmode,
                                  Processmatchfunction processmatch,
                                  Showspecinfo *showspecinfo,
                                  unsigned long bitlength,
                                  bool showbitparallelismfactor,
                                  bool showtime,
                                  GtProgressTimer *timer,
                                  GT_UNUSED GtLogger *logger,
                                  GtError *err)
{
  bool had_err = false;
  GtStackMaxmat4NodeBittab stack;
  Maxmat4NodeBittab *root, current, *child;
  unsigned long resize = 64UL;
  unsigned long rangesize, idx;
  unsigned long offset = 0UL;
  unsigned long i, j;
  unsigned long realcalcquerypos=0, realcalcnodes=0, offsettimes=0;
  unsigned long moveunits;
  unsigned long subjectpositions_size, bwtnumber, subjectpos;
  unsigned long querypos, matchlength, additionalmatchlength, increaseddepth;
  unsigned long stackdepth, stackmaxdepth = 0UL, stackdepth_idx;

  GtUchar alphasize;
  unsigned int numofchars;
  /* GtAlphabet<=gt_encseq_alphabet, is different from MRAEnc alphabet */
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  alphasize = (GtUchar) numofchars;

  GtUchar cc;
  Symbol curSym;
  const MRAEnc *alphabet;
  GtUlongPair seqpospair;
  gt_assert(bwtSeq);
  const FMindex *fmindex = (const FMindex *)bwtSeq;
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  const Mbtab *mbptr;
  GtCodetype code = 0;
  struct matchBound bwtbound;
  unsigned long *rangeOccs;
  rangeOccs = gt_calloc((size_t) GT_MULT2(alphasize), sizeof (*rangeOccs));

  /* tmpcurrentbittab store the bittab of current node temporarily,
   * tmpbittab has 2 functions, first, it was used to assist to generate
   * eqsvector, secondly, it was used to store bittab of child node temporarily
   */
  GtBittab *tmpcurrentbittab =gt_bittab_new(bitlength),
            *tmpbittab = gt_bittab_new(bitlength);
  GtBittab **eqsvector;
  eqsvector = gt_malloc(sizeof (GtBittab*) * alphasize * leastlength);
  for (i=0; i<alphasize*leastlength; i++) {
    eqsvector[i] = gt_bittab_new(bitlength);
  }

  GT_STACK_INIT_WITH_INITFUNC(&stack, resize, initialise_bittab_node);
  while (offset < querylen)
  {
    /* printf("------offset=%lu\n",offset); */
    offsettimes++;
    /* initialize the position of most left query of boundary matches */
    moveunits = offset + bitlength - leastlength + 1;
    /* in last round or querylen < bitlength */
    bool islastround = (querylen-offset <= bitlength);

    if (islastround)
    {
      gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlength,\
                        leastlength,query+offset,querylen-offset);
    }
    else
    {
      gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlength,\
                        leastlength,query+offset,bitlength);
    }

    for (i=0; i<alphasize; i++) {
      gt_bittab_equal(tmpbittab, eqsvector[i*leastlength]);
      for (j=1; j<leastlength; j++) {  /* j corresponds to current.depth */
        gt_bittab_shift_right_equal(tmpbittab);
        gt_bittab_equal(eqsvector[i*leastlength+j], tmpbittab);
      }
    }

    GT_STACK_NEXT_FREE(&stack,root);
    root->depth = 0;
    root->lower = 0;
    root->upper = totallength + 1;
    /* malloc the space if next free was not allocated before */
    if (root->prefixofsuffixbits == NULL)
    {
      root->prefixofsuffixbits  = gt_bittab_new(bitlength);
    }
    for (i = 0; i < bitlength; i++) {
      gt_bittab_set_bit(root->prefixofsuffixbits, i);
    }  /* init all positions as 1 */
    root->code = 0;
    stackdepth = 1UL;

    if (showtime)
    {
      gt_progress_timer_start_new_state(timer,
                                        "start to traverse tree",
                                        stdout);
    }

    while (!GT_STACK_ISEMPTY(&stack))
    {
      if (stackmaxdepth < stackdepth)
      {
        stackmaxdepth = stackdepth;
      }

      current = GT_STACK_POP(&stack);
      stackdepth--;
      gt_assert(current.lower < current.upper);

      /* find nodes whose depth is equal than least length as seeds */
      if (current.depth == leastlength)
      {

        realcalcnodes++;
        if (matchmode == GT_MATCHMODE_MAXMATCH)
        {
          /* save all subject positions in an array */
          subjectpositions_size = current.upper - current.lower;
          unsigned long subjectpositions[subjectpositions_size];

          for (bwtnumber=current.lower, i=0;\
               bwtnumber < current.upper; bwtnumber++, i++)
          {
              subjectpos =
              gt_voidpackedfindfirstmatchconvert(fmindex,
                                                 bwtnumber,
                                                 current.depth);
              subjectpositions[i] = subjectpos;
          }

          for (i=0; i<bitlength; i++)
          {
            /* for every query position */
            if (gt_bittab_bit_is_set(current.prefixofsuffixbits, i))
            {
              realcalcquerypos++;
              if (islastround)
              {
                querypos = querylen - (bitlength-i);
              }
              else
              {
                querypos = offset + i;
              }

              /* calculate match length of every combination
               * between subjectpos and querypos */
              for (j=0; j<subjectpositions_size; j++)
              {
                if ( bittab_isleftmaximal(encseq,subjectpositions[j],\
                                   query,query+querypos) ) {
                  if ((subjectpositions[j] + current.depth) < totallength)
                  {
                    additionalmatchlength =
                                  bittab_lcp(encseq,
                                      subjectpositions[j] + current.depth,
                                      totallength,
                                      query + querypos + current.depth,
                                      query+querylen);         /* qend */
                    matchlength = current.depth + additionalmatchlength;
                  }
                  else if ((subjectpositions[j] + current.depth) == totallength)
                  {
                    matchlength = current.depth;
                  }
                  else
                  {
                    gt_error_set(err, \
                    "The subject position plus \
                    current.depth exceeds totallength.");
                    had_err = true;
                    break;
                  }

                  /*
                   * print or save the result
                   */
                  processmatch(encseq,
                               query,
                               querypos,
                               querylen,
                               matchlength,
                               subjectpositions[j],
                               showspecinfo);
                }
              }
            }
          }
        }
        else  /* MUM or MUMREFERENCE */
        {
          for (i=0; i<bitlength; i++)
          {
            /* for every query position */
            if (gt_bittab_bit_is_set(current.prefixofsuffixbits, i))
            {
              realcalcquerypos++;
              unsigned long querypos;
              if (islastround)
              {
                querypos = querylen - (bitlength-i);
              }
              else
              {
                querypos = offset + i;
              }

              bwtbound.start = current.lower;
              bwtbound.end = current.upper;
              increaseddepth = current.depth;

              while ( (querypos+increaseddepth) < querylen
                        && (bwtbound.start+1) < bwtbound.end )
              {
                cc = *(query+querypos+increaseddepth);
                curSym = MRAEncMapSymbol(alphabet, cc);
                seqpospair = BWTSeqTransformedPosPairOcc(bwtSeq, curSym,
                                            bwtbound.start,bwtbound.end);
                bwtbound.start = bwtSeq->count[curSym] + seqpospair.a;
                bwtbound.end = bwtSeq->count[curSym] + seqpospair.b;
                increaseddepth++;
              }

              /* calculate match length of every combination
               * between subjectpos and querypos */
              if (bwtbound.start+1 == bwtbound.end) {
                subjectpos =
                gt_voidpackedfindfirstmatchconvert(fmindex,
                                                   bwtbound.start,
                                                   increaseddepth);
                if ( bittab_isleftmaximal(encseq,
                                          subjectpos,query,query+querypos) ) {
                  additionalmatchlength =
                                bittab_lcp(encseq,
                                    subjectpos + increaseddepth,
                                    totallength,
                                    query + querypos + increaseddepth,
                                    query + querylen);         /* qend */
                  matchlength = increaseddepth + additionalmatchlength;

                  /*
                   * print or save the result
                   */
                  processmatch(encseq,
                               query,
                               querypos,
                               querylen,
                               matchlength,
                               subjectpos,
                               showspecinfo);
                }
              }
            }
          }
        }
      }
      else if (current.depth < leastlength)
      {
        /* fill data in tmpmbtab from different
         * current.lower and current.upper */
        rangesize
          = MRAEncGetRangeSize(EISGetAlphabet(bwtSeq->seqIdx),0);
        gt_assert(rangesize <= alphasize);
        BWTSeqPosPairRangeOcc(bwtSeq, 0, current.lower,
                                         current.upper,rangeOccs);

        gt_bittab_equal(tmpcurrentbittab, current.prefixofsuffixbits);
        for (idx = 0; idx < rangesize; idx++)
        {
          if ((current.depth+1) <= maxdepth)
          {
            /* use bwt bounds from prebwt file */
            code = current.code * alphasize + idx;
            mbptr = mbtab[current.depth+1] + code;
            bwtbound.start = mbptr->lowerbound;
            bwtbound.end = mbptr->upperbound;
          }
          else
          {
            bwtbound.start = bwtSeq->count[idx] + rangeOccs[idx];
            bwtbound.end = bwtSeq->count[idx] + rangeOccs[rangesize+idx];
          }

          /* in reference with idx is extensible */
          if (bwtbound.start < bwtbound.end)
          {
            if (current.depth > 0UL)
            {
              gt_bittab_and(tmpbittab, tmpcurrentbittab,
                                  eqsvector[idx*leastlength+current.depth]);
            } else
            {
              gt_bittab_equal(tmpbittab,
                                  eqsvector[idx*leastlength]);
            }

            if (gt_bittab_count_set_bits(tmpbittab) > 0) {
               GT_STACK_NEXT_FREE(&stack,child);
              /* malloc the space if next free was not allocated before */
              if (child->prefixofsuffixbits == NULL)
              {
                child->prefixofsuffixbits = gt_bittab_new(bitlength);
              }
              gt_bittab_equal(child->prefixofsuffixbits, tmpbittab);

              /* record match position in reference */
              child->lower = bwtbound.start;
              child->upper = bwtbound.end;
              child->depth = current.depth + 1;  /* record match length */
              child->code = code;
              stackdepth++;
            }
          }
        }
      }
    }

    if (islastround)
    {
      offset = querylen;
    }
    else
    {
      offset = moveunits;
    }
  }

  if (showbitparallelismfactor && querylen!=0)
  {
    printf("# BIT PARALLELISM FACTOR=%3.2lf, nodes/move=%3.2lf, \
    real calculated query positions=%lu, real calculated nodes=%lu\n", \
    (double)realcalcquerypos/querylen, (double)realcalcnodes/(offsettimes+1), \
    realcalcquerypos, realcalcnodes);
  }

  for (stackdepth_idx = 0; stackdepth_idx < stackmaxdepth; stackdepth_idx++)
  {
    gt_bittab_delete(stack.space[stackdepth_idx].prefixofsuffixbits);
  }
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);

  gt_bittab_delete(tmpbittab);
  gt_bittab_delete(tmpcurrentbittab);

  for (i=0; i<alphasize*leastlength; i++) {
    gt_bittab_delete(eqsvector[i]);
  }
  gt_free(eqsvector);
  return had_err?-1:0;
}
#endif
