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

/** the function check if the mapped sequence is left maximal */
static bool isleftmaximal(const GtEncseq *encseq,
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
static unsigned long lcp(const GtEncseq *encseq,
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
  Maxmat4NodeBittab root, current, child;
  unsigned long resize = 64UL; 
  unsigned long rangesize, idx;
  unsigned long offset = 0UL;
  unsigned long i, j;   
	unsigned long realcalcquerypos=0, realcalcnodes=0, offsettimes=0;
	unsigned long 	moveunits;  
    
  GtUchar alphasize;
  unsigned int numofchars;
  /* GtAlphabet<=gt_encseq_alphabet, is different from MRAEnc alphabet */
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));  
  alphasize = (GtUchar) numofchars;  
  
  GtUchar cc;
  Symbol curSym;
  const MRAEnc *alphabet;
  struct GtUlongPair seqpospair;
  gt_assert(bwtSeq);
  const FMindex *fmindex = (const FMindex *)bwtSeq;
  alphabet = BWTSeqGetAlphabet(bwtSeq);
  const Mbtab *mbptr; 
  struct matchBound bwtbound;	  
  unsigned long *rangeOccs;
  rangeOccs = gt_calloc((size_t) GT_MULT2(alphasize), sizeof (*rangeOccs));
    
	GtBittab **eqsvector;
	GtBittab *tmpbittab = gt_bittab_new(bitlength);  /* temporary bittab for every round */	
	eqsvector = gt_malloc(sizeof (GtBittab*) * alphasize * leastlength);
	for (i=0; i<alphasize*leastlength; i++) {
		eqsvector[i] = gt_bittab_new(bitlength);
	}
	     
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
			gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlength,leastlength,query+offset,querylen-offset);
		} 
		else
		{
			gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlength,leastlength,query+offset,bitlength);  			
		}			
		 
		for (i=0; i<alphasize; i++) {
			gt_bittab_equal(tmpbittab, eqsvector[i*leastlength]); 
			for (j=1; j<leastlength; j++) {  /* j corresponds to current.depth */
			  gt_bittab_shift_right_equal(tmpbittab);
        gt_bittab_equal(eqsvector[i*leastlength+j], tmpbittab);
			}
	  }
	
		GT_STACK_INIT(&stack, resize);
		
		root.depth = 0;
		root.lower = 0;
		root.upper = totallength + 1;
		GtBittab *rootprefixofsuffixbits = gt_bittab_new(bitlength);
		for (i = 0; i < bitlength; i++) {
			gt_bittab_set_bit(rootprefixofsuffixbits, 1);
		}  /* init all positions as 1 */
		root.prefixofsuffixbits = rootprefixofsuffixbits;
		root.code = 0;
		GT_STACK_PUSH(&stack,root);

		if (showtime)
		{
			gt_progress_timer_start_new_state(timer,
																				"start to traverse tree",
																				stdout);
		}
    
		while (!GT_STACK_ISEMPTY(&stack))
		{
			current = GT_STACK_POP(&stack);
			gt_assert(current.lower < current.upper);
			
			/* find nodes whose depth is equal than least length as seeds */
			if (current.depth == leastlength)
			{			
				realcalcnodes++;
			  if (matchmode == GT_MATCHMODE_MAXMATCH) 
			  {
					/* save all subject positions in an array */
					const unsigned long subjectpositions_size = current.upper - current.lower;
					unsigned long subjectpositions[subjectpositions_size];

					unsigned long bwtboundthisline;
					for (bwtboundthisline=current.lower, i=0;bwtboundthisline < current.upper; bwtboundthisline++, i++)
					{					
							unsigned long subjectposthisline =
							gt_voidpackedfindfirstmatchconvert(fmindex,
																								 bwtboundthisline,
																								 current.depth);
							subjectpositions[i] = subjectposthisline;
					}					

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

							/* calculate match length of every combination between subjectpos and querypos */
							unsigned long j;
							for(j=0; j<subjectpositions_size; j++)
							{
								if ( isleftmaximal(encseq,subjectpositions[j],query,query+querypos) ) {
									unsigned long matchlength;
									if ((subjectpositions[j] + current.depth) < totallength)
									{
										unsigned long additionalmatchlength = 
																	lcp(encseq,
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
										gt_error_set(err, "The subject position plus current.depth exceeds totallength.");
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
							unsigned long depth = current.depth;

							while ( (querypos+depth) < querylen && (bwtbound.start+1) < bwtbound.end )
							{
								cc = *(query+querypos+depth);
								curSym = MRAEncMapSymbol(alphabet, cc);
								seqpospair = BWTSeqTransformedPosPairOcc(bwtSeq, curSym,
																												 bwtbound.start,bwtbound.end);
								bwtbound.start = bwtSeq->count[curSym] + seqpospair.a;
								bwtbound.end = bwtSeq->count[curSym] + seqpospair.b;		
								depth++;	
							}

							/* calculate match length of every combination between subjectpos and querypos */
							if (bwtbound.start+1 == bwtbound.end) {
								unsigned long subjectpos =
								gt_voidpackedfindfirstmatchconvert(fmindex,
																									 bwtbound.start,
																									 depth);
								if ( isleftmaximal(encseq,subjectpos,query,query+querypos) ) {
									unsigned long additionalmatchlength = 
																lcp(encseq,
																		subjectpos + depth,
																		totallength,
																		query + querypos + depth,
																		query + querylen);         /* qend */								
									unsigned long matchlength = depth + additionalmatchlength;	
									
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
				/* fill data in tmpmbtab from different current.lower and current.upper */
				rangesize
					= MRAEncGetRangeSize(EISGetAlphabet(bwtSeq->seqIdx),0);																							
				gt_assert(rangesize <= alphasize);									
				BWTSeqPosPairRangeOcc(bwtSeq, 0, current.lower, current.upper,rangeOccs);		
				
				for (idx = 0; idx < rangesize; idx++)
				{  
					if ((current.depth+1) <= maxdepth)
					{		   			   
						child.code = current.code * alphasize + idx;	
						mbptr = mbtab[current.depth+1] + child.code;	
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
						/* generate a lot of GtBittab objects for every node */
						GtBittab *prefixofsuffixbits = gt_bittab_new(bitlength); 
						if (current.depth > 0UL)
						{								
							gt_bittab_and(prefixofsuffixbits, current.prefixofsuffixbits, eqsvector[idx*leastlength+current.depth]);
						} else
						{
							gt_bittab_equal(prefixofsuffixbits, eqsvector[idx*leastlength]);   
						}

						if (gt_bittab_count_set_bits(prefixofsuffixbits) > 0) {		
							child.lower = bwtbound.start;  /* record match position in reference */
							child.upper = bwtbound.end;		
							child.depth = current.depth + 1;  			 /* record match length */
							child.prefixofsuffixbits = prefixofsuffixbits;  /*	record match position in query */
							/* child.code was added above in line 467 */	
							GT_STACK_PUSH(&stack,child);						 		
						} 
						else
						{
							gt_bittab_delete(prefixofsuffixbits);			
						}														
					}
				}
			}
			gt_bittab_delete(current.prefixofsuffixbits);		
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
		printf("# BIT PARALLELISM FACTOR=%3.2lf, nodes/move=%3.2lf, real calculated query positions=%lu, real calculated nodes=%lu\n",(double)realcalcquerypos/querylen, (double)realcalcnodes/(offsettimes+1), realcalcquerypos, realcalcnodes);
	}
	
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  
  gt_bittab_delete(tmpbittab);
  			
	for (i=0; i<alphasize*leastlength; i++) {
		gt_bittab_delete(eqsvector[i]);
	}  							
  gt_free(eqsvector);
  return had_err?-1:0;
}
