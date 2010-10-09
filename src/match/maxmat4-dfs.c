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

//#include "core/chardef.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/stack-inlined.h"
#include "core/unused_api.h"

#include "match/eis-voiditf.h"
#include "match/maxmat4-dfs.h"
#include "initeqsvec.h"

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
    /* gt_assert(dbrightbound < totallength); */
    if (dbrightbound < totallength) 
    {
			dbrightboundchar = gt_encseq_get_encoded_char(encseq,
																		 dbrightbound,
																		 GT_READMODE_FORWARD);
		}
  }
  return (unsigned long) (qptr-qnewstart);
}
                            
int gt_pck_bitparallelism(const GtUchar *query,
                          unsigned long querylen,
                          const FMindex *index,
                          const GtEncseq *encseq,
                          unsigned long totallength,
                          GT_UNUSED unsigned long leastlength,
                          //Findmatchfunction findmatchfunction,
                          GT_UNUSED const GtMatchmode matchmode,
                          GT_UNUSED Processmatchfunction processmatch,
                          GT_UNUSED Showspecinfo *showspecinfo,
                          bool showtime,
                          GtProgressTimer *timer,
                          GT_UNUSED GtLogger *logger,
                          GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtStackMaxmat4Node stack;
  Maxmat4Node root, current, child;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs;
  unsigned long resize = 64UL; 
  unsigned long /*maxdepth,*/ rangesize, idx, num_of_rows;
  unsigned long *eqsvector;
  
  unsigned long offset = 0UL;
  unsigned long bitlen = (unsigned long) (CHAR_BIT * sizeof (unsigned long));
  
  //char buffer1[GT_INTWORDSIZE+1], buffer2[GT_INTWORDSIZE+1];	
  
  GtUchar alphasize;
  unsigned int numofchars;
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  alphasize = (GtUchar) numofchars;
  
  rangeOccs = gt_calloc((size_t) GT_MULT2(alphasize), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (alphasize + 3), sizeof (*tmpmbtab ));
  eqsvector = gt_calloc(alphasize, sizeof (*eqsvector));
  
  while (offset < querylen) 
  {
		printf("------offset=%lu\n",offset);
		/* the position of most left query in boundary matches */
		unsigned long mostleftquerypos = offset + bitlen - (leastlength-1);
		/* in last round or querylen < 32 or 64 bits */
		bool inlastround = (querylen-offset <= bitlen);
		// 短于 32 bits 的第一轮 和 长于32 bits 的最后一轮
		if (inlastround)
		{
			gt_initeqsvectorrev(eqsvector,(unsigned long) alphasize,
										query+offset,querylen-offset);  /* why is rev-version used here? */
		} 
		else
		{
			gt_initeqsvectorrev(eqsvector,(unsigned long) alphasize,
										query+offset,bitlen);  			
		}			
		
		GT_STACK_INIT(&stack, resize);
		
		root.depth = 0;
		root.lower = 0;
		root.upper = totallength + 1;
		root.prefixofsuffixbits = ~0UL;
		GT_STACK_PUSH(&stack,root);

		if (showtime)
		{
			gt_progress_timer_start_new_state(timer,
																				"start to traverse tree",
																				stdout);
		}

		//maxdepth = 0;
		while (!GT_STACK_ISEMPTY(&stack))
		{
			current = GT_STACK_POP(&stack);
			//printf("---pop---current.lower=%lu, current.upper=%lu, current.depth=%lu\n",current.lower,current.upper,current.depth);
			//////if (maxdepth < current.depth)
			//////	maxdepth = current.depth;

			gt_assert(current.lower < current.upper);
			num_of_rows = current.upper - current.lower;

			// 唯一的区别是 current.lower and current.upper 值不同
			rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
																										rangeOccs,
																										index,
																										current.lower,
																										current.upper);
			gt_assert(rangesize <= alphasize);
			


			//printf("range size=%lu,current.lower=%lu,current.upper=%lu, alphazise=%u\n",rangesize,current.lower, current.upper,alphasize);
		 // printf("------cc=%u\n",cc);
			
			
			
			for (idx = 0; idx < rangesize; idx++)
			//for (idx = rangesize-1; idx >= 0; idx--)
			{    
				gt_assert (tmpmbtab[idx].lowerbound <= tmpmbtab[idx].upperbound);
				gt_assert ((tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound) <=
									num_of_rows);
				num_of_rows -= (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound);
				
				//printf("------tmpmbtab[%lu].lowerbound=%lu, tmpmbtab[idx].upperbound=%lu\n",idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound);


        /* only the nodes that has child is allowed to enter */
				if ( (tmpmbtab[idx].lowerbound != tmpmbtab[idx].upperbound) )  /* in reference with idx is extensible */
				{  

												/* we found a leave on parent */
												if ( (current.lower + 1 == current.upper) && (current.depth >= leastlength) ) {
									  //printf("-------------###normal###lowerbound=%lu,upperbound=%lu,leaf.depth=%lu--------------\n",current.lower,current.upper,current.depth);
													//printf("next(%s,%lu,depth=%lu)->%s\n",buffer1,idx,child.depth,buffer2);
													unsigned long subjectpos =\
													gt_voidpackedfindfirstmatchconvert(index,
																														 current.lower,
																														 current.depth);
																														 
													/* comparing the unique subject position with possible more than 1 query start position	*/									
													unsigned int i;
													GtBitsequence mask;
													for (i=0, mask = GT_FIRSTBIT;
															 i < (unsigned int) GT_INTWORDSIZE;
															 i++, mask >>= 1)
													{
														if (current.prefixofsuffixbits & mask) {
															unsigned long querypos;
															if (inlastround)
															{
																querypos = querylen - (GT_INTWORDSIZE-i); /* (GT_INTWORDSIZE-i) is length from match start to end of the snippt */ 
															} 
															else
															{
																querypos = offset + i;
															}
															//unsigned long querypos = querylen - (GT_INTWORDSIZE-i);
															if ( isleftmaximal(encseq,subjectpos,query,query+querypos) ) {                                   
																////printf("---leaf---lowerbound=%lu,upperbound=%lu,leaf.depth=%lu\n",current.lower,current.upper,current.depth);
																if (querypos+current.depth == offset+bitlen)
																{
																	if (querypos < mostleftquerypos) 
																	{
																		//printf("------querypos=%lu, mostleftquerypos=%lu\n",querypos, mostleftquerypos);
																		mostleftquerypos = querypos; 	
																	}	
																} 
																else
																{
                                    unsigned long additionalmatchlength = 
																									lcp(encseq,
																											subjectpos + current.depth,
																											totallength,
																											query+querypos+current.depth,
																											query+querylen);  /* qend */								
                                    unsigned long matchlength = current.depth + additionalmatchlength;
																		/*
																		 * print or save the result
																		 */
																		processmatch(encseq,
																								 query,
																								 querypos,  /* querypos */
																								 querylen,
																								 matchlength,
																								 subjectpos,
																								 showspecinfo);
																}
															}
														}
													}												
											}
					else
					{					
							/* tmpmbtab[idx] is a branch of current node, 
							 * that is, tmpmbtab[idx] can be the only child of current node
							 * or one of children of current onde */           
							if (current.depth > 0UL)
							{
								child.prefixofsuffixbits
									= current.prefixofsuffixbits &
										(eqsvector[(GtUchar)idx] << current.depth);
							} else
							{
								child.prefixofsuffixbits = eqsvector[(GtUchar)idx];
							}
							
							if (child.prefixofsuffixbits != 0) {		
								child.lower = tmpmbtab[idx].lowerbound;  /* record match position in reference */
							  child.upper = tmpmbtab[idx].upperbound;		
							  child.depth = current.depth + 1;  			 /* record match length */
								GT_STACK_PUSH(&stack,child);						 /*	record match position in query */		
							}
					}

											

						 

							
				}
				else
				{

				  //printf("-------------###normal###lowerbound=%lu,upperbound=%lu,leaf.depth=%lu--------------",current.lower,current.upper,current.depth);
				  //printf("------tmpmbtab[%lu].lowerbound=%lu, tmpmbtab[idx].upperbound=%lu\n",idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound);
				}
			}
			
		}
		
		if (inlastround)
		{
			offset = querylen;
		} 
		else
		{
			offset = mostleftquerypos;						
		}	
	}

  // gt_logger_log(logger, "max stack depth = %lu", maxdepth);
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_free(eqsvector);
  return had_err;
}
