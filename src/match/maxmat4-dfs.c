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
#include "maxmat4-initeqsvec.h"

#define BIT_LENGTH      1024  //NEW

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
  //unsigned long *eqsvector; NEW
  GtBittab **eqsvector;
  //GtBittab *eqsvector;
  //gt_bittab_unit_test(err);
  
  unsigned long offset = 0UL;
  unsigned long bitlen = BIT_LENGTH/*(unsigned long) (CHAR_BIT * sizeof (unsigned long))*/;  // NEW
  
  //char buffer1[GT_INTWORDSIZE+1], buffer2[GT_INTWORDSIZE+1];	
  
  GtUchar alphasize;
  unsigned int numofchars;
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  alphasize = (GtUchar) numofchars;
  
  rangeOccs = gt_calloc((size_t) GT_MULT2(alphasize), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (alphasize + 3), sizeof (*tmpmbtab ));
  // eqsvector = gt_calloc(alphasize, sizeof (*eqsvector));    NEW
  //unsigned long i, j, size, bit, counter;
  eqsvector = gt_malloc(sizeof (GtBittab*) * alphasize);
  //unsigned long size = gt_rand_max(BIT_LENGTH) + 1;
  unsigned long i;
  for (i = 0; i < alphasize; i++) {
    eqsvector[(GtUchar)i] = gt_bittab_new(BIT_LENGTH);
  }
  
  while (offset < querylen) 
  {
		//printf("------offset=%lu\n",offset);
		/* the position of most left query in boundary matches */
		unsigned long mostleftquerypos = offset + bitlen - (leastlength-1);
		/* in last round or querylen < 32 or 64 bits */
		bool inlastround = (querylen-offset <= bitlen);
		// 短于 32 bits 的第一轮 和 长于32 bits 的最后一轮  NEW
		if (inlastround)
		{  
			gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlen,query+offset,querylen-offset);  /* why is rev-version used here? */
		} 
		else
		{
			gt_maxmat4_initeqsvectorrev(eqsvector,alphasize,bitlen,query+offset,bitlen);  			
		}			
		
		GT_STACK_INIT(&stack, resize);
		
		root.depth = 0;
		root.lower = 0;
		root.upper = totallength + 1;
		GtBittab *rootprefixofsuffixbits = gt_bittab_new(BIT_LENGTH);
		long position;
		for (position = 0; position < BIT_LENGTH; position++) {
			gt_bittab_set_bit(rootprefixofsuffixbits, position);
		}  /* init all positions as 1 */
		root.prefixofsuffixbits = rootprefixofsuffixbits;
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
													unsigned long i;
													//GtBitsequence mask;
													//for (i=0, mask = GT_FIRSTBIT;
															 //i < (unsigned int) GT_INTWORDSIZE;
															 //i++, mask >>= 1)
													//{
													for (i=0; i<BIT_LENGTH; i++)
													{
														if (gt_bittab_bit_is_set(current.prefixofsuffixbits, i)) {
															unsigned long querypos;
															if (inlastround)
															{
																querypos = querylen - (BIT_LENGTH-i); /* (GT_INTWORDSIZE-i) is length from match start to end of the snippt */ 
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
							GtBittab *prefixofsuffixbits = gt_bittab_new(BIT_LENGTH);
              //////gt_bittab_show(current.prefixofsuffixbits, stdout);
							if (current.depth > 0UL)
							{
							  GtBittab *tmp = gt_bittab_new(BIT_LENGTH);
                gt_bittab_equal(tmp, eqsvector[(GtUchar)idx]);           
		            //memcpy(tmp, C[(GtUchar)idx], sizeof(C[(GtUchar)idx])); 
								unsigned long i;
								for (i=0;i < current.depth;i++) {
								  gt_bittab_shift_right_equal(tmp);
								}
								//gt_bittab_shift_right_equal_with_units(tmp, 10);  
								//////gt_bittab_show(tmp, stdout);
								gt_bittab_and(prefixofsuffixbits, current.prefixofsuffixbits, tmp);
								gt_bittab_delete(tmp);
							} else
							{
								gt_bittab_equal(prefixofsuffixbits, eqsvector[(GtUchar)idx]);   
							}
							//////printf("current.lower=%lu,current.upper=%lu,current.depth=%lu--%lu-->child.lower=%lu,child.upper=%lu,child.depth=%lu\n",current.lower,current.upper,current.depth,idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound,current.depth+1);
							//////gt_bittab_show(prefixofsuffixbits, stdout);
							//////printf("------------------------\n");
							if (gt_bittab_count_set_bits(prefixofsuffixbits) > 0) {		
								child.lower = tmpmbtab[idx].lowerbound;  /* record match position in reference */
							  child.upper = tmpmbtab[idx].upperbound;		
							  child.depth = current.depth + 1;  			 /* record match length */
							  child.prefixofsuffixbits = prefixofsuffixbits;  /*	record match position in query */
								GT_STACK_PUSH(&stack,child);						 		
							} 
							else
							{
								//printf("delete\n");
								gt_bittab_delete(prefixofsuffixbits);			
							}
								
					}

											

						 

							
				}
				else
				{

				  //printf("-------------###normal###lowerbound=%lu,upperbound=%lu,leaf.depth=%lu--------------",current.lower,current.upper,current.depth);
				  //printf("------tmpmbtab[%lu].lowerbound=%lu, tmpmbtab[idx].upperbound=%lu\n",idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound);
				}
			}
			gt_bittab_delete(current.prefixofsuffixbits);	
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
  //gt_free(eqsvector); NEW
  for (i = 0; i < alphasize; i++) {
    gt_bittab_delete(eqsvector[(GtUchar)i]);
  }
  gt_free(eqsvector);
  return had_err;
}