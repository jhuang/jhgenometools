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
                            
int gt_pck_bitparallelism(const GtUchar *query,
                          unsigned long querylen,
                          const FMindex *index,
                          const GtEncseq *encseq,
                          unsigned long totallength,
                          GT_UNUSED unsigned long leastlength,
                          //Findmatchfunction findmatchfunction,
                          const GtMatchmode matchmode,
                          GT_UNUSED Processmatchfunction processmatch,
                          GT_UNUSED Showspecinfo *showspecinfo,
                          bool showtime,
                          GtProgressTimer *timer,
                          GtLogger *logger,
                          GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtStackMaxmat4Node stack;
  Maxmat4Node root, current, child;
  Mbtab *tmpmbtab;
  unsigned long *rangeOccs;
  unsigned long resize = 64UL; 
  unsigned long maxdepth, rangesize, idx, num_of_rows;
  unsigned long *eqsvector;
  
  char buffer1[GT_INTWORDSIZE+1], buffer2[GT_INTWORDSIZE+1];	
  
  GtUchar alphasize;
  unsigned int numofchars;
  numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  alphasize = (GtUchar) numofchars;
  
  rangeOccs = gt_calloc((size_t) GT_MULT2(alphasize), sizeof (*rangeOccs));
  tmpmbtab = gt_calloc((size_t) (alphasize + 3), sizeof (*tmpmbtab ));
  eqsvector = gt_calloc(alphasize, sizeof (*eqsvector));
  gt_initeqsvectorrev(eqsvector,(unsigned long) alphasize,
                query,querylen);  /* why is rev-version used here? */
  
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

  maxdepth = 0;
  while (!GT_STACK_ISEMPTY(&stack))
  {
    current = GT_STACK_POP(&stack);
    //printf("---pop---current.lower=%lu, current.upper=%lu, current.depth=%lu\n",current.lower,current.upper,current.depth);
    if (maxdepth < current.depth)
      maxdepth = current.depth;

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
    
		/* node current is now at the end of reference */
		bool currentatend = ( (current.prefixofsuffixbits != 0) && (totallength==(gt_voidpackedfindfirstmatchconvert(index,current.lower,current.depth)+current.depth)) );
		if (currentatend) {
			////printf("-------------###currentatend###lowerbound=%lu,upperbound=%lu,leaf.depth=%lu------------------\n",current.lower,current.upper,current.depth);
							/* process for pop is following */
							if (matchmode == GT_MATCHMODE_MAXMATCH) 
							{
								/* for option maxmatch */								
								if ( current.depth>=leastlength ) {
									                                  
                  /* 三次turn, 每次6个, 中间那一次全军覆没，因为那一次的 child.prefixofsuffixbits ！＝0， 有子 */
									unsigned long bwtboundthisline;
									for (bwtboundthisline=current.lower;bwtboundthisline < current.upper; bwtboundthisline++)
									{	
										////printf("bwtboundthisline=%lu\n",bwtboundthisline );
										if ( SEPARATOR == gt_bwtseqgetsymbol(bwtboundthisline, index) )	
										{					
											unsigned long subjectposthisline =
											gt_voidpackedfindfirstmatchconvert(index,
																												 bwtboundthisline,
																												 current.depth);
																											 
											/* comparing the unique subject position with possible more than 1 query start position */										
											unsigned int i;
											GtBitsequence mask;
											for (i=0, mask = GT_FIRSTBIT;
													 i < (unsigned int) GT_INTWORDSIZE;
													 i++, mask >>= 1)
											{
												if (current.prefixofsuffixbits & mask) {
													unsigned long querypos = querylen - (GT_INTWORDSIZE-i);
													if ( isleftmaximal(encseq,subjectposthisline,query,query+querypos) ) {                                   
														////printf("---leaf---lowerbound=%lu,upperbound=%lu,leaf.depth=%lu\n",current.lower,current.upper,current.depth);

														/*
														 * print the result
														 */
														processmatch(encseq,
																				 query,
																				 querypos,  /* querypos */
																				 querylen,
																				 current.depth,
																				 subjectposthisline,
																				 showspecinfo);
													}
												}
											}
									  }											
									}												
								}
							} 
							else
							{
								/* for option mumreference or mum */
							  if ( (current.depth >= leastlength) && (current.lower + 1 == current.upper) && (SEPARATOR==gt_bwtseqgetsymbol(current.lower,index))	) {

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
											unsigned long querypos = querylen - (GT_INTWORDSIZE-i);
											if ( isleftmaximal(encseq,subjectpos,query,query+querypos) ) {                                   
												////printf("---leaf---lowerbound=%lu,upperbound=%lu,leaf.depth=%lu\n",current.lower,current.upper,current.depth);
												/*
												 * print or save the result
												 */
												processmatch(encseq,
																		 query,
																		 querypos,  /* querypos */
																		 querylen,
																		 current.depth,
																		 subjectpos,
																		 showspecinfo);
											}
										}
									}												
								}
						  }
		}
    for (idx = 0; idx < rangesize; idx++)
    {    
      gt_assert (tmpmbtab[idx].lowerbound <= tmpmbtab[idx].upperbound);
      gt_assert ((tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound) <=
                num_of_rows);
      num_of_rows -= (tmpmbtab[idx].upperbound - tmpmbtab[idx].lowerbound);
      
      //printf("------tmpmbtab[%lu].lowerbound=%lu, tmpmbtab[idx].upperbound=%lu\n",idx,tmpmbtab[idx].lowerbound,tmpmbtab[idx].upperbound);


      if ( (tmpmbtab[idx].lowerbound != tmpmbtab[idx].upperbound) )  /* in reference with idx is extensible */
      {      
				    // 肯定有孩子，关于孩子是不是正确的， 就要看 bitpara 了。比如 atgct是他的孩子 
				    /* tmpmbtab[idx] is a branch of parent node */           
            child.lower = tmpmbtab[idx].lowerbound;
            child.upper = tmpmbtab[idx].upperbound;
            child.depth = current.depth + 1;  
            child.idx = idx;
						if (child.depth > 1UL)
						{
							child.prefixofsuffixbits
								= current.prefixofsuffixbits &
									(eqsvector[(GtUchar)idx] << (child.depth-1));
						} else
						{
							child.prefixofsuffixbits = eqsvector[(GtUchar)idx];
						}
										
						gt_bitsequence_tostring(buffer1,(GtBitsequence) current.prefixofsuffixbits);
						gt_bitsequence_tostring(buffer2,(GtBitsequence) child.prefixofsuffixbits);
						//printf("next(%s,%lu,depth=%lu)->%s\n",buffer1,idx,child.depth,buffer2);
					 
						/* if no corresponding match finds in query or it reaches the end of reference, print the results out */
						if ( (child.prefixofsuffixbits != 0) && (totallength==(gt_voidpackedfindfirstmatchconvert(index,child.lower,child.depth)+child.depth)) ) {
							
							GT_STACK_PUSH(&stack,child);		
							//printf("---push1---child.lower=%lu, child.upper=%lu, child.depth=%lu\n",child.lower,child.upper,child.depth);	
								
            } else if (child.prefixofsuffixbits == 0) {  /* in query with idx is not extensible */
							/* process for pop is following */
							
							////printf("-------------###normal###lowerbound=%lu,upperbound=%lu,leaf.depth=%lu--------------\n",current.lower,current.upper,current.depth);
							if (matchmode == GT_MATCHMODE_MAXMATCH) 
							{
								/* for option maxmatch */								
								if ( current.depth>=leastlength ) {
									////printf("next(%s,%lu,depth=%lu)->%s\n",buffer1,idx,child.depth,buffer2);
									                              
                  /* 三次turn, 每次6个, 中间那一次全军覆没，因为那一次的 child.prefixofsuffixbits ！＝0， 有子 */
									unsigned long bwtboundthisline;
									for (bwtboundthisline=current.lower;bwtboundthisline < current.upper; bwtboundthisline++)
									{	
										////printf("bwtboundthisline=%lu\n",bwtboundthisline );
										if ( child.idx == gt_bwtseqgetsymbol(bwtboundthisline, index) && child.prefixofsuffixbits == 0 )	
										{					
											unsigned long subjectposthisline =
											gt_voidpackedfindfirstmatchconvert(index,
																												 bwtboundthisline,
																												 current.depth);
																											 
											/* comparing the unique subject position with possible more than 1 query start position */								
											unsigned int i;
											GtBitsequence mask;
											for (i=0, mask = GT_FIRSTBIT;
													 i < (unsigned int) GT_INTWORDSIZE;
													 i++, mask >>= 1)
											{
												if (current.prefixofsuffixbits & mask) {
													unsigned long querypos = querylen - (GT_INTWORDSIZE-i);
													if ( isleftmaximal(encseq,subjectposthisline,query,query+querypos) ) {                                   
														////printf("---leaf---lowerbound=%lu,upperbound=%lu,leaf.depth=%lu\n",current.lower,current.upper,current.depth);
														/*
														 * print the result
														 */
														processmatch(encseq,
																				 query,
																				 querypos,  /* querypos */
																				 querylen,
																				 current.depth,
																				 subjectposthisline,
																				 showspecinfo);
													}
												}
											}
									  }											
									}											
								}
							} 
							else
							{
								/* for option mumreference or mum */
								if ( (current.depth >= leastlength) && (current.lower + 1 == current.upper) ) {

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
											unsigned long querypos = querylen - (GT_INTWORDSIZE-i);
											if ( isleftmaximal(encseq,subjectpos,query,query+querypos) ) {                                   
												////printf("---leaf---lowerbound=%lu,upperbound=%lu,leaf.depth=%lu\n",current.lower,current.upper,current.depth);
												/*
												 * print or save the result
												 */
												processmatch(encseq,
																		 query,
																		 querypos,  /* querypos */
																		 querylen,
																		 current.depth,
																		 subjectpos,
																		 showspecinfo);
											}
										}
									}												
								}
						  }
						} else {
							GT_STACK_PUSH(&stack,child);
              //printf("---push---child.lower=%lu, child.upper=%lu, child.depth=%lu\n",child.lower,child.upper,child.depth);
						}
						
      }
    }


  }
  //printf("---maxdepth---maxdepth=%lu\n",maxdepth);
  //    gt_alphabet_decode_seq_to_fp(gt_encseq_alphabet(encseq),stdout,
  //                               query,4);
  //  (void) putchar('\n');
  gt_logger_log(logger, "max stack depth = %lu", maxdepth);
  GT_STACK_DELETE(&stack);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  gt_free(eqsvector);
  return had_err;
}
