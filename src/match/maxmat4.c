/*
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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


//#include <inttypes.h>
//#include <string.h>
//#include <stdbool.h>
//#include "core/alphabet.h"
//#include "core/error.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/unused_api.h"
//#include "core/defined-types.h"
#include "spacedef.h"
//#include "optionargmode.h"
#include "core/format64.h"
//#include "maxmat4.h"
//#include "intcode-def.h"
//#include "core/encodedsequence.h"
//#include "initbasepower.h"
#include "match/eis-bwtseq.h"
//#include "match/esa-mmsearch.c"
#include "extended/reverse.h"


typedef void (*Preprocessmatchfunction)(uint64_t,
                                       unsigned long querylen,
                                       const char *,
                                       void *);
typedef bool (*Processmatchfunction)(const BWTSeq *,
                                const GtEncseq *,
                                unsigned long,
                                unsigned long,
                                       //unsigned long *,  
                                       const GtUchar *,   
                                       const GtUchar *,   
                                       const GtUchar *,    
                                       GtArray *);
                                
typedef void (*Postprocessmatchfunction)(void *,
                                      const GtUchar *,
                                      unsigned long);

typedef void (*Showmatchfunction)(const GtEncseq *encseq,
                                const GtAlphabet *,
                                const GtUchar *,
                                unsigned long,
                                unsigned long,
                                unsigned long,
                                unsigned long,
                                void *);
                                
   

typedef struct
{
  bool nucleotidesonly,                                      
       showstring,                           
       showreversepositions,                     
       showsequencelengths;   
  GtReadmode queryreadmode;
} Showspecinfo;
                                
typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  GtMatchmode matchmode;
  const GtAlphabet *alphabet;
  Preprocessmatchfunction preprocessmatchfunction;
  Processmatchfunction processmatchfunction;
  Postprocessmatchfunction postprocessmatchfunction;
  Showmatchfunction showmatchfunction;
  //void *processinfo;
  const GtEncseq *encseq;
  
  GtArray *mumcandtab;
  GtArray *maximalmatchtab;  
  Definedunsignedlong leastlength;
  Showspecinfo *showspecinfo;
} Matchprocessinfo;  


                                                                           
                                        
/*
	The following function compares two MUM-candidates. The MUM-candidate
	with smaller dbstart-value comes first.
	If both MUMs have the same dbstart-value, then the MUM-candidate
	with the larger length comes first.
*/
static int compareMUMcandidates(MUMcandidate *p,MUMcandidate *q)
{
	if(p->dbstart == q->dbstart)
	{
		return (p->mumlength < q->mumlength) ? 1 : -1;
	}
	return (p->dbstart > q->dbstart) ? 1 : -1;
}
//int gt_range_compare(const GtRange *range_a, const GtRange *range_b)
//{
  //gt_assert(range_a->start <= range_a->end && range_b->start <= range_b->end);

  //if ((range_a->start == range_b->start) && (range_a->end == range_b->end))
    //return 0; /* range_a == range_b */

  //if ((range_a->start < range_b->start) ||
      //((range_a->start == range_b->start) && (range_a->end < range_b->end)))
    //return -1; /* range_a < range_b */

  //return 1; /* range_a > range_b */
//}


static void matchposinsinglesequence(Matchprocessinfo *matchprocessinfo,
                                      uint64_t unitnum,
                                      const GtUchar *query,
                                      unsigned long querylen,
                                      const char *querydesc)
{
  const GtUchar *qptr;
  unsigned long remaining;
  unsigned long subjectpos, *sptr;
  bool hasmatch;

  if (matchprocessinfo->preprocessmatchfunction != NULL)
  {
    matchprocessinfo->preprocessmatchfunction(unitnum,
                                          querylen,
                                          querydesc,
                                          matchprocessinfo->showspecinfo);
  }
  if (matchprocessinfo->encseq != NULL)
  {
    sptr = &subjectpos;
  } 
  		
    
  //if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM) 
	//{ 
		//for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
		//{
			//hasmatch = matchprocessinfo->processmatchfunction((const BWTSeq *) matchprocessinfo->genericindex,
																			 //matchprocessinfo->encseq,
																			 //matchprocessinfo->totallength,
																			 //(matchprocessinfo->leastlength).valueunsignedlong,
																			 //query,                  // *query      absolute query start position
																			 //qptr,                   // *qstart     point position in query (qptr will be variable from the point) 
																			 //query+querylen,         // *qend       absolute query end position
																			 //matchprocessinfo->mumcandtab
																			 //);	
			//if ( hasmatch )		
			//{
				//matchprocessinfo->postprocessmatchfunction(matchprocessinfo->alphabet,
																					 //query,
																					 //(unsigned long) (qptr-query),
																					 //querylen,
																					 //matchprocessinfo);
			//}																	
		//}
	//}
	//else 
	//if ( (matchprocessinfo->matchmode == GT_MATCHMODE_MUMREFERENCE) || (matchprocessinfo->matchmode == GT_MATCHMODE_MAXMATCH) ) 
	//{  
		// 对于每一条 query 有一个循环 ＝> 这个query的所有后缀
		for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
		{
			hasmatch = matchprocessinfo->processmatchfunction((const BWTSeq *) matchprocessinfo->genericindex,
																			 matchprocessinfo->encseq,
																			 matchprocessinfo->totallength,
																			 (matchprocessinfo->leastlength).valueunsignedlong,
																			 //sptr,                   // *subjectpos
																			 query,                  // *query      absolute query start position
																			 qptr,                   // *qstart     point position in query (qptr will be variable from the point) 
																			 query+querylen,         // *qend       absolute query end position
																			 matchprocessinfo->maximalmatchtab
																			 );	
																			 
																			 

				if ( hasmatch )		
				{
					// in case GT_MATCHMODE_MUMREFERENCE(>=1 records) or GT_MATCHMODE_MUMREFERENCE and GT_MATCHMODE_MUM(==1 record) the maximalmatchtab is not empty 
					while (gt_array_size(matchprocessinfo->maximalmatchtab)!=0) {
						Maximalmatch *mm = (Maximalmatch *)gt_array_pop(matchprocessinfo->maximalmatchtab); 
						if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM)
		        {	    
							/*
								The following code stores the information about a MUM-candidate
								in the next free position of the dynamic array matchprocessinfo->mumcandtab.
							*/
							MUMcandidate mumcandidate;
							mumcandidate.mumlength = mm->matchlength;
							mumcandidate.dbstart = mm->dbstart;
							//mumcandidate.queryseq = unitnum;
							mumcandidate.querystart = qptr;					    
							gt_array_add(matchprocessinfo->mumcandtab, mumcandidate);
				    }
				    else
				    {
							//unsigned long matchlength = mm->matchlength;
							//unsigned long subjectpos = mm->dbstart;     // 名字不统一 dbstart另一个是subjectpos
							matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
																		 matchprocessinfo->alphabet,
																		 query,
																		 (unsigned long) (qptr-query),
																		 querylen,
																		 mm->matchlength,
																		 mm->dbstart,
																		 matchprocessinfo->showspecinfo);
					  }
					
				  }
		    }																	
		}
		
		
		//printf ("# size=%lu \n", gt_array_size(matchprocessinfo->mumcandtab));
		if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM) 
	  {	
		    matchprocessinfo->postprocessmatchfunction(matchprocessinfo, query, querylen);
		    gt_array_reset(matchprocessinfo->mumcandtab);
		}


	


}

/*
  The following function shows the query sequence description of sequence
  number 'unitnum'.
*/
static void showquerydesc(GT_UNUSED uint64_t unitnum,
                        unsigned long querylen,
                        const char *querydesc,
                        void *info)  
{
	Showspecinfo *showspecinfo = (Showspecinfo *) info;
  //printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
	////char *pch = strchr(querydesc,' ');
	////unsigned long querydesclength = (unsigned long)(pch-querydesc);
  if (querydesc != NULL && querydesc[0] != '\0')
  {
		////char *buf = gt_calloc(1, sizeof (char) * (querydesclength +1));  
    ////(void) strncpy(buf, querydesc, querydesclength);
    
		if (showspecinfo->showsequencelengths)
		{
			if (showspecinfo->queryreadmode==GT_READMODE_FORWARD)
		  {
			  printf("> %s  Len = %lu",querydesc,querylen);
			}
			else
			{
				printf("> %s Reverse  Len = %lu",querydesc,querylen);
			}
			//printf("> %s  Len = %lu",buf,querylen); // first line is query sequence
		} 
		else
		{
			printf("> %s",querydesc);
		} 
		////gt_free(buf);
  }
  printf("\n");
}


static void showmaximalmatch(const GtEncseq *encseq,
                                const GtAlphabet *alphabet,
                                const GtUchar *start,                 
                                unsigned long querypos,
                                unsigned long querylength,
                                unsigned long matchlength,
                                unsigned long subjectpos,
                                void *info)
{		
	  Showspecinfo *showspecinfo = (Showspecinfo *) info;
		unsigned long seqnum = gt_encseq_seqnum(encseq, subjectpos);
		subjectpos = subjectpos - gt_encseq_seqstartpos(encseq, seqnum);
		unsigned long seqtotalnum = gt_encseq_num_of_sequences(encseq);			    

		const char *referencedesc;   
		unsigned long referencedesclength;
		referencedesc = gt_encseq_description(encseq, &referencedesclength, seqnum);
		char *pch = strchr(referencedesc,' ');
		referencedesclength = (unsigned long)(pch-referencedesc);
							
		if (referencedesc != NULL && referencedesc[0] != '\0' && seqtotalnum!=1)
		{
			char *buf = gt_calloc(1, sizeof (char) * (referencedesclength +1));  
			(void) strncpy(buf, referencedesc, referencedesclength);
			printf("  %s",buf);   
			gt_free(buf);
		}

		printf("   %8lu  ",subjectpos+1);
		if (showspecinfo->showreversepositions)
		{	
			if (showspecinfo->queryreadmode==GT_READMODE_REVCOMPL)
			{
				printf("%8lu  ",querylength-querypos);
			}
			else
			{
				printf("%8lu  ",querypos+1);
			}
		} 
		else
		{
			printf("%8lu  ",querypos+1);
		}
		printf("%8lu\n",matchlength);
		if (showspecinfo->showstring)
		{
			gt_alphabet_decode_seq_to_fp(alphabet,stdout,start + querypos,
																	 matchlength);
			(void) putchar('\n');															 
		}
}

/*
	Output all MUM candidates that are unique in the query sequence.
	These are the MUMs. The MUM-candidates are stored in table
	mumcandtab. The MUM is processed further by the function
	showmatchfunction.
*/
static void mumuniqueinquery(void *info,                                       
                                      const GtUchar *query,
                                      unsigned long querylen)
{	
	  Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
		if(gt_array_size(matchprocessinfo->mumcandtab) > 0)
		{
			unsigned int currentright, dbright = 0;
			MUMcandidate *mumcandptr;
			bool ignorecurrent, ignoreprevious = false;
			
			/*
				Sort all MUM-candidates according by increasing dbstart-value
				and decreasing length.
			*/
			gt_array_sort_stable(matchprocessinfo->mumcandtab, (GtCompare)compareMUMcandidates); 
			int i;			
			for (i = 0; i < gt_array_size(matchprocessinfo->mumcandtab); i++) 
			{
				mumcandptr = (MUMcandidate *)gt_array_get(matchprocessinfo->mumcandtab, i);

				ignorecurrent = false;
				currentright = mumcandptr->dbstart + mumcandptr->mumlength - 1;
				if(dbright > currentright)
				{
					// 排除真包含的情况
					ignorecurrent = true;
				} else
				{
					if(dbright == currentright)
					{
						ignorecurrent = true;
						// ( (如果前一个没有被 排除 ) && (两个mumcand 对应同一段 dbsubstring=>这两个mumcand重复了) ) => 排除前一个
						if(!ignoreprevious && (mumcandptr-1)->dbstart == mumcandptr->dbstart)
						{
							ignoreprevious = true;
						}
					} else
					{
						dbright = currentright;
					}
				}
				if( (mumcandptr > (MUMcandidate *)gt_array_get_first(matchprocessinfo->mumcandtab)) && !ignoreprevious)  // mumcandptr > mumcand->spaceMUMcandidate && 
				{
					matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
																 matchprocessinfo->alphabet,
																 query,
																 (unsigned long) ((mumcandptr-1)->querystart-query),
																 querylen,
																 (mumcandptr-1)->mumlength,
																 (mumcandptr-1)->dbstart,
																 matchprocessinfo->showspecinfo);
				}
				ignoreprevious = ignorecurrent;
			}
			
			// 最后一个
			if(!ignoreprevious)
			{
				mumcandptr = (MUMcandidate *)gt_array_get_last(matchprocessinfo->mumcandtab);											
				matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
																 matchprocessinfo->alphabet,
																 query,
																 (unsigned long) (mumcandptr->querystart-query),
																 querylen,
																 mumcandptr->mumlength,
																 mumcandptr->dbstart,
																 matchprocessinfo->showspecinfo);			
			}
								
		}	
}	


int gt_findmum(const GtEncseq *encseq,
                              const void *genericindex,
                              unsigned long totallength,
                              const GtAlphabet *alphabet,
                              const GtStrArray *queryfilenames,
                              const GtMatchmode matchmode,
                              Definedunsignedlong leastlength,
                              bool nucleotidesonly,                     
                              bool bothdirections,                    
                              bool reversecomplement,                    
                              bool showstring,                           
                              bool showreversepositions,                     
                              bool showsequencelengths,   
                              GT_UNUSED bool verbose,                  
                              GtError *err)
{
  Matchprocessinfo matchprocessinfo;
  Showspecinfo showspecinfo;
  bool haserr = false;
  GtSeqIterator *seqit;
  const GtUchar *query;
  unsigned long querylen;
  char *querydesc = NULL;
  int retval;
  uint64_t unitnum;

  gt_error_check(err);
  matchprocessinfo.genericindex = genericindex;
  matchprocessinfo.totallength = totallength;
  matchprocessinfo.matchmode = matchmode;
  
  matchprocessinfo.preprocessmatchfunction = showquerydesc;
  if (matchprocessinfo.matchmode == GT_MATCHMODE_MUM) 
	{ 
      matchprocessinfo.processmatchfunction = gt_packedindexmumreference;
      matchprocessinfo.postprocessmatchfunction = mumuniqueinquery;
	}
	else if (matchprocessinfo.matchmode == GT_MATCHMODE_MUMREFERENCE) 
	{
		  matchprocessinfo.processmatchfunction = gt_packedindexmumreference;
		  matchprocessinfo.postprocessmatchfunction = NULL;
	}  
	else if (matchprocessinfo.matchmode == GT_MATCHMODE_MAXMATCH) 
	{
		  matchprocessinfo.processmatchfunction = gt_packedindexmaxmatch;
		  matchprocessinfo.postprocessmatchfunction = NULL;
  }
  matchprocessinfo.showmatchfunction = showmaximalmatch;
  
  
  matchprocessinfo.alphabet = alphabet;
  //matchprocessinfo.processinfo = &rangespecinfo;
  matchprocessinfo.encseq = encseq;
  

  GtArray *mumcandtab = gt_array_new(sizeof (MUMcandidate));
  GtArray *maximalmatchtab = gt_array_new(sizeof (Maximalmatch));
  matchprocessinfo.mumcandtab = mumcandtab;
  matchprocessinfo.maximalmatchtab = maximalmatchtab;
  matchprocessinfo.leastlength = leastlength;
  matchprocessinfo.showspecinfo = &showspecinfo;
  
  
  showspecinfo.nucleotidesonly = nucleotidesonly;                  
  //rangespecinfo.bothdirections = bothdirections;                  
  //rangespecinfo.reversecomplement = reversecomplement;                   
  showspecinfo.showstring = showstring;                         
  showspecinfo.showreversepositions = showreversepositions;                   
  showspecinfo.showsequencelengths = showsequencelengths;
  // fill the original model of 3 functions, it actually misses 2 functions, the other 2 are processgmatchlength, postprocessgmatchlength

  
  // analyse the queryfilenames with seqiterator
  seqit = gt_seqiterator_sequence_buffer_new(queryfilenames, err);
  if (!seqit)
    haserr = true;
  if (!haserr)
  {
    gt_seqiterator_set_symbolmap(seqit, gt_alphabet_symbolmap(alphabet));
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = gt_seqiterator_next(seqit,
                                &query,
                                &querylen,
                                &querydesc,
                                err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      

			if ( !reversecomplement ) 
			{ 
				showspecinfo.queryreadmode = GT_READMODE_FORWARD;
				matchposinsinglesequence(&matchprocessinfo,
													unitnum,
													query,
													querylen,
													querydesc);
			}
			if ( reversecomplement || bothdirections ) 
			{
				GtUchar *revcompquery = NULL;
				revcompquery = gt_calloc(querylen+1, sizeof (GtUchar)); 
				// copy GtUchar from query to revcompquery
        memcpy(revcompquery, query, querylen*sizeof(GtUchar));
                    
        char *temp_char = gt_calloc(querylen+1, sizeof (char));
        // transformation from GtUchar to char  (revcomquery --> temp_char)
        gt_alphabet_decode_seq_to_cstr(alphabet,temp_char,revcompquery,querylen);
        
        // in situ char-replacement
        haserr = gt_reverse_complement(temp_char, querylen, err);
        if (haserr)
          break;
          
        // transformation from char to GtUchar (temp_char -> revcomquery) 
        gt_alphabet_encode_seq(alphabet, revcompquery, temp_char, querylen);

        gt_free(temp_char);
        showspecinfo.queryreadmode = GT_READMODE_REVCOMPL;
				matchposinsinglesequence(&matchprocessinfo,
													unitnum,
													revcompquery,
													querylen,
													querydesc);
													
				gt_free(revcompquery);	
			} 
      
      FREESPACE(querydesc);
    }
    gt_seqiterator_delete(seqit);
    //gt_array_reset(mumcandtab);
    gt_array_reset(maximalmatchtab);
  }
  gt_array_delete(mumcandtab);
  gt_array_delete(maximalmatchtab);
  return haserr ? -1 : 0;
}

