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
typedef void (*Postprocessmatchfunction)(const GtAlphabet *,
                                const GtUchar *,
                                //unsigned long,
                                unsigned long,
                                unsigned long,
                                //unsigned long,
                                void *);
                                
   


                                
typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  GtMatchmode matchmode;
  const GtAlphabet *alphabet;
  Preprocessmatchfunction preprocessmatchfunction;
  Processmatchfunction processmatchfunction;
  Postprocessmatchfunction postprocessmatchfunction;
  //void *processinfo;
  const GtEncseq *encseq;
  
  GtArray *mumcandtab;
  GtArray *maximalmatchtab;
  bool nucleotidesonly,                                      
       showstring,                           
       showreversepositions,                     
       showsequencelengths;   
  Definedunsignedlong leastlength;
  GtReadmode queryreadmode;
} Matchprocessinfo;  



/*
  Functions processing a maximal match are of the following type.
*/

//typedef Sint (*Processmatchfunction)
             //(void *,Uint,Uint,Uint,Uint); 

                                        
                                        
                                        



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
                                          matchprocessinfo);
  }
  if (matchprocessinfo->encseq != NULL)
  {
    sptr = &subjectpos;
  } 
  		
    
  if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM) 
	{ 
		printf ("%s \n", "match mode 'mum' is still in work");
	}
	else if ( (matchprocessinfo->matchmode == GT_MATCHMODE_MUMREFERENCE) || (matchprocessinfo->matchmode == GT_MATCHMODE_MAXMATCH) ) 
	{  
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
		//printf ("# matchlength=%lu \n", matchlength);
		//if (matchlength >= (matchprocessinfo->leastlength).valueunsignedlong)
			if ( hasmatch )		
			{
				matchprocessinfo->postprocessmatchfunction(matchprocessinfo->alphabet,
																					 query,
																					 //matchlength,
																					 (unsigned long) (qptr-query),
																					 querylen,
																					 //subjectpos,
																					 matchprocessinfo);
			}																	
		}
	} 
}

static void showunitnum(GT_UNUSED uint64_t unitnum,
                        unsigned long querylen,
                        const char *querydesc,
                        void *info)  
{
	Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
  //printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
	////char *pch = strchr(querydesc,' ');
	////unsigned long querydesclength = (unsigned long)(pch-querydesc);
  if (querydesc != NULL && querydesc[0] != '\0')
  {
		////char *buf = gt_calloc(1, sizeof (char) * (querydesclength +1));  
    ////(void) strncpy(buf, querydesc, querydesclength);
		if (matchprocessinfo->showsequencelengths)
		{
			if (matchprocessinfo->queryreadmode==GT_READMODE_FORWARD)
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


static void output(const GtAlphabet *alphabet,
                                const GtUchar *start,
                                //unsigned long matchlength,
                                unsigned long querypos,
                                unsigned long querylength,
                                //unsigned long subjectpos,
                                void *info)
{
	Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
  //int i;						
	//for (i = 0; i < gt_array_size(matchprocessinfo->maximalmatchtab); i++) {
	while (gt_array_size(matchprocessinfo->maximalmatchtab)!=0) {
		//unsigned long matchlength = ((Maximalmatch *)gt_array_get(matchprocessinfo->maximalmatchtab, i))->matchlength;
		//unsigned long subjectpos = ((Maximalmatch *)gt_array_get(matchprocessinfo->maximalmatchtab, i))->dbstart;  // 名字不统一 dbstart另一个是subjectpos
		Maximalmatch *mm = (Maximalmatch *)gt_array_pop(matchprocessinfo->maximalmatchtab);
		unsigned long matchlength = mm->matchlength;
		unsigned long subjectpos = mm->dbstart;
		
		unsigned long seqnum = gt_encseq_seqnum(matchprocessinfo->encseq, subjectpos);
		subjectpos = subjectpos - gt_encseq_seqstartpos(matchprocessinfo->encseq, seqnum);
		unsigned long seqtotalnum = gt_encseq_num_of_sequences(matchprocessinfo->encseq);			    

		const char *referencedesc;   
		unsigned long referencedesclength;
		referencedesc = gt_encseq_description(matchprocessinfo->encseq, &referencedesclength, seqnum);
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
		if (matchprocessinfo->showreversepositions)
		{	
			if (matchprocessinfo->queryreadmode==GT_READMODE_REVCOMPL)
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
		if (matchprocessinfo->showstring)
		{
			gt_alphabet_decode_seq_to_fp(alphabet,stdout,start + querypos,
																	 matchlength);
			(void) putchar('\n');															 
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
  //Rangespecinfo rangespecinfo;
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
  matchprocessinfo.preprocessmatchfunction = showunitnum;
  if (matchprocessinfo.matchmode == GT_MATCHMODE_MUM) 
	{ 
      //matchprocessinfo.processmatchfunction = gt_packedindexmum;
	}
	else if (matchprocessinfo.matchmode == GT_MATCHMODE_MUMREFERENCE) 
	{
		  matchprocessinfo.processmatchfunction = gt_packedindexmumreference;
	}  
	else if (matchprocessinfo.matchmode == GT_MATCHMODE_MAXMATCH) 
	{
		  matchprocessinfo.processmatchfunction = gt_packedindexmaxmatch;
  }
  matchprocessinfo.postprocessmatchfunction = output;
  
  matchprocessinfo.alphabet = alphabet;
  //matchprocessinfo.processinfo = &rangespecinfo;
  matchprocessinfo.encseq = encseq;
  

  GtArray *mumcandtab = gt_array_new(sizeof (MUMcandidate));
  GtArray *maximalmatchtab = gt_array_new(sizeof (Maximalmatch));
  matchprocessinfo.mumcandtab = mumcandtab;
  matchprocessinfo.maximalmatchtab = maximalmatchtab;
  matchprocessinfo.leastlength = leastlength;
  matchprocessinfo.nucleotidesonly = nucleotidesonly;                  
  //rangespecinfo.bothdirections = bothdirections;                  
  //rangespecinfo.reversecomplement = reversecomplement;                   
  matchprocessinfo.showstring = showstring;                         
  matchprocessinfo.showreversepositions = showreversepositions;                   
  matchprocessinfo.showsequencelengths = showsequencelengths;
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
				matchprocessinfo.queryreadmode = GT_READMODE_FORWARD;
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
        matchprocessinfo.queryreadmode = GT_READMODE_REVCOMPL;
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

