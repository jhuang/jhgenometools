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



typedef void (*Preprocessgmatchlength)(uint64_t,
                                       unsigned long querylen,
                                       const char *,
                                       void *);

typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  GtMatchmode matchmode;
  const GtAlphabet *alphabet;
  Preprocessgmatchlength preprocessgmatchlength;
  void *processinfo;
  const GtEncseq *encseq;
} Substringinfo;


static void matchposinsinglesequence(Substringinfo *substringinfo,
                                      uint64_t unitnum,
                                      const GtUchar *query,
                                      unsigned long querylen,
                                      const char *querydesc)
{
  const GtUchar *qptr;
  unsigned long gmatchlength, remaining;
  unsigned long witnessposition, *wptr;

  if (substringinfo->preprocessgmatchlength != NULL)
  {
    substringinfo->preprocessgmatchlength(unitnum,
                                          querylen,
                                          querydesc,
                                          substringinfo->processinfo);
  }
  if (substringinfo->encseq != NULL)
  {
    wptr = &witnessposition;
  } 
  		
    
  if (substringinfo->matchmode == GT_MATCHMODE_MUM) 
	{ 
		printf ("%s \n", "match mode 'mum' is still in work");
	}
	else if (substringinfo->matchmode == GT_MATCHMODE_MUMREFERENCE) 
	{  
		for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
		{
			gmatchlength = gt_packedindexmumreference((const BWTSeq *) substringinfo->genericindex,
																			 substringinfo->encseq,
																			 substringinfo->alphabet,
																			 substringinfo->totallength,
																			 wptr,
																			 query,
																			 qptr,
																			 query+querylen,
																			 //substringinfo->queryreadmode,
																			 substringinfo->processinfo);
																		
		}
	}
	else if (substringinfo->matchmode == GT_MATCHMODE_MAXMATCH) 
	{
		for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
		{
			gmatchlength = gt_packedindexmaxmatch((const BWTSeq *) substringinfo->genericindex,
																	 substringinfo->encseq,
																	 substringinfo->alphabet,
																	 substringinfo->totallength,
																	 wptr,
																	 query,
																	 qptr,
																	 query+querylen,
																	 //substringinfo->queryreadmode,
																	 substringinfo->processinfo);
		}
	}  
}

static void showunitnum(GT_UNUSED uint64_t unitnum,
                        unsigned long querylen,
                        const char *querydesc,
                        GT_UNUSED void *info)  // info is type Rangespecinfo 
{
	Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;
  //printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
	////char *pch = strchr(querydesc,' ');
	////unsigned long querydesclength = (unsigned long)(pch-querydesc);
  if (querydesc != NULL && querydesc[0] != '\0')
  {
		////char *buf = gt_calloc(1, sizeof (char) * (querydesclength +1));  
    ////(void) strncpy(buf, querydesc, querydesclength);
		if (rangespecinfo->showsequencelengths)
		{
			if (rangespecinfo->queryreadmode==GT_READMODE_FORWARD)
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
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;
  bool haserr = false;
  GtSeqIterator *seqit;
  const GtUchar *query;
  unsigned long querylen;
  char *querydesc = NULL;
  int retval;
  uint64_t unitnum;

  gt_error_check(err);
  substringinfo.genericindex = genericindex;
  substringinfo.totallength = totallength;
  substringinfo.matchmode = matchmode;
  rangespecinfo.leastlength = leastlength;
  rangespecinfo.nucleotidesonly = nucleotidesonly;                  
  //rangespecinfo.bothdirections = bothdirections;                  
  //rangespecinfo.reversecomplement = reversecomplement;                   
  rangespecinfo.showstring = showstring;                         
  rangespecinfo.showreversepositions = showreversepositions;                   
  rangespecinfo.showsequencelengths = showsequencelengths;
  // fill the original model of 3 functions, it actually misses 2 functions, the other 2 are processgmatchlength, postprocessgmatchlength
  substringinfo.preprocessgmatchlength = showunitnum;
  
  substringinfo.alphabet = alphabet;
  substringinfo.processinfo = &rangespecinfo;
  substringinfo.encseq = encseq;
  
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
				rangespecinfo.queryreadmode = GT_READMODE_FORWARD;
				matchposinsinglesequence(&substringinfo,
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
        rangespecinfo.queryreadmode = GT_READMODE_REVCOMPL;
				matchposinsinglesequence(&substringinfo,
													unitnum,
													revcompquery,
													querylen,
													querydesc);
													
				gt_free(revcompquery);	
			} 
      
      FREESPACE(querydesc);
    }
    gt_seqiterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}

