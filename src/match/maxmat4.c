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


typedef struct
{
  bool nucleotidesonly,                     
       //bothdirections,                    
       //reversecomplement,                    
       showstring,                           
       showreversepositions,                     
       showsequencelengths;   
  Definedunsignedlong leastlength;
} Rangespecinfo;

typedef void (*Preprocessgmatchlength)(uint64_t,
                                       const char *,
                                       void *);
typedef void (*Processgmatchlength)(const GtAlphabet *,
                                    const GtUchar *,
                                    unsigned long,
                                    unsigned long,
                                    unsigned long,
                                    void *);
//typedef void (*Postprocessgmatchlength)(const GtAlphabet *,
                                        //uint64_t,
                                        //const char *,
                                        //const GtUchar *,
                                        //unsigned long,
                                        //void *);

typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  const GtAlphabet *alphabet;
  Preprocessgmatchlength preprocessgmatchlength;
  Processgmatchlength processgmatchlength;
  //Postprocessgmatchlength postprocessgmatchlength;
  void *processinfo;
  const GtEncseq *encseq;
} Substringinfo;


// fill the substringinfo after the greedy match pos. alle dbsequence has already merged in one single sequence
static void matchposinsinglesequence(Substringinfo *substringinfo,
                                      uint64_t unitnum,
                                      const GtUchar *query,
                                      unsigned long querylen,
                                      const char *desc,                                 
                                      GtReadmode readmode,
                                      unsigned long leastlength)
{
  const GtUchar *qptr;
  unsigned long gmatchlength, remaining;
  unsigned long witnessposition, *wptr;

  if (substringinfo->preprocessgmatchlength != NULL)
  {
    substringinfo->preprocessgmatchlength(unitnum,
                                          desc,
                                          substringinfo->processinfo);
  }
  //if (substringinfo->encseq != NULL)
  //{
    wptr = &witnessposition;
  //} 
  
  /*
  // why don't have to free memory for the struct?
  Querysubstring querysubstring;
  querysubstring.queryrep = query;
  for (querysubstring.offset = 0;
       querysubstring.offset <= queryrep->length - leastlength;
       querysubstring.offset++)  */
      
  printf ("%s \n", "querystart, matchlength, subjectpos, sequecnce");     
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    gmatchlength = gt_packedindexmum((const BWTSeq *) substringinfo->genericindex,
                                   substringinfo->encseq,
                                   substringinfo->alphabet,
                                   substringinfo->totallength,
                                  wptr,
                                  query,
                                  qptr,
                                  query+querylen,
                                  readmode,
                                  leastlength);
                                  
    //if (gmatchlength > 0)
    //{
      //substringinfo->processgmatchlength(substringinfo->alphabet,
                                         //query,
                                         //gmatchlength,
                                         //(unsigned long) (qptr-query),  // querystart
                                         //witnessposition,
                                         //substringinfo->processinfo);
    //}
  }
  //if (substringinfo->postprocessgmatchlength != NULL)
  //{
    //substringinfo->postprocessgmatchlength(substringinfo->alphabet,
                                           //unitnum,
                                           //desc,
                                           //query,
                                           //querylen,
                                           //substringinfo->processinfo);
  //}
}

static void showunitnum(uint64_t unitnum,
                        const char *desc,
                        GT_UNUSED void *info)
{
  printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
  if (desc != NULL && desc[0] != '\0')
  {
    printf(" (%s)",desc);
  }
  printf("\n");
}

static void showifinlengthrange(const GtAlphabet *alphabet,
                                const GtUchar *start,
                                unsigned long gmatchlength,
                                unsigned long querystart,
                                unsigned long subjectpos,
                                void *info)
{
  Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;
  if (gmatchlength >= rangespecinfo->leastlength.valueunsignedlong) 
  {
	printf ("%s \n", "querystart, matchlength, subjectpos, sequecnce");
    //if (rangespecinfo->showstring)
    //{
      printf("%lu ",querystart);
    //}
    printf("%lu",gmatchlength);
    //if (rangespecinfo->showreversepositions)
    //{
      printf(" %lu",subjectpos);
    //}
    //if (rangespecinfo->showsequencelengths)
    //{
      (void) putchar(' ');
      gt_alphabet_decode_seq_to_fp(alphabet,stdout,start + querystart,
                                   gmatchlength);
    //}
    (void) putchar('\n');
  }
}

int gt_findmum(const GtEncseq *encseq,
                              const void *genericindex,
                              unsigned long totallength,
                              const GtAlphabet *alphabet,
                              const GtStrArray *queryfilenames,
                              GtReadmode readmode,
                              Definedunsignedlong leastlength,
                              bool nucleotidesonly,                     
                              //bool bothdirections,                    
                              //bool reversecomplement,                    
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
  char *desc = NULL;
  int retval;
  uint64_t unitnum;

  gt_error_check(err);
  substringinfo.genericindex = genericindex;
  substringinfo.totallength = totallength;
  rangespecinfo.leastlength = leastlength;
  rangespecinfo.nucleotidesonly = nucleotidesonly;                  
  //rangespecinfo.bothdirections = bothdirections;                  
  //rangespecinfo.reversecomplement = reversecomplement;                   
  rangespecinfo.showstring = showstring;                         
  rangespecinfo.showreversepositions = showreversepositions;                   
  rangespecinfo.showsequencelengths = showsequencelengths;
  // fill the original model of 3 functions
  substringinfo.preprocessgmatchlength = showunitnum;
  substringinfo.processgmatchlength = showifinlengthrange;
  //substringinfo.postprocessgmatchlength = NULL;
  
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
                                &desc,
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
      matchposinsinglesequence(&substringinfo,
                                unitnum,
                                query,
                                querylen,
                                desc,
                                readmode,
                                leastlength.valueunsignedlong );
      FREESPACE(desc);
    }
    gt_seqiterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}

