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

#include "core/seqiterator_sequence_buffer.h"
#include "core/unused_api.h"
#include "spacedef.h"
#include "core/format64.h"
#include "match/eis-voiditf.h"
#include "extended/reverse.h"
#include "match/maxmat4def.h"

/*
  This file contains functions to appropriately call the function 
  gt_voidpackedindexmumreference and gt_voidpackedindexmaxmatch
  and to process their result according to the options given by the user.
*/

/* 
  The following structure stores MUM candidates. That is, maximal matches
  which are unique in the subject-sequence but not necessarily in the
  query sequence. 
*/
typedef struct
{
  unsigned long mumlength,    /* length of the mum */
                subjectpos;   /* start position in the subject-sequence */
  const GtUchar *qstart;      /* start position in the query sequence */
} MUMcandidate;
                                                                       
/*
  The following is the type for the functions to finding maximal matches or
  MUM candidates.
*/
typedef bool (*Findmatchfunction)(const GtUchar *query,
																		const GtUchar *qstart,
																		const GtUchar *qend,
																		Processmatchfunction processmatch,
																		Matchprocessinfo *info);
                                     
/*
  The following function is imported from eis-voiditf.h
*/
bool gt_voidpackedindexmumcandidates(const GtUchar *query,
																		const GtUchar *qstart,
																		const GtUchar *qend,
																		Processmatchfunction processmatch,
																		Matchprocessinfo *info);

/*
  The following function is imported from eis-voiditf.h
*/
bool gt_voidpackedindexmaxmatches(const GtUchar *query,
																	const GtUchar *qstart,
																	const GtUchar *qend,
																	Processmatchfunction processmatch,
																	Matchprocessinfo *info);
																          
/*
  The following function shows the query sequence description of sequence
  number 'unitnum'.
*/
static void showquerydesc(GT_UNUSED uint64_t unitnum,
                        unsigned long querylen,
                        const char *querydesc,
                        Showspecinfo *showspecinfo)
{
  //Showspecinfo *showspecinfo = (Showspecinfo *) info;
  /* printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum)); */
  if (querydesc != NULL && querydesc[0] != '\0')
  {
    if (showspecinfo->showsequencelengths)
    {
      /* first line is query sequence */
      if (showspecinfo->queryreadmode==GT_READMODE_FORWARD)
      {
        printf("> %s  Len = %lu",querydesc,querylen);
      }
      else
      {
        printf("> %s Reverse  Len = %lu",querydesc,querylen);
      }
    }
    else
    {
      printf("> %s",querydesc);
    }
  }
  printf("\n");
}  

/*
  The following is a function to further process mum-candidate,
  maxmatch or mum specfied by its length, its start position
  in the subject sequence as well as the start of the match in the query.
  Then it shows the relevant information as a triple of three integers.
*/
static short int showmatch(const GtEncseq *encseq,
                             const GtUchar *query,
                             unsigned long querypos,
                             unsigned long querylength,
                             unsigned long matchlength,
                             unsigned long subjectpos,
                             Showspecinfo *showspecinfo)
{
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
    gt_alphabet_decode_seq_to_fp(gt_encseq_alphabet(encseq),stdout,
                                 query + querypos,matchlength);
    (void) putchar('\n');
  }
  return 0;
}

/*
  The following function stores the information about a MUM-candidate
  in the next free position of the dynamic array mumcandtab.
*/
static short int storeMUMcandidate (GT_UNUSED const GtEncseq *encseq,
																	 const GtUchar *query,
																	 unsigned long querypos,
																	 GT_UNUSED unsigned long querylength,
																	 unsigned long matchlength,
																	 unsigned long subjectpos,
																	 Showspecinfo *showspecinfo)
{
  //Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;

	/*
		printf("storeMUMcandiate %lu %lu %lu\n",
							matchlength,
							subjectpos,
							//seqnum,
							query);
	*/
	/*
	 * save the result in the dynamic array mumcandtab.
	 */
	MUMcandidate mumcand;;
  mumcand.mumlength = matchlength;
  mumcand.subjectpos = subjectpos;
  /*mumcand.queryseq = unitnum;*/
  mumcand.qstart = query + querypos;
  gt_array_add(showspecinfo->mumcandtab, mumcand);
  return 0;
}
                                   
/*
  The following function compares two MUM-candidates. The MUM-candidate
  with smaller subjectpos-value comes first.
  If both MUMs have the same subjectpos-value, then the MUM-candidate
  with the larger length comes first.
*/
static int compareMUMcandidates(MUMcandidate *p,MUMcandidate *q)
{
  if (p->subjectpos == q->subjectpos)
  {
    return (p->mumlength < q->mumlength) ? 1 : -1;
  }
  return (p->subjectpos > q->subjectpos) ? 1 : -1;
}

/*
  Output all MUM candidates that are unique in the query sequence.
  These are the MUMs. The MUM-candidates are stored in table
  mumcandtab. The MUM is processed further by the function
  showmatchfunction.
*/
static void mumuniqueinquery(Matchprocessinfo *matchprocessinfo,
                             const GtUchar *query,
                             unsigned long querylen)
{
  //Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
  if (gt_array_size(matchprocessinfo->showspecinfo->mumcandtab) > 0)
  {
    unsigned int currentright, dbright = 0;
    MUMcandidate *mumcandptr;
    bool ignorecurrent, ignoreprevious = false;

    /*
      Sort all MUM-candidates according by increasing subjectpos-value
      and decreasing length.
    */
    gt_array_sort_stable(matchprocessinfo->showspecinfo->mumcandtab,
                         (GtCompare)compareMUMcandidates);
    int i;
    for (i = 0; i < gt_array_size(matchprocessinfo->showspecinfo->mumcandtab); i++)
    {
      mumcandptr = (MUMcandidate *)gt_array_get(matchprocessinfo->showspecinfo->mumcandtab,i);

      ignorecurrent = false;
      currentright = mumcandptr->subjectpos + mumcandptr->mumlength - 1;
      if (dbright > currentright)
      {
        ignorecurrent = true;
      } else
      {
        if (dbright == currentright)
        {
          ignorecurrent = true;
          if (!ignoreprevious &&
              (mumcandptr-1)->subjectpos == mumcandptr->subjectpos)
          {
            ignoreprevious = true;
          }
        } else
        {
          dbright = currentright;
        }
      }
      if ( (mumcandptr>\
            (MUMcandidate *)gt_array_get_first(matchprocessinfo->showspecinfo->mumcandtab))
            && !ignoreprevious)
      {
        showmatch(matchprocessinfo->encseq,
								 query,
								 /* qstart - query = querypos <=> querypos + query = qstart */
								 (unsigned long) ((mumcandptr-1)->qstart-query),
								 querylen,
								 (mumcandptr-1)->mumlength,
								 (mumcandptr-1)->subjectpos,
								 matchprocessinfo->showspecinfo);
      }
      ignoreprevious = ignorecurrent;
    }

    /* the last one */
    if (!ignoreprevious)
    {
      mumcandptr =\
          (MUMcandidate *)gt_array_get_last(matchprocessinfo->showspecinfo->mumcandtab);
      showmatch(matchprocessinfo->encseq,
							 query,
							 (unsigned long) (mumcandptr->qstart-query),
							 querylen,
							 mumcandptr->mumlength,
							 mumcandptr->subjectpos,
							 matchprocessinfo->showspecinfo);
    }

  }
}

static void matchposinsinglesequence(uint64_t unitnum,
                                      const GtUchar *query,
                                      unsigned long querylen,
                                      const char *querydesc,
                                      Findmatchfunction findmatchfunction,
                                      Processmatchfunction processmatch,
                                      Matchprocessinfo *matchprocessinfo)
{
  const GtUchar *qptr;
  unsigned long remaining;
  //unsigned long subjectpos;
  //bool hasmatch;

  /* for every query, show its description at the beginning */
  showquerydesc(unitnum,
								querylen,
								querydesc,
								matchprocessinfo->showspecinfo);

  /* take every suffix of query, qptr is a alias name of qstart */
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {		
    /* hasmatch = */findmatchfunction(
                   query,           /* absolute query start position */
                   qptr,            /* variable position in query */
                   query+querylen,  /* absolute query end position */
                   processmatch,
                   matchprocessinfo);
  }

  //printf ("# size=%lu \n", gt_array_size(matchprocessinfo->showspecinfo->mumcandtab)); 
  if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM)
  {
      mumuniqueinquery(matchprocessinfo,
											query,
											querylen);
      gt_array_reset(matchprocessinfo->showspecinfo->mumcandtab);
  }
}

int gt_findmum(const GtEncseq *encseq,
								const void *genericindex,
								unsigned long totallength,
								const GtAlphabet *alphabet,
								const GtStrArray *queryfilenames,
								const GtMatchmode matchmode,
								Definedunsignedlong leastlength,
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
  Findmatchfunction findmatchfunction = NULL;
  Processmatchfunction processmatch = NULL;
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

  //  preprocessmatchfunction = showquerydesc;
  //if (matchprocessinfo.matchmode == GT_MATCHMODE_MUM)
  //{
    //matchprocessinfo.processmatchfunction = gt_voidpackedindexmumreference;
    //matchprocessinfo.postprocessmatchfunction = mumuniqueinquery;
  //}
  //else if (matchprocessinfo.matchmode == GT_MATCHMODE_MUMREFERENCE)
  //{
    //matchprocessinfo.processmatchfunction = gt_voidpackedindexmumreference;
    //matchprocessinfo.postprocessmatchfunction = NULL;
  //}
  //else if (matchprocessinfo.matchmode == GT_MATCHMODE_MAXMATCH)
  //{
    //matchprocessinfo.processmatchfunction = gt_voidpackedindexmaxmatch;
    //matchprocessinfo.postprocessmatchfunction = NULL;
  //}
  //matchprocessinfo.showmatchfunction = showmaximalmatch;
  if(matchprocessinfo.matchmode == GT_MATCHMODE_MUM)
  {
    findmatchfunction = gt_voidpackedindexmumcandidates;
    processmatch = storeMUMcandidate;
	}
  else if (matchprocessinfo.matchmode == GT_MATCHMODE_MUMREFERENCE)
  {
    findmatchfunction = gt_voidpackedindexmumcandidates;
    processmatch = showmatch;
  }
  else if (matchprocessinfo.matchmode == GT_MATCHMODE_MAXMATCH)
  {
		findmatchfunction = gt_voidpackedindexmaxmatches;
    processmatch = showmatch;
	}

  matchprocessinfo.alphabet = alphabet;
  matchprocessinfo.encseq = encseq;
  matchprocessinfo.showspecinfo = &showspecinfo;

  GtArray *mumcandtab = gt_array_new(sizeof (MUMcandidate));
  matchprocessinfo.leastlength = leastlength;

  showspecinfo.showstring = showstring;
  showspecinfo.showreversepositions = showreversepositions;
  showspecinfo.showsequencelengths = showsequencelengths;
  showspecinfo.mumcandtab = mumcandtab;

  /* analyse the queryfilenames with seqiterator */
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
        matchposinsinglesequence(unitnum,
                                 query,
                                 querylen,
                                 querydesc,
                                 findmatchfunction,
                                 processmatch,
                                 &matchprocessinfo);
      }
      if ( reversecomplement || bothdirections )
      {
        GtUchar *revcompquery = NULL;
        revcompquery = gt_calloc(querylen+1, sizeof (GtUchar));
        /* copy GtUchar from query to revcompquery */
        memcpy(revcompquery, query, querylen*sizeof (GtUchar));

        char *temp_char = gt_calloc(querylen+1, sizeof (char));
        /* transformation from GtUchar to char  (revcomquery --> temp_char) */
        gt_alphabet_decode_seq_to_cstr(alphabet,
                                       temp_char,
                                       revcompquery,
                                       querylen);

        /* in situ char-replacement */
        haserr = gt_reverse_complement(temp_char, querylen, err);
        if (haserr)
          break;

        /* transformation from char to GtUchar (temp_char -> revcomquery) */
        gt_alphabet_encode_seq(alphabet, revcompquery, temp_char, querylen);

        gt_free(temp_char);
        showspecinfo.queryreadmode = GT_READMODE_REVCOMPL;
        matchposinsinglesequence(unitnum,
                                 revcompquery,
                                 querylen,
                                 querydesc,
                                 findmatchfunction,
                                 processmatch,
                                 &matchprocessinfo);

        gt_free(revcompquery);
      }

      FREESPACE(querydesc);
    }
    gt_seqiterator_delete(seqit);
  }
  gt_array_delete(mumcandtab);
  return haserr ? -1 : 0;
}
