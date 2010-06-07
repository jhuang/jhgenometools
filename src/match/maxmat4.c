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
#include "match/eis-bwtseq.h"
#include "extended/reverse.h"

typedef void (*Preprocessmatchfunction)(uint64_t,
                                       unsigned long querylen,
                                       const char *,
                                       void *);

typedef bool (*Processmatchfunction)(const BWTSeq *,
                                     const GtEncseq *,
                                     unsigned long,
                                     unsigned long,
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
  bool showstring,
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
  const GtEncseq *encseq;

  GtArray *mumcandtab;
  GtArray *maximalmatchtab;
  Definedunsignedlong leastlength;
  Showspecinfo *showspecinfo;
} Matchprocessinfo;

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

  /* take every suffix of query, qptr is a alias name of qstart */
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    hasmatch = matchprocessinfo->processmatchfunction(
                   (const BWTSeq *) matchprocessinfo->genericindex,
                   matchprocessinfo->encseq,
                   matchprocessinfo->totallength,
                   (matchprocessinfo->leastlength).valueunsignedlong,
                   query,           /* absolute query start position */
                   qptr,            /* variable position in query */
                   query+querylen,  /* absolute query end position */
                   matchprocessinfo->maximalmatchtab);

      if ( hasmatch )
      {
        /*
         * in case GT_MATCHMODE_MAXMATCH theare are >=1 records,
         * in case GT_MATCHMODE_MUMREFERENCE and GT_MATCHMODE_MUM
         * there is only 1 record in matchprocessinfo->maximalmatchtab
         */
        while (gt_array_size(matchprocessinfo->maximalmatchtab)!=0) {
          Maximalmatch *mm =\
              (Maximalmatch *)gt_array_pop(matchprocessinfo->maximalmatchtab);
          if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM)
          {
            /*
              The following code stores the information about a MUM-candidate
              in the next free position of the dynamic array
              matchprocessinfo->mumcandtab.
            */
            MUMcandidate mumcandidate;
            mumcandidate.mumlength = mm->matchlength;
            mumcandidate.subjectpos = mm->subjectpos;
            /*mumcandidate.queryseq = unitnum;*/
            mumcandidate.qstart = qptr;
            gt_array_add(matchprocessinfo->mumcandtab, mumcandidate);
          }
          else
          {
            /*
             * in case GT_MATCHMODE_MAXMATCH and GT_MATCHMODE_MUMREFERENCE
             * the Maximalmatch is directly printed
             */
            matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
                                   matchprocessinfo->alphabet,
                                   query,
                                   (unsigned long) (qptr-query), /* querypos */
                                   querylen,
                                   mm->matchlength,
                                   mm->subjectpos,
                                   matchprocessinfo->showspecinfo);
          }

        }
      }
  }

  /* printf ("# size=%lu \n", gt_array_size(matchprocessinfo->mumcandtab)); */
  if (matchprocessinfo->matchmode == GT_MATCHMODE_MUM)
  {
      matchprocessinfo->postprocessmatchfunction(matchprocessinfo,
                                                 query,
                                                 querylen);
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
static void showmaximalmatch(const GtEncseq *encseq,
                             const GtAlphabet *alphabet,
                             const GtUchar *query,
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
    gt_alphabet_decode_seq_to_fp(alphabet,stdout,query + querypos,
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
  if (gt_array_size(matchprocessinfo->mumcandtab) > 0)
  {
    unsigned int currentright, dbright = 0;
    MUMcandidate *mumcandptr;
    bool ignorecurrent, ignoreprevious = false;

    /*
      Sort all MUM-candidates according by increasing subjectpos-value
      and decreasing length.
    */
    gt_array_sort_stable(matchprocessinfo->mumcandtab,
                         (GtCompare)compareMUMcandidates);
    int i;
    for (i = 0; i < gt_array_size(matchprocessinfo->mumcandtab); i++)
    {
      mumcandptr = (MUMcandidate *)gt_array_get(matchprocessinfo->mumcandtab,
                                                i);

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
            (MUMcandidate *)gt_array_get_first(matchprocessinfo->mumcandtab))
            && !ignoreprevious)
      {
        matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
                               matchprocessinfo->alphabet,
                               query,
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
          (MUMcandidate *)gt_array_get_last(matchprocessinfo->mumcandtab);
      matchprocessinfo->showmatchfunction(matchprocessinfo->encseq,
                               matchprocessinfo->alphabet,
                               query,
                               (unsigned long) (mumcandptr->qstart-query),
                               querylen,
                               mumcandptr->mumlength,
                               mumcandptr->subjectpos,
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
  matchprocessinfo.encseq = encseq;
  matchprocessinfo.showspecinfo = &showspecinfo;

  GtArray *mumcandtab = gt_array_new(sizeof (MUMcandidate));
  GtArray *maximalmatchtab = gt_array_new(sizeof (Maximalmatch));
  matchprocessinfo.mumcandtab = mumcandtab;
  matchprocessinfo.maximalmatchtab = maximalmatchtab;
  matchprocessinfo.leastlength = leastlength;

  showspecinfo.showstring = showstring;
  showspecinfo.showreversepositions = showreversepositions;
  showspecinfo.showsequencelengths = showsequencelengths;

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
  }
  gt_array_delete(mumcandtab);
  gt_array_delete(maximalmatchtab);
  return haserr ? -1 : 0;
}
