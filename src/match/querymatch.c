/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/symboldef.h"
#include "readmode-def.h"
#include "seqpos-def.h"
#include "querymatch.h"
#include "format64.h"

struct Querymatch
{
   Seqpos len,
          dbstart,
          querystart,
          querytotallength;
   bool selfmatch;
   uint64_t queryseqnum;
   Readmode readmode;
};

Querymatch *querymatch_new(void)
{
  return gt_malloc(sizeof(Querymatch));
}

void querymatch_fill(Querymatch *querymatch,
                     Seqpos len,
                     Seqpos dbstart,
                     Readmode readmode,
                     bool selfmatch,
                     uint64_t queryseqnum,
                     Seqpos querystart,
                     Seqpos querytotallength)
{
  querymatch->len = len;
  querymatch->dbstart = dbstart;
  querymatch->readmode = readmode;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querystart = querystart;
  querymatch->querytotallength = querytotallength;
}

void querymatch_delete(Querymatch *querymatch)
{
  if (querymatch != NULL)
  {
    gt_free(querymatch);
  }
}

#ifdef VERIFY
static void verifymatch(const Encodedsequence *encseq,
                        Seqpos len,
                        Seqpos pos1,
                        uint64_t seqnum2,
                        Seqpos pos2,
                        Readmode readmode)
{
  if (readmode == Reversemode)
  {
    Seqinfo seqinfo;
    Seqpos offset, totallength = getencseqtotallength(encseq);
    GtUchar cc1, cc2;

    getencseqSeqinfo(&seqinfo,encseq,(Seqpos) seqnum2);
    pos2 += seqinfo.seqstartpos;
    for (offset = 0; offset < len; offset++)
    {
      gt_assert(pos1 + len - 1 < totallength);
      gt_assert(pos2 + len - 1 < totallength);
      cc1 = getencodedchar(encseq,pos1+offset,Forwardmode);
      cc2 = getencodedchar(encseq,pos2+len-1-offset,Forwardmode);
      gt_assert(cc1 == cc2 && ISNOTSPECIAL(cc1));
    }
    if (pos1 + len < totallength)
    {
      cc1 = getencodedchar(encseq,pos1+len,Forwardmode);
    } else
    {
      cc1 = SEPARATOR;
    }
    if (pos2 > 0)
    {
      cc2 = getencodedchar(encseq,pos2-1,Forwardmode);
    } else
    {
      cc2 = SEPARATOR;
    }
    gt_assert(cc1 != cc2 || ISSPECIAL(cc1));
  }
}
#endif

int querymatch_output(GT_UNUSED void *info,
                      const Encodedsequence *encseq,
                      const Querymatch *querymatch,
                      GT_UNUSED GtError *err)
{
  const char *outflag = "FRCP";
  Seqinfo seqinfo;
  unsigned long dbseqnum;
  Seqpos querystart, dbstart_relative;

  gt_assert(encseq != NULL);
  dbseqnum = getencseqfrompos2seqnum(encseq,querymatch->dbstart);
  getencseqSeqinfo(&seqinfo,encseq,dbseqnum);
  gt_assert((int) querymatch->readmode < 4);
  if (querymatch->readmode == Reversemode ||
      querymatch->readmode == Reversecomplementmode)
  {
    gt_assert(querymatch->querystart + querymatch->len <=
              querymatch->querytotallength);
    querystart = querymatch->querytotallength -
                 querymatch->querystart - querymatch->len;
  } else
  {
    querystart = querymatch->querystart;
  }
  gt_assert(querymatch->dbstart >= seqinfo.seqstartpos);
  dbstart_relative = querymatch->dbstart - seqinfo.seqstartpos;
  if (!querymatch->selfmatch ||
      (uint64_t) dbseqnum != querymatch->queryseqnum ||
      dbstart_relative <= querystart)
  {
#ifdef VERIFY
    verifymatch(encseq,
                querymatch->len,
                querymatch->dbstart,
                querymatch->queryseqnum,
                querystart,
                querymatch->readmode);
#endif
    printf(FormatSeqpos " %lu " FormatSeqpos " %c " FormatSeqpos
           " " Formatuint64_t " " FormatSeqpos "\n",
           PRINTSeqposcast(querymatch->len),
           dbseqnum,
           PRINTSeqposcast(dbstart_relative),
           outflag[querymatch->readmode],
           PRINTSeqposcast(querymatch->len),
           PRINTuint64_tcast(querymatch->queryseqnum),
           PRINTSeqposcast(querystart));
  }
  return 0;
}

Seqpos querymatch_len(const Querymatch *querymatch)
{
  return querymatch->len;
}

Seqpos querymatch_dbstart(const Querymatch *querymatch)
{
  return querymatch->dbstart;
}

Seqpos querymatch_querystart(const Querymatch *querymatch)
{
  return querymatch->querystart;
}

uint64_t querymatch_queryseqnum(const Querymatch *querymatch)
{
  return querymatch->queryseqnum;
}