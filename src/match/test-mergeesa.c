/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#include <errno.h>
#include <limits.h>
#include <string.h>
#include "core/fa.h"
#include "sarr-def.h"
#include "spacedef.h"
#include "emimergeesa.h"
#include "esa-fileend.h"
#include "core/logger.h"
#include "lcpoverflow.h"

#include "encseq2offset.pr"

typedef struct
{
  GtStr *outfilename;
  FILE *fp;
} NameandFILE;

typedef struct
{
  NameandFILE outsuf,
              outlcp,
              outllv;
  unsigned long currentlcpindex,
         absstartpostable[SIZEOFMERGERESULTBUFFER];
} Mergeoutinfo;

static int initNameandFILE(NameandFILE *nf,
                            const GtStr *outindex,
                            const char *suffix,
                            GtError *err)
{
  gt_error_check(err);
  nf->outfilename = gt_str_clone(outindex);
  gt_str_append_cstr(nf->outfilename,suffix);
  nf->fp = gt_fa_fopen(gt_str_get(nf->outfilename),"wb",err);
  if (nf->fp == NULL)
  {
    return -1;
  }
  return 0;
}

static void freeNameandFILE(NameandFILE *nf)
{
  gt_fa_xfclose(nf->fp);
  gt_str_delete(nf->outfilename);
}

static int outputsuflcpllv(void *processinfo,
                           const unsigned long *sequenceoffsettable,
                           const Suflcpbuffer *buf,
                           GtError *err)
{
  Mergeoutinfo *mergeoutinfo = (Mergeoutinfo *) processinfo;

  unsigned int i, lastindex;
  unsigned long lcpvalue;
  Largelcpvalue currentexception;
  GtUchar smallvalue;
  bool haserr = false;

  gt_error_check(err);
  for (i=0; i<buf->nextstoreidx; i++)
  {
    mergeoutinfo->absstartpostable[i]
      = sequenceoffsettable[buf->suftabstore[i].idx] +
        buf->suftabstore[i].startpos;
  }
  if (fwrite(mergeoutinfo->absstartpostable,
            sizeof (unsigned long),
            (size_t) buf->nextstoreidx,
            mergeoutinfo->outsuf.fp)
         != (size_t) buf->nextstoreidx)
  {
    gt_error_set(err,"fwrite(%s) of %u unsigned long-value failed: %s",
                  gt_str_get(mergeoutinfo->outsuf.outfilename),
                  buf->nextstoreidx,strerror(errno));
    haserr = true;
  }
  if (!haserr)
  {
    if (buf->lastpage)
    {
      lastindex = buf->nextstoreidx - 1;
    } else
    {
      lastindex = buf->nextstoreidx;
    }
    for (i=0; i<lastindex; i++)
    {
      lcpvalue = buf->lcptabstore[i];
      if (lcpvalue < (unsigned long) LCPOVERFLOW)
      {
        smallvalue = (GtUchar) lcpvalue;
      } else
      {
        currentexception.position = mergeoutinfo->currentlcpindex;
        currentexception.value = lcpvalue;
        if (fwrite(&currentexception,sizeof (Largelcpvalue),
                 (size_t) 1,mergeoutinfo->outllv.fp) != (size_t) 1)
        {
          gt_error_set(err,"fwrite(%s) of Largelcpvalue failed: %s",
                        gt_str_get(mergeoutinfo->outllv.outfilename),
                        strerror(errno));
          haserr = true;
          break;
        }
        smallvalue = LCPOVERFLOW;
      }
      if (fwrite(&smallvalue,sizeof (GtUchar),(size_t) 1,
                mergeoutinfo->outlcp.fp) != (size_t) 1)
      {
        gt_error_set(err,"fwrite(%s) of GtUchar failed: %s",
                       gt_str_get(mergeoutinfo->outlcp.outfilename),
                       strerror(errno));
        haserr = true;
        break;
      }
      mergeoutinfo->currentlcpindex++;
    }
  }
  return haserr ? -1 : 0;
}

static int mergeandstoreindex(const GtStr *storeindex,
                              Emissionmergedesa *emmesa,
                              GtError *err)
{
  Mergeoutinfo mergeoutinfo;
  GtUchar smalllcpvalue;
  GtSpecialcharinfo specialcharinfo;
  unsigned long *sequenceoffsettable, totallength;
  bool haserr = false;

  gt_error_check(err);
  if (initNameandFILE(&mergeoutinfo.outsuf,storeindex,SUFTABSUFFIX,err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (initNameandFILE(&mergeoutinfo.outlcp,storeindex,LCPTABSUFFIX,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (initNameandFILE(&mergeoutinfo.outllv,storeindex,LARGELCPTABSUFFIX,
                       err) != 0)
    {
      haserr = true;
    }
  }
  smalllcpvalue = 0;
  if (!haserr && fwrite(&smalllcpvalue,sizeof (GtUchar),(size_t) 1,
                mergeoutinfo.outlcp.fp) != (size_t) 1)
  {
    gt_error_set(err,"fwrite(%s) failed: %s",
                  gt_str_get(mergeoutinfo.outlcp.outfilename),
                  strerror(errno));
    haserr = true;
  }
  if (!haserr)
  {
    mergeoutinfo.currentlcpindex = (unsigned long) 1;
    sequenceoffsettable = gt_encseqtable2sequenceoffsets(&totallength,
                                                      &specialcharinfo,
                                                      emmesa->suffixarraytable,
                                                      emmesa->numofindexes);
    gt_assert(sequenceoffsettable != NULL);
    while (emmesa->numofentries > 0)
    {
      if (gt_emissionmergedesa_stepdeleteandinsertothersuffixes(emmesa,
                                                                err) != 0)
      {
        haserr = true;
        break;
      }
      if (outputsuflcpllv(&mergeoutinfo,
                         sequenceoffsettable,
                         &emmesa->buf,
                         err) != 0)
      {
        haserr = true;
        break;
      }
    }
    FREESPACE(sequenceoffsettable);
  }
  freeNameandFILE(&mergeoutinfo.outsuf);
  freeNameandFILE(&mergeoutinfo.outlcp);
  freeNameandFILE(&mergeoutinfo.outllv);
  return haserr ? -1 : 0;
}

int gt_performtheindexmerging(const GtStr *storeindex,
                           const GtStrArray *indexnametab,
                           GtLogger *logger,
                           GtError *err)
{
  Emissionmergedesa emmesa;
  unsigned int demand = SARR_ESQTAB | SARR_SUFTAB | SARR_LCPTAB;
  bool haserr = false;

  gt_error_check(err);
  if (gt_emissionmergedesa_init(&emmesa,
                             indexnametab,
                             demand,
                             logger,
                             err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_str_array_size(indexnametab) > 1UL)
    {
      if (mergeandstoreindex(storeindex,&emmesa,err) != 0)
      {
        haserr = true;
      }
    } else
    {
      gt_error_set(err,"merging requires more than one index");
      haserr = true;
    }
  }
  gt_emissionmergedesa_wrap(&emmesa);
  return haserr ? -1 : 0;
}
