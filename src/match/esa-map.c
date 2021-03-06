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

#include <stdio.h>
#include <errno.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include "core/array.h"
#include "core/endianess_api.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/format64.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "sarr-def.h"
#include "esa-fileend.h"
#include "sfx-ri-def.h"
#include "spacedef.h"
#include "bcktab.h"

#define DBFILEKEY "dbfile="

#define INITBufferedfile(INDEXNAME,STREAM,TYPE,SUFFIX)\
        (STREAM)->fp = gt_fa_fopen_with_suffix(INDEXNAME,SUFFIX,"rb",err);\
        if ((STREAM)->fp == NULL)\
        {\
          haserr = true;\
          (STREAM)->bufferedfilespace = NULL;\
        } else\
        {\
          (STREAM)->nextread = 0;\
          (STREAM)->nextfree = 0;\
          ALLOCASSIGNSPACE((STREAM)->bufferedfilespace,NULL,TYPE,\
                           FILEBUFFERSIZE);\
        }

static int scanprjfileuintkeysviafileptr(Suffixarray *suffixarray,
                                         const char *indexname,
                                         GtLogger *logger,
                                         FILE *fpin,
                                         GtError *err)
{
  uint32_t integersize, littleendian, readmodeint;
  unsigned int linenum;
  unsigned long currentlinelength;

  Definedunsignedlong maxbranchdepth;
  size_t dbfilelen = strlen(DBFILEKEY);
  bool haserr = false;
  GtArray *riktab;
  GtStr *currentline;
  /* the following five variables are local as the parsed values are
     not required: they are determined by reading the encseq */
  GtSpecialcharinfo specialcharinfo;
  unsigned long totallength,
                numofsequences,
                numofdbsequences,
                numofquerysequences;

  gt_error_check(err);
  riktab = gt_array_new(gt_sizeofReadintkeys());
  SETREADINTKEYS("totallength",&totallength,NULL);
  SETREADINTKEYS("specialcharacters",
                 &specialcharinfo.specialcharacters,NULL);
  SETREADINTKEYS("specialranges",
                 &specialcharinfo.specialranges,NULL);
  SETREADINTKEYS("realspecialranges",
                 &specialcharinfo.realspecialranges,NULL);
  SETREADINTKEYS("lengthofspecialprefix",
                 &specialcharinfo.lengthofspecialprefix,NULL);
  SETREADINTKEYS("lengthofspecialsuffix",
                 &specialcharinfo.lengthofspecialsuffix,NULL);
  SETREADINTKEYS("numofsequences",&numofsequences,NULL);
  SETREADINTKEYS("numofdbsequences",&numofdbsequences,NULL);
  gt_setreadintkeys(riktab,"numofquerysequences",&numofquerysequences,0,NULL);
  SETREADINTKEYS("longest",&suffixarray->longest.valueunsignedlong,
                           &suffixarray->longest.defined);
  SETREADINTKEYS("prefixlength",&suffixarray->prefixlength,NULL);
  SETREADINTKEYS("largelcpvalues",
                 &suffixarray->numoflargelcpvalues.valueunsignedlong,
                 &suffixarray->numoflargelcpvalues.defined);
  SETREADINTKEYS("maxbranchdepth",&maxbranchdepth.valueunsignedlong,
                 &maxbranchdepth.defined);
  SETREADINTKEYS("integersize",&integersize,NULL);
  SETREADINTKEYS("littleendian",&littleendian,NULL);
  SETREADINTKEYS("readmode",&readmodeint,NULL);
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, fpin) != EOF; linenum++)
  {
    currentlinelength = gt_str_length(currentline);
    if (dbfilelen <= (size_t) currentlinelength &&
       memcmp(DBFILEKEY,gt_str_get(currentline),dbfilelen) == 0)
    {
      /* Nothing */
    } else
    {
      if (gt_analyzeuintline(indexname,
                         PROJECTFILESUFFIX,
                         linenum,
                         gt_str_get(currentline),
                         currentlinelength,
                         riktab,
                         err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_str_reset(currentline);
  }
  gt_str_delete(currentline);
  if (!haserr && gt_allkeysdefined(indexname,PROJECTFILESUFFIX,riktab,
                                logger,err) != 0)
  {
    haserr = true;
  }
  if (!haserr &&
      integersize != (uint32_t) 32 &&
      integersize != (uint32_t) 64)
  {
    gt_error_set(err,"%s%s contains illegal line defining the integer size",
                 indexname,PROJECTFILESUFFIX);
    haserr = true;
  }
  if (!haserr && integersize != (uint32_t) (sizeof (unsigned long) * CHAR_BIT))
  {
    gt_error_set(err,"index was generated for %u-bit integers while "
                      "this program uses %u-bit integers",
                      (unsigned int) integersize,
                      (unsigned int) (sizeof (unsigned long) * CHAR_BIT));
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_is_little_endian())
    {
      if (littleendian != (uint32_t) 1)
      {
        gt_error_set(err,"computer has little endian byte order, while index "
                      "was build on computer with big endian byte order");
        haserr = true;
      }
    } else
    {
      if (littleendian == (uint32_t) 1)
      {
        gt_error_set(err,"computer has big endian byte order, while index "
                      "was build on computer with little endian byte "
                      "order");
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    if (readmodeint > (uint32_t) 3)
    {
      gt_error_set(err,"illegal readmode %u",(unsigned int) readmodeint);
      haserr = true;
    }
    suffixarray->readmode = (GtReadmode) readmodeint;
  }
  gt_array_delete(riktab);
  return haserr ? -1 : 0;
}

static void initsuffixarray(Suffixarray *suffixarray)
{
  suffixarray->encseq = NULL;
  suffixarray->suftab = NULL;
  suffixarray->lcptab = NULL;
  suffixarray->llvtab = NULL;
  suffixarray->bwttab = NULL;
  suffixarray->bcktab = NULL;
  suffixarray->bwttabstream.fp = NULL;
  suffixarray->suftabstream.fp = NULL;
  suffixarray->llvtabstream.fp = NULL;
  suffixarray->lcptabstream.fp = NULL;
  suffixarray->suftabstream.bufferedfilespace = NULL;
  suffixarray->lcptabstream.bufferedfilespace = NULL;
  suffixarray->llvtabstream.bufferedfilespace = NULL;
  suffixarray->bwttabstream.bufferedfilespace = NULL;
}

static bool scanprjfileuintkeys(Suffixarray *suffixarray,
                                const char *indexname,
                                GtLogger *logger,
                                GtError *err)
{
  bool haserr = false;
  FILE *fp;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,PROJECTFILESUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr && scanprjfileuintkeysviafileptr(suffixarray,
                                               indexname,logger,
                                               fp,err) != 0)
  {
    haserr = true;
  }
  gt_fa_xfclose(fp);
  return haserr;
}

void gt_freesuffixarray(Suffixarray *suffixarray)
{
  gt_fa_xmunmap((void *) suffixarray->suftab);
  suffixarray->suftab = NULL;
  gt_fa_xmunmap((void *) suffixarray->lcptab);
  suffixarray->lcptab = NULL;
  gt_fa_xmunmap((void *) suffixarray->llvtab);
  suffixarray->llvtab = NULL;
  gt_fa_xmunmap((void *) suffixarray->bwttab);
  suffixarray->bwttab = NULL;
  gt_fa_xfclose(suffixarray->suftabstream.fp);
  suffixarray->suftabstream.fp = NULL;
  FREESPACE(suffixarray->suftabstream.bufferedfilespace);
  gt_fa_xfclose(suffixarray->lcptabstream.fp);
  suffixarray->lcptabstream.fp = NULL;
  FREESPACE(suffixarray->lcptabstream.bufferedfilespace);
  gt_fa_xfclose(suffixarray->llvtabstream.fp);
  suffixarray->llvtabstream.fp = NULL;
  FREESPACE(suffixarray->llvtabstream.bufferedfilespace);
  gt_fa_xfclose(suffixarray->bwttabstream.fp);
  suffixarray->bwttabstream.fp = NULL;
  FREESPACE(suffixarray->bwttabstream.bufferedfilespace);
  gt_encseq_delete(suffixarray->encseq);
  suffixarray->encseq = NULL;
  if (suffixarray->bcktab != NULL)
  {
    gt_bcktab_delete(suffixarray->bcktab);
    suffixarray->bcktab = NULL;
  }
}

static int inputsuffixarray(bool map,
                            Suffixarray *suffixarray,
                            unsigned int demand,
                            const char *indexname,
                            GtLogger *logger,
                            GtError *err)
{
  bool haserr = false;
  GtEncseqLoader *el;
  unsigned long totallength = 0;

  gt_error_check(err);
  initsuffixarray(suffixarray);
  el = gt_encseq_loader_new();
  if (!(demand & SARR_DESTAB))
    gt_encseq_loader_do_not_require_des_tab(el);
  if (!(demand & SARR_SDSTAB))
    gt_encseq_loader_do_not_require_sds_tab(el);
  if (!(demand & SARR_SSPTAB))
    gt_encseq_loader_do_not_require_ssp_tab(el);
  gt_encseq_loader_set_logger(el, logger);
  suffixarray->encseq = gt_encseq_loader_load(el, indexname, err);
  gt_encseq_loader_delete(el);
  if (suffixarray->encseq == NULL)
  {
    haserr = true;
  } else
  {
    totallength = gt_encseq_total_length(suffixarray->encseq);
  }
  if (!haserr)
  {
    haserr = scanprjfileuintkeys(suffixarray,indexname,logger,err);
  }
  if (!haserr && (demand & SARR_SUFTAB))
  {
    if (map)
    {
      suffixarray->suftab
        = gt_mmap_check_size_with_suffix(indexname,
                                         SUFTABSUFFIX,
                                         totallength+1,
                                         sizeof (*suffixarray->suftab),
                                         err);
      if (suffixarray->suftab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->suftabstream,unsigned long,
                       SUFTABSUFFIX);
    }
    if (!haserr && !suffixarray->longest.defined)
    {
      gt_error_set(err,"longest not defined");
      haserr = true;
    }
  }
  if (!haserr && (demand & SARR_LCPTAB))
  {
    if (map)
    {
      suffixarray->lcptab = gt_mmap_check_size_with_suffix(indexname,
                                            LCPTABSUFFIX,
                                            totallength+1,
                                            sizeof (*suffixarray->lcptab),
                                            err);
      if (suffixarray->lcptab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->lcptabstream,GtUchar,
                       LCPTABSUFFIX);
      if (!haserr &&
          fseek(suffixarray->lcptabstream.fp,(long) sizeof (GtUchar),SEEK_SET))
      {
        gt_error_set(err,"fseek(esastream) failed: %s",strerror(errno));
        haserr = true;
      }
    }
    if (!haserr && !suffixarray->numoflargelcpvalues.defined)
    {
      gt_error_set(err,"numoflargelcpvalues not defined");
      haserr = true;
    }
    if (!haserr && suffixarray->numoflargelcpvalues.valueunsignedlong > 0)
    {
      if (map)
      {
        suffixarray->llvtab
          = gt_mmap_check_size_with_suffix(indexname,
                                           LARGELCPTABSUFFIX,
                                           (unsigned long)
                                           suffixarray->numoflargelcpvalues.
                                           valueunsignedlong,
                                           sizeof (*suffixarray->llvtab),
                                           err);
        if (suffixarray->llvtab == NULL)
        {
          haserr = true;
        }
      } else
      {
        INITBufferedfile(indexname,&suffixarray->llvtabstream,Largelcpvalue,
                         LARGELCPTABSUFFIX);
      }
    }
  }
  if (!haserr && (demand & SARR_BWTTAB))
  {
    if (map)
    {
      suffixarray->bwttab
        = gt_mmap_check_size_with_suffix(indexname,
                                         BWTTABSUFFIX,
                                         totallength+1,
                                         sizeof (*suffixarray->bwttab),
                                         err);
      if (suffixarray->bwttab == NULL)
      {
        haserr = true;
      }
    } else
    {
      INITBufferedfile(indexname,&suffixarray->bwttabstream,GtUchar,
                       BWTTABSUFFIX);
    }
  }
  if (!haserr && (demand & SARR_BCKTAB))
  {
    if (map)
    {
      suffixarray->bcktab = gt_mapbcktab(indexname,
                     gt_alphabet_num_of_chars(
                              gt_encseq_alphabet(suffixarray->encseq)),
                     suffixarray->prefixlength,
                     err);
      if (suffixarray->bcktab == NULL)
      {
        haserr = true;
      }
    } else
    {
      gt_error_set(err,"cannot stream bcktab");
      haserr = true;
    }
  }
  if (haserr)
  {
    gt_freesuffixarray(suffixarray);
  }
  return haserr ? -1 : 0;
}

int streamsuffixarray(Suffixarray *suffixarray,
                      unsigned int demand,
                      const char *indexname,
                      GtLogger *logger,
                      GtError *err)
{
  gt_error_check(err);
  return inputsuffixarray(false,
                          suffixarray,
                          demand,
                          indexname,
                          logger,
                          err);
}

int gt_mapsuffixarray(Suffixarray *suffixarray,
                      unsigned int demand,
                      const char *indexname,
                      GtLogger *logger,
                      GtError *err)
{
  gt_error_check(err);
  /* printf("sizeof (Suffixarray)=%lu\n",sizeof (Suffixarray)); */
  return inputsuffixarray(true,
                          suffixarray,
                          demand,
                          indexname,
                          logger,
                          err);
}
