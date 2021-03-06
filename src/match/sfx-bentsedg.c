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

#include <limits.h>
#include <stdio.h>
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/xansi_api.h"
#include "core/fa.h"
#include "core/arraydef.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/encseq.h"
#include "turnwheels.h"
#include "esa-fileend.h"
#include "bcktab.h"
#include "kmer2string.h"
#include "lcpoverflow.h"
#include "sfx-bltrie.h"
#include "sfx-remainsort.h"
#include "sfx-copysort.h"
#include "sfx-bentsedg.h"
#include "sfx-suffixgetset.h"
#include "stamp.h"

#define ACCESSCHARRAND(POS)    gt_encseq_get_encoded_char(bsr->encseq,\
                                                          POS,bsr->readmode)
#define ACCESSCHARSEQ(ESR)     gt_encseq_reader_next_encoded_char(ESR)
#define ISNOTEND(POS)          ((POS) < bsr->totallength &&\
                                ISNOTSPECIAL(ACCESSCHARRAND(POS)))

#define DEREFSTOPPOSSEQ(VAR,POS,STOPPOS,ESR)\
        (((POS) < (STOPPOS) && ISNOTSPECIAL(VAR = ACCESSCHARSEQ(ESR))) ?\
        ((unsigned long) VAR) : GT_UNIQUEINT(POS))

#define DEREFSEQ(VAR,POS,ESR) DEREFSTOPPOSSEQ(VAR,POS,bsr->totallength,ESR)

#define BS_SWAPARRAY(TMP,SUBBUCKETLEFT,IDX1,IDX2)\
        if ((IDX1) != (IDX2))\
        {\
          TMP = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX1);\
          gt_suffixsortspace_set(bsr->sssp,SUBBUCKETLEFT,IDX1,\
                                 gt_suffixsortspace_get(bsr->sssp,\
                                                        SUBBUCKETLEFT,\
                                                        IDX2));\
          gt_suffixsortspace_set(bsr->sssp,SUBBUCKETLEFT,IDX2,TMP);\
        }

#define STACKTOP\
        bsr->mkvauxstack.spaceMKVstack[bsr->mkvauxstack.nextfreeMKVstack]

#define UPDATELCP(MINVAL,MAXVAL)\
        gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);\
        if ((MINVAL) > commonunits.common)\
        {\
          MINVAL = commonunits.common;\
        }\
        if ((MAXVAL) < commonunits.common)\
        {\
          MAXVAL = commonunits.common;\
        }

GT_DECLAREARRAYSTRUCT(Largelcpvalue);

typedef struct
{
  void *reservoir;
  size_t sizereservoir;
  unsigned long *bucketoflcpvalues, /* pointer into reservoir */
                maxbranchdepth,
                numoflargelcpvalues,
                totalnumoflargelcpvalues,
                countoutputlcpvalues;
  uint8_t *smalllcpvalues; /* pointer into reservoir */
  const Compressedtable *completelcpvalues;
  GtArrayLargelcpvalue largelcpvalues;
} Lcpsubtab;

typedef struct
{
  bool defined;
  GtCodetype code;
  unsigned int prefixindex;
#undef SKDEBUG
#ifdef SKDEBUG
  unsigned long startpos;
#endif
} Suffixwithcode;

struct Outlcpinfo
{
  FILE *outfplcptab,
       *outfpllvtab;
  unsigned long totallength;
  Turningwheel *tw;
  unsigned int minchanged;
  Lcpsubtab lcpsubtab;
  Suffixwithcode previoussuffix;
  bool previousbucketwasempty,
       assideeffect;
};

#define CMPCHARBYCHARPTR2INT(VAR,SUBBUCKETLEFT,TMPVAR,IDX)\
        VAR = (((cptr = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX)+\
                        depth)\
                < bsr->totallength &&\
                ISNOTSPECIAL(TMPVAR = ACCESSCHARRAND(cptr)))\
                    ? ((unsigned long) TMPVAR) : GT_UNIQUEINT(cptr))

typedef GtEndofTwobitencoding Sfxcmp;

#define PTR2INT(VAR,SUBBUCKETLEFT,IDX)\
        {\
          unsigned long pos\
            = gt_suffixsortspace_get(bsr->sssp,SUBBUCKETLEFT,IDX);\
          if (pos + depth < bsr->totallength)\
          {\
            pos += depth;\
            gt_encseq_extract2bitencwithtwobitencodingstoppos(&(VAR),\
                                                              bsr->esr1,\
                                                              bsr->encseq,\
                                                              bsr->readmode,\
                                                              pos);\
          } else\
          {\
            VAR.tbe = 0;\
            VAR.unitsnotspecial = 0;\
            VAR.position = pos;\
          }\
        }

#define Sfxdocompare(COMMONUNITS,X,Y)\
        ret##X##Y = gt_encseq_compare_pairof_twobitencodings(bsr->fwd,\
                                                             bsr->complement,\
                                                             COMMONUNITS,&X,&Y)

#define SfxcmpEQUAL(X,Y)      (ret##X##Y == 0)
#define SfxcmpSMALLER(X,Y)    (ret##X##Y < 0)
#define SfxcmpGREATER(X,Y)    (ret##X##Y > 0)

typedef struct
{
  unsigned long subbucketleft,
                width,
                depth;
} MKVstack;

typedef struct
{
  GtEndofTwobitencoding etbe;
  unsigned long suftaboffset;
} Medianinfo;

typedef Medianinfo MedianElem;

typedef struct
{
  unsigned long suffix;
  unsigned char lcpwithpivot;
  char cmpresult;
} Countingsortinfo;

GT_DECLAREARRAYSTRUCT(MKVstack);

typedef struct
{
  const GtEncseq *encseq;
  GtEncseqReader *esr1, /* XXX be carefull with threads */
                 *esr2;
  GtReadmode readmode;
  bool fwd, complement, assideeffect;
  unsigned long totallength;
  GtArrayMKVstack mkvauxstack; /* XXX be carefull with treads */
  Lcpsubtab *lcpsubtab;
  Medianinfo *medianinfospace;
  Countingsortinfo *countingsortinfo;
  const Sfxstrategy *sfxstrategy;
  Blindtrie *blindtrie;
  Rmnsufinfo *rmnsufinfo;
  unsigned long leftlcpdist[GT_UNITSIN2BITENC],
                rightlcpdist[GT_UNITSIN2BITENC];
  GtSuffixsortspace *sssp;
  Dc_processunsortedrange dc_processunsortedrange;
  void *voiddcov;
  bool *equalwithprevious;
  unsigned long countinsertionsort,
                countqsort,
                countcountingsort,
                countbltriesort;
} Bentsedgresources;

#ifdef SKDEBUG
static unsigned long baseptr;

static void showsuffixrange(const Bentsedgresources *bsr,
                            unsigned long subbucketleft,
                            unsigned long width,
                            unsigned long depth)
{
  unsigned long pi;

  if (bsr->lcpsubtab == NULL)
  {
    printf("of %lu suffixes at depth %lu:\n",width,depth);
  } else
  {
    printf("of %lu suffixes [%lu,%lu] at depth %lu:\n",
           width,
           bsr->sssp->bucketleftidx + subbucketleft,
           bsr->sssp->bucketleftidx + subbucketleft + width,
           depth);
  }
  for (pi = 0; pi <= width; pi++)
  {
    unsigned long pos = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pi);
    printf("suffix %lu:",pos);
    gt_encseq_showatstartpos(stdout,bsr->fwd,bsr->complement,bsr->encseq,pos);
  }
}
#endif

#ifdef WITHCHECKSTARTPOINTER
static unsigned int checkstartpointorder(const unsigned long *left,
                                         const unsigned long *right)
{
  const unsigned long *ptr;
  bool ascending;

  gt_assert(left < right);
  gt_assert(*left != *(left+1));
  ascending = (*left < *(left+1)) ? true : false;
  for (ptr = left+1; ptr < right; ptr++)
  {
    gt_assert(*ptr != *(ptr+1));
    if (*ptr < *(ptr+1))
    {
      if (!ascending)
      {
        return 0;
      }
    } else
    {
      if (*ptr > *(ptr+1))
      {
        if (ascending)
        {
          return 0;
        }
      }
    }
  }
  return ascending ? 1U : 2U;
}
#endif

static unsigned long medianof3cmpcharbychar(const Bentsedgresources *bsr,
                                            unsigned long subbucketleft,
                                            unsigned long depth,
                                            unsigned long a,
                                            unsigned long b,
                                            unsigned long c)
{
  unsigned long vala, valb, valc, cptr;
  GtUchar tmpavar, tmpbvar;

  CMPCHARBYCHARPTR2INT(vala,subbucketleft,tmpavar,a);
  CMPCHARBYCHARPTR2INT(valb,subbucketleft,tmpbvar,b);
  if (vala == valb)
  {
    return a;
  }
  CMPCHARBYCHARPTR2INT(valc,subbucketleft,tmpavar,c);
  if (vala == valc || valb == valc)
  {
    return c;
  }
  return vala < valb ?
        (valb < valc ? b : (vala < valc ? c : a))
      : (valb > valc ? b : (vala < valc ? a : c));
}

static unsigned long medianof3(const Bentsedgresources *bsr,
                               unsigned long subbucketleft,
                               unsigned long depth,
                               unsigned long a,
                               unsigned long b,
                               unsigned long c)
{
  Sfxcmp vala, valb, valc;
  GtCommonunits commonunits;
  int retvalavalb, retvalavalc, retvalbvalc;

  PTR2INT(vala,subbucketleft,a);
  PTR2INT(valb,subbucketleft,b);
  Sfxdocompare(&commonunits,vala,valb);
  if (SfxcmpEQUAL(vala,valb))
  {
    return a;
  }
  PTR2INT(valc,subbucketleft,c);
  Sfxdocompare(&commonunits,vala,valc);
  if (SfxcmpEQUAL(vala,valc))
  {
    return c;
  }
  Sfxdocompare(&commonunits,valb,valc);
  if (SfxcmpEQUAL(valb,valc))
  {
    return c;
  }
  return SfxcmpSMALLER(vala,valb) ?
        (SfxcmpSMALLER(valb,valc) ? b : (SfxcmpSMALLER(vala,valc) ? c : a))
      : (SfxcmpGREATER(valb,valc) ? b : (SfxcmpSMALLER(vala,valc) ? a : c));
}

static void updatelcpvalue(Bentsedgresources *bsr,
                           unsigned long idx,
                           unsigned long value)
{
  if (value >= (unsigned long) LCPOVERFLOW)
  {
    bsr->lcpsubtab->numoflargelcpvalues++; /* this may overcount as
                                              there may be some value
                                              which was already overflowing */
  }
  bsr->lcpsubtab->bucketoflcpvalues[idx] = value;
}

static void bs_insertionsort(Bentsedgresources *bsr,
                             unsigned long subbucketleft,
                             unsigned long width,
                             unsigned long offset)
{
  unsigned long sval1, sval2, pi, pj, startpos1, startpos2, temp,
                lcpindex, lcplen = 0;
  int retval;
  GtCommonunits commonunits;

#ifdef SKDEBUG
  printf("insertion sort ");
  showsuffixrange(bsr,subbucketleft,width,offset);
#endif
  bsr->countinsertionsort++;
  for (pi = 1UL; pi < width; pi++)
  {
    for (pj = pi; pj > 0; pj--)
    {
      sval1 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pj-1);
      sval2 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pj);
      if (bsr->sfxstrategy->cmpcharbychar)
      {
        startpos1 = sval1 + offset;
        if (startpos1 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr1, bsr->encseq,
                                                bsr->readmode, startpos1);
        }
        startpos2 = sval2 + offset;
        if (startpos2 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr2, bsr->encseq,
                                                bsr->readmode, startpos2);
        }
        for (;;)
        {
          unsigned long ccs, cct;
          GtUchar tmp1, tmp2;

          ccs = DEREFSEQ(tmp1,startpos1,bsr->esr1);
          cct = DEREFSEQ(tmp2,startpos2,bsr->esr2);
          if (ccs != cct)
          {
            lcplen = startpos2 - sval2;
            retval = (ccs < cct) ? -1 : 1;
            break;
          }
          startpos1++;
          startpos2++;
        }
      } else
      {
#ifdef SKDEBUG
        printf("gt_encseq_compare_viatwobitencoding[%lu,%lu] "
               "at offset %lu\n",
                sval1,sval2,offset);
        gt_encseq_showatstartpos(stdout,
                                 bsr->fwd,
                                 bsr->complement,
                                 bsr->encseq,
                                 sval1);
        gt_encseq_showatstartpos(stdout,
                                 bsr->fwd,
                                 bsr->complement,
                                 bsr->encseq,
                                 sval2);
#endif
        retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     bsr->encseq,
                                                     bsr->readmode,
                                                     bsr->esr1,
                                                     bsr->esr2,
                                                     sval1,
                                                     sval2,
                                                     offset,
                                                     0);
        lcplen = commonunits.finaldepth;
      }
      gt_assert(retval != 0);
      if (bsr->lcpsubtab != NULL && bsr->assideeffect)
      {
        lcpindex = subbucketleft+pj;
        if (pj < pi && retval > 0)
        {
          updatelcpvalue(bsr,lcpindex+1,
                         bsr->lcpsubtab->bucketoflcpvalues[lcpindex]);
        }
        updatelcpvalue(bsr,lcpindex,lcplen);
      }
      if (retval < 0)
      {
        break;
      }
      BS_SWAPARRAY(temp,subbucketleft,pj,pj-1);
    }
  }
}

static void bs_insertionsortmaxdepth(Bentsedgresources *bsr,
                                     unsigned long subbucketleft,
                                     unsigned long width,
                                     unsigned long offset,
                                     unsigned long maxdepth)
{
  unsigned long sval1, sval2, pi, pj, startpos1, startpos2, temp,
                lcpindex, lcplen = 0, idx = 0;
  int retval;
  bool tempb;
  GtCommonunits commonunits;

#ifdef SKDEBUG
  printf("insertion sort (offset=%lu,maxdepth=%lu)\n",offset,maxdepth);
  showsuffixrange(bsr,subbucketleft,width,offset);
#endif
  bsr->countinsertionsort++;
  for (pi = 1UL; pi < width; pi++)
  {
    for (pj = pi; pj > 0; pj--)
    {
      sval1 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pj-1);
      sval2 = gt_suffixsortspace_get(bsr->sssp,subbucketleft,pj);
      if (bsr->sfxstrategy->cmpcharbychar)
      {
        unsigned long endpos1, endpos2;

        endpos1 = sval1+maxdepth;
        if (endpos1 > bsr->totallength)
        {
          endpos1 = bsr->totallength;
        }
        endpos2 = sval2+maxdepth;
        if (endpos2 > bsr->totallength)
        {
          endpos2 = bsr->totallength;
        }
        startpos1 = sval1+offset;
        if (startpos1 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr1, bsr->encseq,
                                                bsr->readmode, startpos1);
        }
        startpos2 = sval2+offset;
        if (startpos2 < bsr->totallength)
        {
          gt_encseq_reader_reinit_with_readmode(bsr->esr2, bsr->encseq,
                                                bsr->readmode, startpos2);
        }
        for (;;)
        {
          unsigned long ccs, cct;
          GtUchar tmp1, tmp2;

          ccs = DEREFSTOPPOSSEQ(tmp1,startpos1,endpos1,bsr->esr1);
          cct = DEREFSTOPPOSSEQ(tmp2,startpos2,endpos2,bsr->esr2);
          lcplen = startpos2 - sval2;
          if (lcplen == maxdepth)
          {
            retval = 0;
            break;
          }
          gt_assert(lcplen < maxdepth);
          if (ccs != cct)
          {
            retval = (ccs < cct) ? -1 : 1;
            break;
          }
          startpos1++;
          startpos2++;
        }
      } else
      {
        gt_assert(offset < maxdepth);
        retval = gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     bsr->encseq,
                                                     bsr->readmode,
                                                     bsr->esr1,
                                                     bsr->esr2,
                                                     sval1,
                                                     sval2,
                                                     offset,
                                                     maxdepth);
        lcplen = commonunits.finaldepth;
        gt_assert(lcplen <= maxdepth);
        if (lcplen == maxdepth)
        {
          gt_assert(retval == 0);
        }
      }
#ifdef SKDEBUG
      printf("cmp %lu and %lu: retval = %d, lcplen = %lu\n",
             sval1, sval2, retval, (unsigned long) lcplen);
#endif
      if (retval != 0 && bsr->lcpsubtab != NULL && bsr->assideeffect)
      {
        lcpindex = subbucketleft + pj;
        if (pj < pi && retval > 0)
        {
          updatelcpvalue(bsr,lcpindex+1,
                         bsr->lcpsubtab->bucketoflcpvalues[lcpindex]);
        }
        updatelcpvalue(bsr,lcpindex,lcplen);
      }
      if (retval < 0)
      {
        break;
      }
      idx = pj;
      if (retval == 0)
      {
        gt_assert(idx > 0);
        bsr->equalwithprevious[idx] = true;
        break;
      }
      BS_SWAPARRAY(temp,subbucketleft,pj,pj-1);
      tempb = bsr->equalwithprevious[idx-1];
      bsr->equalwithprevious[idx-1] = bsr->equalwithprevious[idx];
      bsr->equalwithprevious[idx] = tempb;
    }
  }
  if (idx > 0)
  {
    unsigned long equalsrangewidth = 0;
    unsigned long bucketleftidx
     = gt_suffixsortspace_bucketleftidx_get(bsr->sssp);
#ifdef SKDEBUG
    printf("ordered suffix %lu\n",gt_suffixsortspace_get(bsr->sssp,
                                                         subbucketleft,0));
#endif
    for (idx = 1UL; idx < width; idx++)
    {
#ifdef SKDEBUG
      printf("ordered suffix %lu, equalwithprevious=%s\n",
              gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx),
              bsr->equalwithprevious[idx] ? "true" : "false");
#endif
      if (bsr->equalwithprevious[idx])
      {
        bsr->equalwithprevious[idx] = false;
        equalsrangewidth++;
      } else
      {
        if (equalsrangewidth > 0)
        {
#ifdef SKDEBUG
          printf("process interval of width %lu\n",
                 equalsrangewidth + 1);
#endif
          bsr->dc_processunsortedrange(
                              bsr->voiddcov,
                              bucketleftidx + subbucketleft + idx - 1
                                            - equalsrangewidth,
                              equalsrangewidth + 1, maxdepth);
          equalsrangewidth = 0;
        }
      }
    }
    if (equalsrangewidth > 0)
    {
#ifdef SKDEBUG
      printf("process interval of width %lu\n",
             equalsrangewidth + 1);
#endif
      bsr->dc_processunsortedrange(
                           bsr->voiddcov,
                           bucketleftidx + subbucketleft + width - 1
                                         - equalsrangewidth,
                           equalsrangewidth + 1, maxdepth);
    }
  }
}

#define DOMEDIANCOMPARE(A,B)\
        gt_encseq_compare_pairof_twobitencodings(fwd,complement,&commonunits,\
                                                 &((A)->etbe),&((B)->etbe))

#define MedianElemGREATER(A,B)  (DOMEDIANCOMPARE(A,B) > 0)

#define MedianElemSWAP(A,B)     {\
                                  register MedianElem tmp = *(A);\
                                                      *(A) = *(B);\
                                                      *(B) = tmp;\
                                }

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

static MedianElem *quickmedian (bool fwd,bool complement,
                                MedianElem *arr,unsigned long width)
{
  MedianElem *low, *high, *median, *middle, *ll, *hh;
  GtCommonunits commonunits;

  gt_assert(width > 0);
  low = arr;
  high = arr + width - 1;
  median = low + GT_DIV2(width);
  for (;;)
  {
    if (high <= low)                   /* One element only */
    {
      return median;
    }
    if (high == low + 1)
    {                                  /* Two elements only */
      if (MedianElemGREATER(low,high))
      {
        MedianElemSWAP (low, high);
      }
      return median;
    }

    /* Find median of low, middle and high items; swap into position low */
    middle = low + GT_DIV2(high - low + 1);
    if (MedianElemGREATER(middle,high))
    {
      MedianElemSWAP (middle, high);
    }
    if (MedianElemGREATER(low,high))
    {
      MedianElemSWAP (low, high);
    }
    if (MedianElemGREATER(middle,low))
    {
      MedianElemSWAP (middle, low);
    }
    /* Swap low item (now in position middle) into position (low+1) */
    MedianElemSWAP (middle, low + 1);

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;)
    {
      do
      {
        ll++;
      } while (MedianElemGREATER(low,ll));
      do
      {
        hh--;
      } while  (MedianElemGREATER(hh,low));
      if (hh < ll)
      {
        break;
      }
      MedianElemSWAP (ll, hh);
    }

    /* Swap middle item (in position low) back into correct position */
    MedianElemSWAP (low, hh);

    /* Re-set active partition */
    if (hh <= median)
    {
      low = ll;
    }
    if (hh >= median)
    {
      high = hh - 1;
    }
  }
}

#ifdef WITHcheckmedian

static void checkmedian(bool fwd,
                        bool complement,
                        const Medianinfo *median,
                        const Medianinfo *space,
                        unsigned long width)
{
  unsigned long sum1, sum2, idx, smaller = 0, larger = 0, equal = 0, equalpart;
  unsigned int commonunits;
  int cmp;

  for (idx = 0; idx < width; idx++)
  {
    cmp = DOMEDIANCOMPARE(space + idx,median);
    if (cmp > 0)
    {
     larger++;
    } else
    {
      if (cmp < 0)
      {
        smaller++;
      } else
      {
        equal++;
      }
    }
  }
  if (smaller == larger)
  {
    return;
  }
  for (equalpart = 0; equalpart < equal; equalpart++)
  {
    sum1 = smaller + equalpart;
    sum2 = larger + (equal-1) - equalpart;
    if (sum1 < sum2)
    {
      if (sum1 + 1 == sum2)
      {
        return;
      }
    } else
    {
      if (sum1 > sum2)
      {
        if (sum1 == sum2 + 1)
        {
          return;
        }
      } else
      {
        return;
      }
    }
  }
  fprintf(stderr,"problem with equal=%lu,smaller=%lu,larger=%lu\n",
                  equal,smaller,larger);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}
#endif

static unsigned long realmedian(const Bentsedgresources *bsr,
                                unsigned long subbucketleft,
                                unsigned long width,
                                unsigned long depth)
{
  Medianinfo *medianptr;
  unsigned long idx;

  for (idx = 0; idx < width; idx++)
  {
    bsr->medianinfospace[idx].suftaboffset = idx;
    PTR2INT(bsr->medianinfospace[idx].etbe,subbucketleft,idx);
  }
  medianptr = quickmedian(bsr->fwd,bsr->complement,bsr->medianinfospace,width);
/*
  checkmedian(bsr->fwd,bsr->complement,medianptr,medianinfospace,width);
*/
  gt_assert(medianptr != NULL);
  return medianptr->suftaboffset;
}

#define MINMEDIANOF9WIDTH 31UL

static unsigned long cmpcharbychardelivermedian(const Bentsedgresources *bsr,
                                                unsigned long subbucketleft,
                                                unsigned long width,
                                                unsigned long depth)
{
  unsigned long pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  { /* On big arrays, pseudomedian of 9 */
    unsigned long offset, doubleoffset;
    offset = GT_DIV8(width);
    doubleoffset = GT_MULT2(offset);
    pl = medianof3cmpcharbychar(bsr,subbucketleft,depth,pl,pl+offset,
                                pl+doubleoffset);
    pm = medianof3cmpcharbychar(bsr,subbucketleft,depth,pm-offset,
                                pm,pm+offset);
    pr = medianof3cmpcharbychar(bsr,subbucketleft,depth,
                                pr-doubleoffset,pr-offset,
                                pr);
  }
  return medianof3cmpcharbychar(bsr,subbucketleft,depth,pl,pm,pr);
}

static unsigned long blockcmpdelivermedian(const Bentsedgresources *bsr,
                                           unsigned long subbucketleft,
                                           unsigned long width,
                                           unsigned long depth,
                                           unsigned long maxwidthrealmedian)
{
  unsigned long pl = 0,
                pm = GT_DIV2(width),
                pr = width - 1;

  if (width >= MINMEDIANOF9WIDTH)
  {
    if (width > maxwidthrealmedian)
    { /* On big arrays, pseudomedian of 9 */
      unsigned long offset, doubleoffset;
      offset = GT_DIV8(width);
      doubleoffset = GT_MULT2(offset);
      pl = medianof3(bsr,subbucketleft,depth,pl,pl+offset,
                     pl+doubleoffset);
      pm = medianof3(bsr,subbucketleft,depth,pm-offset,pm,pm+offset);
      pr = medianof3(bsr,subbucketleft,depth,pr-doubleoffset,
                     pr-offset,pr);
      pm = medianof3(bsr,subbucketleft,depth,pl,pm,pr);
    } else /* width <= maxwidthrealmedian */
    {
      pm = realmedian(bsr, subbucketleft,width, depth);
    }
  } else
  {
    pm = medianof3(bsr,subbucketleft,depth,pl,pm,pr);
  }
  return pm;
}

/*
static void showcountingsortinfo(const Countingsortinfo *countingsortinfo,
                              unsigned long idx)
{
  printf("countingsortinfo[%lu]=(%lu,",idx,
          (unsigned long) countingsortinfo[idx].suffix);
  printf("%lu,",(unsigned long) countingsortinfo[idx].lcpwithpivot);
  printf("%d)\n",countingsortinfo[idx].cmpresult);
}
*/

static bool comparisonsort(Bentsedgresources *bsr,
                           unsigned long subbucketleft,
                           unsigned long width,
                           unsigned long depth)
{
  gt_assert(width > 1UL);
  gt_assert(bsr->sfxstrategy->maxinsertionsort <=
            bsr->sfxstrategy->maxbltriesort);
  if (width <= bsr->sfxstrategy->maxinsertionsort)
  {
    bs_insertionsort(bsr,subbucketleft,width,depth);
    return true;
  }
  if (width <= bsr->sfxstrategy->maxbltriesort)
  {
    unsigned long numoflargelcpvalues;

    gt_assert(bsr->sfxstrategy->differencecover == 0);
    numoflargelcpvalues
      = gt_blindtrie_suffixsort(bsr->blindtrie,
                                subbucketleft,
                                bsr->lcpsubtab == NULL
                                  ? NULL
                                  : bsr->lcpsubtab->bucketoflcpvalues +
                                    subbucketleft,
                                width,
                                depth,
                                (unsigned long)
                                   bsr->sfxstrategy->differencecover,
                                NULL,
                                NULL);
    if (bsr->lcpsubtab != NULL)
    {
      bsr->lcpsubtab->numoflargelcpvalues += numoflargelcpvalues;
    }
    bsr->countbltriesort++;
    return true;
  }
  return false;
}

static void subsort_bentleysedgewick(Bentsedgresources *bsr,
                                     unsigned long subbucketleft,
                                     unsigned long width,
                                     unsigned long depth)
{
  if (width > 1UL)
  {
    if (bsr->sfxstrategy->ssortmaxdepth.defined)
    {
      if (depth >=
               (unsigned long) bsr->sfxstrategy->ssortmaxdepth.valueunsignedint)
      {
        unsigned long leftindex
          = gt_suffixsortspace_bucketleftidx_get(bsr->sssp) + subbucketleft;
        gt_rmnsufinfo_addunsortedrange(bsr->rmnsufinfo,
                                       leftindex,
                                       leftindex + width - 1,
                                       depth);
        return;
      }
    } else
    {
      if (bsr->sfxstrategy->differencecover > 0)
      {
        if (depth >= (unsigned long) bsr->sfxstrategy->differencecover)
        {
          bsr->dc_processunsortedrange(bsr->voiddcov,subbucketleft,width,depth);
          return;
        }
        if (width <= bsr->sfxstrategy->maxinsertionsort)
        {
          bs_insertionsortmaxdepth(bsr,subbucketleft,width,depth,
                                   (unsigned long)
                                   bsr->sfxstrategy->differencecover);
          return;
        }
        if (width <= bsr->sfxstrategy->maxbltriesort)
        {
          unsigned long numoflargelcpvalues;

          numoflargelcpvalues
            = gt_blindtrie_suffixsort(bsr->blindtrie,
                                      subbucketleft,
                                      bsr->lcpsubtab == NULL
                                        ? NULL
                                        : bsr->lcpsubtab->bucketoflcpvalues +
                                          subbucketleft,
                                      width,
                                      depth,
                                      (unsigned long)
                                        bsr->sfxstrategy->differencecover,
                                      bsr->voiddcov,
                                      bsr->dc_processunsortedrange);
          if (bsr->lcpsubtab != NULL)
          {
            bsr->lcpsubtab->numoflargelcpvalues += numoflargelcpvalues;
          }
          bsr->countbltriesort++;
          return;
        }
      } else
      {
        if (comparisonsort(bsr,subbucketleft,width,depth))
        {
          return;
        }
      }
    }
    /* push */
    GT_CHECKARRAYSPACE(&bsr->mkvauxstack,MKVstack,1024);
    STACKTOP.subbucketleft = subbucketleft;
    STACKTOP.width = width;
    STACKTOP.depth = depth;
    bsr->mkvauxstack.nextfreeMKVstack++;
  }
}

static void sarrcountingsort(Bentsedgresources *bsr,
                             unsigned long subbucketleft,
                             unsigned long width,
                             const Sfxcmp *pivotcmpbits,
                             unsigned long pivotidx,
                             unsigned long depth)
{
  int cmp;
  unsigned int maxsmallerwithlcp = 0, maxlargerwithlcp = 0;
  GtCommonunits commonunits;
  GtEndofTwobitencoding etbecurrent;
  unsigned long idx, smaller = 0, larger = 0,
                insertindex, end, equaloffset, currentwidth;
  Countingsortinfo *csiptr;
  /* const bool cmpcharbychar = false; */

  bsr->countcountingsort++;
  for (idx = 0; idx < width; idx++)
  {
    if (idx != pivotidx)
    {
      PTR2INT(etbecurrent,subbucketleft,idx);
      cmp = gt_encseq_compare_pairof_twobitencodings(bsr->fwd,
                                                     bsr->complement,
                                                     &commonunits,
                                                     &etbecurrent,
                                                     pivotcmpbits);
      bsr->countingsortinfo[idx].suffix
        = gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx);
      gt_assert(commonunits.common <= (unsigned int) GT_UNITSIN2BITENC);
      bsr->countingsortinfo[idx].lcpwithpivot = commonunits.common;
      if (cmp > 0)
      {
        gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);
        bsr->rightlcpdist[commonunits.common]++;
        if (maxlargerwithlcp < commonunits.common)
        {
          maxlargerwithlcp = commonunits.common;
        }
        bsr->countingsortinfo[idx].cmpresult = (char) 1;
        larger++;
      } else
      {
        if (cmp < 0)
        {
          gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);
          bsr->leftlcpdist[commonunits.common]++;
          if (maxsmallerwithlcp < commonunits.common)
          {
            maxsmallerwithlcp = commonunits.common;
          }
          bsr->countingsortinfo[idx].cmpresult = (char) -1;
          smaller++;
        } else
        {
          gt_assert(commonunits.common == (unsigned int) GT_UNITSIN2BITENC);
          bsr->countingsortinfo[idx].cmpresult = 0;
        }
      }
    } else
    {
      bsr->countingsortinfo[idx].suffix
        = gt_suffixsortspace_get(bsr->sssp,subbucketleft,idx);
      bsr->countingsortinfo[idx].lcpwithpivot = (unsigned char)
                                                GT_UNITSIN2BITENC;
      bsr->countingsortinfo[idx].cmpresult = (char) 0;
    }
  }
  for (idx = 1UL; idx <= (unsigned long) maxsmallerwithlcp; idx++)
  {
    bsr->leftlcpdist[idx] += bsr->leftlcpdist[idx-1];
  }
  for (idx = 1UL; idx <= (unsigned long) maxlargerwithlcp; idx++)
  {
    bsr->rightlcpdist[idx] += bsr->rightlcpdist[idx-1];
  }
  equaloffset = width - larger;
  for (csiptr = bsr->countingsortinfo + width -1;
       csiptr >= bsr->countingsortinfo;
       csiptr--)
  {
    switch (csiptr->cmpresult)
    {
      case -1:
        insertindex = --(bsr->leftlcpdist[csiptr->lcpwithpivot]);
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,insertindex,
                               csiptr->suffix);
        break;
      case 0:
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,--equaloffset,
                               csiptr->suffix);
        break;
      case 1:
        insertindex = --(bsr->rightlcpdist[csiptr->lcpwithpivot]);
        gt_suffixsortspace_set(bsr->sssp,subbucketleft,
                               width - 1 - insertindex,
                               csiptr->suffix);
        break;
    }
  }
  for (idx = 0; idx <= (unsigned long) maxsmallerwithlcp; idx++)
  {
    if (idx < (unsigned long) maxsmallerwithlcp)
    {
      end = bsr->leftlcpdist[idx+1];
    } else
    {
      end = smaller;
    }
    if (bsr->leftlcpdist[idx] + 1 < end) /* at least two elements */
    {
      currentwidth = end - bsr->leftlcpdist[idx];
      subsort_bentleysedgewick(bsr,
                               subbucketleft + bsr->leftlcpdist[idx],
                               currentwidth,
                               depth + idx);
    }
    if (bsr->lcpsubtab != NULL && bsr->assideeffect &&
        bsr->leftlcpdist[idx] < end)
    { /* at least one element */
      updatelcpvalue(bsr,subbucketleft+end,depth + idx);
    }
    bsr->leftlcpdist[idx] = 0;
  }
  if (width - smaller - larger > 1UL)
  {
    currentwidth = width - smaller - larger;
    subsort_bentleysedgewick(bsr,
                             subbucketleft + smaller,
                             currentwidth,
                             depth + GT_UNITSIN2BITENC);
  }
  for (idx = 0; idx <= (unsigned long) maxlargerwithlcp; idx++)
  {
    if (idx < (unsigned long) maxlargerwithlcp)
    {
      end = bsr->rightlcpdist[idx+1];
    } else
    {
      end = larger;
    }
    if (bsr->rightlcpdist[idx] + 1 < end) /* at least two elements */
    {
      currentwidth = end - bsr->rightlcpdist[idx];
      subsort_bentleysedgewick(bsr,
                               subbucketleft + width - end,
                               currentwidth,
                               depth + idx);
    }
    if (bsr->lcpsubtab != NULL && bsr->assideeffect &&
        bsr->rightlcpdist[idx] < end)
    { /* at least one element */
      updatelcpvalue(bsr,subbucketleft + width - end,depth + idx);
    }
    bsr->rightlcpdist[idx] = 0;
  }
}

static inline void vectorswap(GtSuffixsortspace *sssp,
                              unsigned long subbucketleft1,
                              unsigned long subbucketleft2,
                              unsigned long width)
{
  unsigned long idx, tmp;

  for (idx = 0; idx < width; idx++)
  {
    tmp = gt_suffixsortspace_get(sssp,subbucketleft1,idx);
    gt_suffixsortspace_set(sssp,subbucketleft1,idx,
                           gt_suffixsortspace_get(sssp,subbucketleft2,idx));
    gt_suffixsortspace_set(sssp,subbucketleft2,idx,tmp);
  }
}

static void bentleysedgewick(Bentsedgresources *bsr,
                             unsigned long width,
                             unsigned long depth)
{
  bsr->mkvauxstack.nextfreeMKVstack = 0;
  subsort_bentleysedgewick(bsr, 0, width, depth);
  while (bsr->mkvauxstack.nextfreeMKVstack > 0)
  {
    unsigned long leftplusw, pa, pb, pc, pd, pm, bucketright,
                  cptr, temp, pivotcmpcharbychar = 0, valcmpcharbychar,
                  wtmp, subbucketleft;
    unsigned int smallermaxlcp, greatermaxlcp, smallerminlcp, greaterminlcp;
    Sfxcmp pivotcmpbits, val;
    int retvalpivotcmpbits;
    GtUchar tmpvar;
    GtCommonunits commonunits;
    const int commonunitsequal = bsr->sfxstrategy->cmpcharbychar
                                 ? 1
                                 : GT_UNITSIN2BITENC;

    /* pop */
    bsr->mkvauxstack.nextfreeMKVstack--;
    subbucketleft = STACKTOP.subbucketleft;
    width = STACKTOP.width;
    depth = STACKTOP.depth;
    bucketright = width - 1;

    if (bsr->sfxstrategy->cmpcharbychar)
    {
      pm = cmpcharbychardelivermedian(bsr, subbucketleft, width, depth);
      BS_SWAPARRAY(temp, subbucketleft, 0, pm);
      CMPCHARBYCHARPTR2INT(pivotcmpcharbychar, subbucketleft,tmpvar,0);
    } else
    {
      pm = blockcmpdelivermedian(bsr,
                                 subbucketleft,
                                 width,
                                 depth,
                                 bsr->sfxstrategy->maxwidthrealmedian);
      if (width <= bsr->sfxstrategy->maxcountingsort &&
          width >= MINMEDIANOF9WIDTH)
      {
        PTR2INT(pivotcmpbits,subbucketleft,pm);
        sarrcountingsort(bsr,
                         subbucketleft,
                         width,
                         &pivotcmpbits,
                         pm,
                         depth);
        /* new values for subbucketleft, bucketright, depth */
        continue;
      }
      BS_SWAPARRAY(temp, subbucketleft, 0, pm);
      PTR2INT(pivotcmpbits,subbucketleft,0);
    }
    bsr->countqsort++;
    /* now pivot element is at index subbucketleft */
    /* all elements to be compared are between pb and pc */
    /* pa is the position at which the next element smaller than the
       pivot element is inserted at */
    /* pd is the position at which the next element greater than the
       pivot element is inserted at */
    pa = pb = 1UL;
    pc = pd = bucketright;
    if (bsr->sfxstrategy->cmpcharbychar)
    {
      smallerminlcp = greaterminlcp = smallermaxlcp = greatermaxlcp = 0;
      for (;;)
      {
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,subbucketleft,tmpvar,pb);
          if (valcmpcharbychar > pivotcmpcharbychar)
          {
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            BS_SWAPARRAY(temp, subbucketleft, pa, pb);
            pa++;
          }
          pb++;
        }
        while (pb <= pc)
        {
          CMPCHARBYCHARPTR2INT(valcmpcharbychar,subbucketleft,tmpvar,pc);
          if (valcmpcharbychar < pivotcmpcharbychar)
          { /* stop for elements < pivot */
            break;
          }
          if (valcmpcharbychar == pivotcmpcharbychar)
          {
            /* exchange equal element and element at index pd */
            BS_SWAPARRAY(temp, subbucketleft, pc, pd);
            pd--;
          }
          pc--;
        }
        if (pb > pc)
        { /* no elements to compare to pivot */
          break;
        }
        BS_SWAPARRAY(temp, subbucketleft, pb, pc);
        pb++;
        pc--;
      }
    } else
    {
      smallermaxlcp = greatermaxlcp = 0;
      smallerminlcp = greaterminlcp = (unsigned int) GT_UNITSIN2BITENC;
      for (;;)
      {
        /* look for elements identical or smaller than pivot from left */
        while (pb <= pc)
        {
          PTR2INT(val,subbucketleft,pb);
          Sfxdocompare(&commonunits,val,pivotcmpbits);
          if (SfxcmpGREATER(val,pivotcmpbits))
          { /* stop for elements val > pivot */
            UPDATELCP(greaterminlcp,greatermaxlcp);
            break;
          }
          if (SfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pa, pb);
            pa++;
          } else /* smaller */
          {
            UPDATELCP(smallerminlcp,smallermaxlcp);
          }
          pb++;
        }
        /* look for elements identical or greater than pivot from right */
        while (pb <= pc)
        {
          PTR2INT(val,subbucketleft,pc);
          Sfxdocompare(&commonunits,val,pivotcmpbits);
          if (SfxcmpSMALLER(val,pivotcmpbits))
          { /* stop for elements val < pivot */
            UPDATELCP(smallerminlcp,smallermaxlcp);
            break;
          }
          if (SfxcmpEQUAL(val,pivotcmpbits))
          {
            /* exchange equal element and element at index pa */
            BS_SWAPARRAY(temp, subbucketleft, pc, pd);
            pd--;
          } else /* greater */
          {
            UPDATELCP(greaterminlcp,greatermaxlcp);
          }
          pc--;
        }
        if (pb > pc)
        { /* interval is empty */
          break;
        }
        BS_SWAPARRAY(temp, subbucketleft, pb, pc);
        pb++;
        pc--;
      }
    }
    gt_assert(pb >= pa);
    wtmp = MIN(pa,pb-pa);
    /* move w elements at the left to the middle */
    vectorswap(bsr->sssp, subbucketleft, subbucketleft+pb-wtmp, wtmp);
    gt_assert(pd >= pc);
    gt_assert(bucketright >= pd);
    wtmp = MIN(pd-pc, bucketright-pd);
    /* move w elements at the right to the middle */
    vectorswap(bsr->sssp, subbucketleft+pb, subbucketleft+bucketright+1-wtmp,
               wtmp);

    /* all elements equal to the pivot are now in the middle namely in the
       range [subbucketleft + (pb-pa) and bucketright - (pd-pc)] */
    /* hence we have to sort the elements in the intervals
       [subbucketleft..subbucketleft+(pb-pa)-1] and
       [bucketright-(pd-pc)+1..bucketright] */

    gt_assert(pb >= pa);
    if ((wtmp = pb-pa) > 0)
    {
      leftplusw = wtmp;
      if (bsr->lcpsubtab != NULL && bsr->assideeffect)
      {
        /*
          left part has suffix with lcp up to length smallermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the left
          which is at a minimum distance to the pivot and thus to an
          element in the final part of the left side.
        */
        updatelcpvalue(bsr,subbucketleft + leftplusw,depth + smallermaxlcp);
      }
      subsort_bentleysedgewick(bsr,
                               subbucketleft,
                               wtmp,
                               depth + smallerminlcp);
    } else
    {
      leftplusw = 0;
    }

    cptr = gt_suffixsortspace_get(bsr->sssp,subbucketleft,leftplusw) + depth;
    if (ISNOTEND(cptr))
    {
      subsort_bentleysedgewick(bsr,
                               subbucketleft + leftplusw,
                               bucketright-(pd-pb)-leftplusw,
                               depth+commonunitsequal);
    }

    gt_assert(pd >= pc);
    if ((wtmp = (unsigned long) (pd-pc)) > 0)
    {
      if (bsr->lcpsubtab != NULL && bsr->assideeffect)
      {
        /*
          right part has suffix with lcp up to length largermaxlcp w.r.t.
          to the pivot. This lcp belongs to a suffix on the right
          which is at a minimum distance to the pivot and thus to an
          element in the first part of the right side.
        */
        updatelcpvalue(bsr,subbucketleft + bucketright - wtmp + 1,
                       depth + greatermaxlcp);
      }
      subsort_bentleysedgewick(bsr,
                               subbucketleft + bucketright - wtmp + 1,
                               wtmp,
                               depth + greaterminlcp);
    }
  }
}

static unsigned long computelocallcpvalue(const Suffixwithcode *previoussuffix,
                                          const Suffixwithcode *currentsuffix,
                                          unsigned int minchanged)
{
  unsigned int lcpvalue;

  if (previoussuffix->code == currentsuffix->code)
  {
    lcpvalue = MIN(previoussuffix->prefixindex,
                   currentsuffix->prefixindex);
  } else
  {
    gt_assert(previoussuffix->code < currentsuffix->code);
    lcpvalue = MIN(minchanged,MIN(previoussuffix->prefixindex,
                                  currentsuffix->prefixindex));
  }
  return (unsigned long) lcpvalue;
}

static unsigned int bucketends(Outlcpinfo *outlcpinfo,
                               Suffixwithcode *previoussuffix,
                               GT_UNUSED unsigned long firstspecialsuffix,
                               unsigned int minchanged,
                               unsigned long specialsinbucket,
                               GtCodetype code,
                               const Bcktab *bcktab)
{
  unsigned long lcpvalue;
  unsigned int maxprefixindex, minprefixindex;
  Suffixwithcode firstspecialsuffixwithcode;

  /*
     there is at least one element in the bucket. if there is more than
     one element in the bucket, then we insert them using the
     information from the bcktab
  */
  if (specialsinbucket > 1UL)
  {
    maxprefixindex = gt_pfxidx2lcpvalues(&minprefixindex,
                                      outlcpinfo->lcpsubtab.smalllcpvalues,
                                      specialsinbucket,
                                      bcktab,
                                      code);
    if (outlcpinfo->lcpsubtab.maxbranchdepth < (unsigned long) maxprefixindex)
    {
      outlcpinfo->lcpsubtab.maxbranchdepth = (unsigned long) maxprefixindex;
    }
  } else
  {
    minprefixindex = maxprefixindex = gt_singletonmaxprefixindex(bcktab,code);
  }
  firstspecialsuffixwithcode.code = code;
  firstspecialsuffixwithcode.prefixindex = maxprefixindex;
#ifdef SKDEBUG
  firstspecialsuffixwithcode.startpos = firstspecialsuffix;
  /*
  consistencyofsuffix(__LINE__,
                      encseq,readmode,bcktab,numofchars,
                      &firstspecialsuffixwithcode);
  */
#endif
  lcpvalue = computelocallcpvalue(previoussuffix,
                                  &firstspecialsuffixwithcode,
                                  minchanged);
  if (outlcpinfo->lcpsubtab.maxbranchdepth < lcpvalue)
  {
    outlcpinfo->lcpsubtab.maxbranchdepth = lcpvalue;
  }
  outlcpinfo->lcpsubtab.smalllcpvalues[0] = (uint8_t) lcpvalue;
  outlcpinfo->lcpsubtab.countoutputlcpvalues += specialsinbucket;
  gt_xfwrite(outlcpinfo->lcpsubtab.smalllcpvalues,
             sizeof (*outlcpinfo->lcpsubtab.smalllcpvalues),
             (size_t) specialsinbucket,
             outlcpinfo->outfplcptab);
  return minprefixindex;
}

Outlcpinfo *gt_newOutlcpinfo(const char *indexname,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             unsigned long totallength,
                             bool assideeffect,
                             GtError *err)
{
  bool haserr = false;
  Outlcpinfo *outlcpinfo;

  outlcpinfo = gt_malloc(sizeof (*outlcpinfo));
  if (indexname == NULL)
  {
    outlcpinfo->outfplcptab = NULL;
    outlcpinfo->outfpllvtab = NULL;
  } else
  {
    outlcpinfo->outfplcptab = gt_fa_fopen_with_suffix(indexname,LCPTABSUFFIX,
                                                      "wb",err);
    if (outlcpinfo->outfplcptab == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      outlcpinfo->outfpllvtab
        = gt_fa_fopen_with_suffix(indexname,LARGELCPTABSUFFIX,"wb",err);
      if (outlcpinfo->outfpllvtab == NULL)
      {
        haserr = true;
      }
    }
  }
  outlcpinfo->assideeffect = assideeffect;
  outlcpinfo->lcpsubtab.countoutputlcpvalues = 0;
  outlcpinfo->totallength = totallength;
  outlcpinfo->lcpsubtab.totalnumoflargelcpvalues = 0;
  outlcpinfo->lcpsubtab.maxbranchdepth = 0;
  outlcpinfo->lcpsubtab.reservoir = NULL;
  outlcpinfo->lcpsubtab.sizereservoir = 0;
  GT_INITARRAY(&outlcpinfo->lcpsubtab.largelcpvalues,Largelcpvalue);
  outlcpinfo->lcpsubtab.smalllcpvalues = NULL;
  outlcpinfo->minchanged = 0;
  if (assideeffect)
  {
    outlcpinfo->tw = gt_newTurningwheel(prefixlength,numofchars);
  } else
  {
    outlcpinfo->tw = NULL;
  }
#ifdef SKDEBUG
  outlcpinfo->previoussuffix.startpos = 0;
#endif
  outlcpinfo->previoussuffix.code = 0;
  outlcpinfo->previoussuffix.prefixindex = 0;
  outlcpinfo->previoussuffix.defined = false;
  outlcpinfo->previousbucketwasempty = false;
  if (haserr)
  {
    gt_free(outlcpinfo);
    return NULL;
  }
  return outlcpinfo;
}

static void outlcpvalues(Lcpsubtab *lcpsubtab,
                         unsigned long bucketleft,
                         unsigned long bucketright,
                         unsigned long posoffset,
                         FILE *fplcptab,
                         FILE *fpllvtab)
{
  unsigned long idx;
  unsigned long lcpvalue;
  Largelcpvalue *largelcpvalueptr;

  lcpsubtab->largelcpvalues.nextfreeLargelcpvalue = 0;
  if (lcpsubtab->numoflargelcpvalues > 0 &&
      lcpsubtab->numoflargelcpvalues >=
      lcpsubtab->largelcpvalues.allocatedLargelcpvalue)
  {
    lcpsubtab->largelcpvalues.spaceLargelcpvalue
      = gt_realloc(lcpsubtab->largelcpvalues.spaceLargelcpvalue,
                   sizeof (Largelcpvalue) * lcpsubtab->numoflargelcpvalues);
    lcpsubtab->largelcpvalues.allocatedLargelcpvalue
      = lcpsubtab->numoflargelcpvalues;
  }
  for (idx=bucketleft; idx<=bucketright; idx++)
  {
    if (lcpsubtab->bucketoflcpvalues != NULL)
    {
      lcpvalue = lcpsubtab->bucketoflcpvalues[idx];
    } else
    {
      lcpvalue = compressedtable_get(lcpsubtab->completelcpvalues,idx);
    }
    if (lcpsubtab->maxbranchdepth < lcpvalue)
    {
      lcpsubtab->maxbranchdepth = lcpvalue;
    }
    if (lcpvalue < (unsigned long) LCPOVERFLOW)
    {
      lcpsubtab->smalllcpvalues[idx-bucketleft] = (uint8_t) lcpvalue;
    } else
    {
      gt_assert(lcpsubtab->largelcpvalues.nextfreeLargelcpvalue <
                lcpsubtab->largelcpvalues.allocatedLargelcpvalue);
      largelcpvalueptr = lcpsubtab->largelcpvalues.spaceLargelcpvalue +
                         lcpsubtab->largelcpvalues.nextfreeLargelcpvalue++;
      largelcpvalueptr->position = posoffset+idx;
      largelcpvalueptr->value = lcpvalue;
      lcpsubtab->smalllcpvalues[idx-bucketleft] = LCPOVERFLOW;
    }
  }
  lcpsubtab->countoutputlcpvalues += (bucketright - bucketleft + 1);
  gt_xfwrite(lcpsubtab->smalllcpvalues,
             sizeof (*lcpsubtab->smalllcpvalues),
             (size_t) (bucketright - bucketleft + 1),fplcptab);
  if (lcpsubtab->largelcpvalues.nextfreeLargelcpvalue > 0)
  {
    lcpsubtab->totalnumoflargelcpvalues
      += lcpsubtab->largelcpvalues.nextfreeLargelcpvalue;
    gt_xfwrite(lcpsubtab->largelcpvalues.spaceLargelcpvalue,
               sizeof (Largelcpvalue),
               (size_t)
               lcpsubtab->largelcpvalues.nextfreeLargelcpvalue,
               fpllvtab);
  }
}

#define NUMBEROFZEROS 1024

static unsigned long outmany0lcpvalues(unsigned long countoutputlcpvalues,
                                       unsigned long totallength,
                                       FILE *outfplcptab)
{
  unsigned long i, countout, many;
  uint8_t outvalues[NUMBEROFZEROS] = {0};

  many = totallength + 1 - countoutputlcpvalues;
  countout = many/NUMBEROFZEROS;
  for (i=0; i<countout; i++)
  {
    gt_xfwrite(outvalues,sizeof (uint8_t),(size_t) NUMBEROFZEROS,outfplcptab);
  }
  gt_xfwrite(outvalues,sizeof (uint8_t),(size_t) many % NUMBEROFZEROS,
             outfplcptab);
  return many;
}

static void multioutlcpvalues(Lcpsubtab *lcpsubtab,
                              unsigned long totallength,
                              const Compressedtable *lcptab,
                              unsigned long bucketsize,
                              FILE *fplcptab,
                              FILE *fpllvtab)
{
  unsigned long buffersize = 512UL, sizeforsmalllcpvalues;
  unsigned long remaining, left, width;
  bool mallocsmalllcpvalues;

  if (buffersize > totallength + 1)
  {
    buffersize = totallength+1;
  }
  lcpsubtab->numoflargelcpvalues = buffersize;
  lcpsubtab->bucketoflcpvalues = NULL;
  lcpsubtab->completelcpvalues = lcptab;
  sizeforsmalllcpvalues = (unsigned long)
                          sizeof (*lcpsubtab->smalllcpvalues) * buffersize;
  lcpsubtab->smalllcpvalues
    = compressedtable_unusedmem(lcptab,
                                (size_t) sizeforsmalllcpvalues);
  if (lcpsubtab->smalllcpvalues == NULL)
  {
    lcpsubtab->smalllcpvalues = gt_malloc((size_t) sizeforsmalllcpvalues);
    mallocsmalllcpvalues = true;
  } else
  {
    mallocsmalllcpvalues = false;
  }
  remaining = bucketsize;
  left = 0;
  gt_assert(fplcptab != NULL && fpllvtab != NULL);
  while (remaining > 0)
  {
    width = MIN(remaining, buffersize);
    outlcpvalues(lcpsubtab,
                 left,
                 left + width - 1,
                 0,
                 fplcptab,
                 fpllvtab);
    remaining -= width;
    left += width;
  }
  if (mallocsmalllcpvalues)
  {
    gt_free(lcpsubtab->smalllcpvalues);
  }
  lcpsubtab->countoutputlcpvalues = bucketsize;
}

void gt_freeOutlcptab(Outlcpinfo *outlcpinfo)
{
  if (outlcpinfo == NULL)
  {
    return;
  }
  if (outlcpinfo->assideeffect)
  {
    gt_free(outlcpinfo->lcpsubtab.reservoir);
    outlcpinfo->lcpsubtab.reservoir = NULL;
    outlcpinfo->lcpsubtab.sizereservoir = 0;
    if (outlcpinfo->tw != NULL)
    {
      gt_freeTurningwheel(&outlcpinfo->tw);
    }
  }
  if (outlcpinfo->lcpsubtab.countoutputlcpvalues < outlcpinfo->totallength+1)
  {
    outlcpinfo->lcpsubtab.countoutputlcpvalues
      += outmany0lcpvalues(outlcpinfo->lcpsubtab.countoutputlcpvalues,
                           outlcpinfo->totallength,
                           outlcpinfo->outfplcptab);
  }
  gt_assert(outlcpinfo->lcpsubtab.countoutputlcpvalues ==
            outlcpinfo->totallength + 1);
  GT_FREEARRAY(&outlcpinfo->lcpsubtab.largelcpvalues,Largelcpvalue);
  gt_fa_fclose(outlcpinfo->outfplcptab);
  gt_fa_fclose(outlcpinfo->outfpllvtab);
  gt_free(outlcpinfo);
}

unsigned long getnumoflargelcpvalues(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->lcpsubtab.totalnumoflargelcpvalues;
}

unsigned long getmaxbranchdepth(const Outlcpinfo *outlcpinfo)
{
  return outlcpinfo->lcpsubtab.maxbranchdepth;
}

static void initBentsedgresources(Bentsedgresources *bsr,
                                  GtSuffixsortspace *suffixsortspace,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  Bcktab *bcktab,
                                  GtCodetype mincode,
                                  GtCodetype maxcode,
                                  unsigned long partwidth,
                                  unsigned int numofchars,
                                  unsigned int prefixlength,
                                  Outlcpinfo *outlcpinfo,
                                  const Sfxstrategy *sfxstrategy)
{
  unsigned long idx;

  bsr->readmode = readmode;
  bsr->totallength = gt_encseq_total_length(encseq);
  bsr->sfxstrategy = sfxstrategy;
  bsr->sssp = suffixsortspace;
  gt_suffixsortspace_bucketleftidx_set(bsr->sssp,0);
  bsr->encseq = encseq;
  bsr->fwd = GT_ISDIRREVERSE(bsr->readmode) ? false : true;
  bsr->complement = GT_ISDIRCOMPLEMENT(bsr->readmode) ? true : false;
  for (idx = 0; idx < (unsigned long) GT_UNITSIN2BITENC; idx++)
  {
    bsr->leftlcpdist[idx] = bsr->rightlcpdist[idx] = 0;
  }
  if (outlcpinfo != NULL)
  {
    bsr->lcpsubtab = &outlcpinfo->lcpsubtab;
    bsr->assideeffect = outlcpinfo->assideeffect;
  } else
  {
    bsr->lcpsubtab = NULL;
  }
  bsr->esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  bsr->esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  if (bcktab != NULL)
  {
    gt_determinemaxbucketsize(bcktab,
                              mincode,
                              maxcode,
                              partwidth,
                              numofchars,
                              false,
                              0); /* not necesarry as hashexceptions = false */
    /* gt_bcktab_showlog2info(bcktab,logger); */
    if (outlcpinfo != NULL && outlcpinfo->assideeffect)
    {
      size_t sizeforlcpvalues; /* in bytes */

      gt_assert(bsr->lcpsubtab != NULL);
      gt_assert(sizeof (*bsr->lcpsubtab->smalllcpvalues) == sizeof (uint8_t));
      gt_assert(sizeof (*bsr->lcpsubtab->bucketoflcpvalues) ==
                sizeof (unsigned long));
      sizeforlcpvalues = gt_bcktab_sizeforlcpvalues(bcktab);
      if (bsr->lcpsubtab->sizereservoir < sizeforlcpvalues)
      {
        bsr->lcpsubtab->sizereservoir = sizeforlcpvalues;
        bsr->lcpsubtab->reservoir = gt_realloc(bsr->lcpsubtab->reservoir,
                                               bsr->lcpsubtab->sizereservoir);
        /* point to the same area, since this is not used simultaneously */
        /* be careful for the parallel version */
        bsr->lcpsubtab->smalllcpvalues = (uint8_t *) bsr->lcpsubtab->reservoir;
        bsr->lcpsubtab->bucketoflcpvalues
          = (unsigned long *) bsr->lcpsubtab->reservoir;
      }
    }
  }
  GT_INITARRAY(&bsr->mkvauxstack,MKVstack);
  if (sfxstrategy->cmpcharbychar)
  {
    bsr->countingsortinfo = NULL;
    bsr->medianinfospace = NULL;
  } else
  {
    bsr->countingsortinfo = gt_malloc(sizeof (*bsr->countingsortinfo) *
                                      sfxstrategy->maxcountingsort);
    if (sfxstrategy->maxwidthrealmedian >= MINMEDIANOF9WIDTH)
    {
      bsr->medianinfospace = gt_malloc(sizeof (*bsr->medianinfospace) *
                                       sfxstrategy->maxwidthrealmedian);
    } else
    {
      bsr->medianinfospace = NULL;
    }
  }
  if (bcktab != NULL && sfxstrategy->ssortmaxdepth.defined)
  {
    bsr->rmnsufinfo = gt_rmnsufinfo_new(suffixsortspace,
                                        -1,
                                        NULL,
                                        bsr->encseq,
                                        bcktab,
                                        maxcode,
                                        numofchars,
                                        prefixlength,
                                        readmode,
                                        partwidth,
                                        false,
                                        true);
    gt_assert(bsr->rmnsufinfo != NULL);
  } else
  {
    bsr->rmnsufinfo = NULL;
  }
  if (sfxstrategy->ssortmaxdepth.defined)
  {
    bsr->blindtrie = NULL;
  } else
  {
    bsr->blindtrie = gt_blindtrie_new(bsr->sssp,
                                      sfxstrategy->maxbltriesort,
                                      0, /* the nodenumberincrement */
                                      encseq,
                                      sfxstrategy->cmpcharbychar,
                                      bsr->esr1,
                                      bsr->esr2,
                                      readmode);
  }
  bsr->voiddcov = NULL;
  bsr->dc_processunsortedrange = NULL;
  if (bsr->sfxstrategy->ssortmaxdepth.defined ||
      bsr->sfxstrategy->differencecover > 0)
  {
    bsr->equalwithprevious = gt_malloc(sizeof (*bsr->equalwithprevious) *
                                       bsr->sfxstrategy->maxinsertionsort);
    for (idx=0; idx < bsr->sfxstrategy->maxinsertionsort; idx++)
    {
      bsr->equalwithprevious[idx] = false;
    }
  } else
  {
    bsr->equalwithprevious = NULL;
  }
  bsr->countinsertionsort = 0;
  bsr->countqsort = 0;
  bsr->countcountingsort = 0;
  bsr->countbltriesort = 0;
}

static void wrapBentsedgresources(Bentsedgresources *bsr,
                                  unsigned long partwidth,
                                  Lcpsubtab *lcpsubtab,
                                  FILE *outfplcptab,
                                  FILE *outfpllvtab,
                                  GtLogger *logger)
{
  gt_free(bsr->countingsortinfo);
  bsr->countingsortinfo = NULL;
  gt_free(bsr->medianinfospace);
  bsr->medianinfospace = NULL;
  gt_blindtrie_delete(bsr->blindtrie);
  if (bsr->rmnsufinfo != NULL)
  {
    Compressedtable *lcptab;

    lcptab = gt_rmnsufinfo_delete(&bsr->rmnsufinfo,
                                  bsr->lcpsubtab == NULL ? false : true);
    if (lcptab != NULL)
    {
      multioutlcpvalues(lcpsubtab,bsr->totallength,lcptab,partwidth,
                        outfplcptab,outfpllvtab);
      compressedtable_free(lcptab,true);
    }
  }
  if (bsr->esr1 != NULL)
  {
    gt_encseq_reader_delete(bsr->esr1);
  }
  if (bsr->esr2 != NULL)
  {
    gt_encseq_reader_delete(bsr->esr2);
  }
  gt_free(bsr->equalwithprevious);
  GT_FREEARRAY(&bsr->mkvauxstack,MKVstack);
  gt_logger_log(logger,"countinsertionsort=%lu",bsr->countinsertionsort);
  gt_logger_log(logger,"countbltriesort=%lu",bsr->countbltriesort);
  gt_logger_log(logger,"countcountingsort=%lu",bsr->countcountingsort);
  gt_logger_log(logger,"countqsort=%lu",bsr->countqsort);
}

void gt_qsufsort(GtSuffixsortspace *suffixsortspace,
                 unsigned long partwidth,
                 int mmapfiledesc,
                 GtStr *mmapfilename,
                 const GtEncseq *encseq,
                 GtReadmode readmode,
                 GT_UNUSED GtCodetype mincode,
                 GtCodetype maxcode,
                 Bcktab *bcktab,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 bool hashexceptions,
                 bool absoluteinversesuftab,
                 Outlcpinfo *outlcpinfo)
{
  Rmnsufinfo *rmnsufinfo;
  Compressedtable *lcptab;

  gt_assert(mincode == 0);
  rmnsufinfo = gt_rmnsufinfo_new(suffixsortspace,
                                 mmapfiledesc,
                                 mmapfilename,
                                 encseq,
                                 bcktab,
                                 maxcode,
                                 numofchars,
                                 prefixlength,
                                 readmode,
                                 partwidth,
                                 hashexceptions,
                                 absoluteinversesuftab);
  gt_rmnsufinfo_bcktab2firstlevelintervals(rmnsufinfo);
  lcptab = gt_rmnsufinfo_delete(&rmnsufinfo,
                                outlcpinfo == NULL ? false : true);
  if (lcptab != NULL)
  {
    gt_assert(outlcpinfo != NULL);
    gt_assert(outlcpinfo->outfplcptab != NULL);
    gt_assert(outlcpinfo->outfpllvtab != NULL);

    multioutlcpvalues(&outlcpinfo->lcpsubtab,gt_encseq_total_length(encseq),
                      lcptab,partwidth,outlcpinfo->outfplcptab,
                      outlcpinfo->outfpllvtab);
    compressedtable_free(lcptab,true);
  }
}

/*
  The following function is called  in sfxsuffixer.c sorts all buckets by
  different suffix comparison methods without the help of other sorting
  information. GtSuffixsortspace contains the sortspace which is accessed
  by some negative offset.
*/

void gt_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       GtBucketspec2 *bucketspec2,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       unsigned long partwidth,
                       Bcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       Outlcpinfo *outlcpinfo,
                       const Sfxstrategy *sfxstrategy,
                       unsigned long long *bucketiterstep,
                       GtLogger *logger)
{
  GtCodetype code;
  unsigned int rightchar = (unsigned int) (mincode % numofchars),
               minprefixindex;
  Bucketspecification bucketspec;
  unsigned long lcpvalue;
  Suffixwithcode firstsuffixofbucket;
  Bentsedgresources bsr;

  initBentsedgresources(&bsr,
                        suffixsortspace,
                        encseq,
                        readmode,
                        bcktab,
                        mincode,
                        maxcode,
                        partwidth,
                        numofchars,
                        prefixlength,
                        outlcpinfo,
                        sfxstrategy);
  for (code = mincode; code <= maxcode; code++)
  {
    if (bucketspec2 != NULL)
    {
      if (gt_copysort_checkhardwork(bucketspec2,code))
      {
        rightchar = (unsigned int) (code % numofchars);
      } else
      {
        continue;
      }
    }
    (*bucketiterstep)++;
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                         bcktab,
                                         code,
                                         maxcode,
                                         partwidth,
                                         rightchar,
                                         numofchars);
    if (outlcpinfo != NULL && outlcpinfo->assideeffect)
    {
      bsr.lcpsubtab->numoflargelcpvalues = 0;
      if (code > 0)
      {
        (void) gt_nextTurningwheel(outlcpinfo->tw);
        if (outlcpinfo->previousbucketwasempty)
        {
          outlcpinfo->minchanged = MIN(outlcpinfo->minchanged,
                                     gt_minchangedTurningwheel(outlcpinfo->tw));
        } else
        {
          outlcpinfo->minchanged = gt_minchangedTurningwheel(outlcpinfo->tw);
        }
      }
    }
    if (bucketspec.nonspecialsinbucket > 0)
    {
      if (bucketspec.nonspecialsinbucket > 1UL)
      {
        if (outlcpinfo != NULL && outlcpinfo->assideeffect)
        {
          gt_assert(bsr.lcpsubtab != NULL);
        }
        gt_suffixsortspace_bucketleftidx_set(bsr.sssp,bucketspec.left);
        bentleysedgewick(&bsr,
                         bucketspec.nonspecialsinbucket,
                         (unsigned long) prefixlength);
        gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
      }
      if (outlcpinfo != NULL && outlcpinfo->assideeffect)
      {
        if (outlcpinfo->previoussuffix.defined)
        {
          /* compute lcpvalue of first element of bucket with
             last element of previous bucket */
          firstsuffixofbucket.code = code;
          firstsuffixofbucket.prefixindex = prefixlength;
#ifdef SKDEBUG
          firstsuffixofbucket.startpos
            = gt_suffixsortspace_get(bsr.sssp,0,bucketspec.left);
          /*
          consistencyofsuffix(__LINE__,
                              encseq,readmode,bcktab,numofchars,
                              &firstsuffixofbucket);
          */
#endif
          lcpvalue = computelocallcpvalue(&outlcpinfo->previoussuffix,
                                          &firstsuffixofbucket,
                                          outlcpinfo->minchanged);
        } else
        {
          /* first part first code */
          lcpvalue = 0;
        }
        gt_assert(bsr.lcpsubtab != NULL);
#ifdef SKDEBUG
        baseptr = bucketspec.left;
#endif
        updatelcpvalue(&bsr,0,lcpvalue);
        /* all other lcp-values are computed and they can be output */
        outlcpvalues(&outlcpinfo->lcpsubtab,
                     0,
                     bucketspec.nonspecialsinbucket-1,
                     bucketspec.left,
                     outlcpinfo->outfplcptab,
                     outlcpinfo->outfpllvtab);
        /* previoussuffix becomes last nonspecial element in current bucket */
        outlcpinfo->previoussuffix.code = code;
        outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
        outlcpinfo->previoussuffix.startpos
          = gt_suffixsortspace_get(bsr.sssp,0,
                                   bucketspec.left
                                     + bucketspec.nonspecialsinbucket - 1);
        /*
        consistencyofsuffix(__LINE__,
                            encseq,readmode,bcktab,numofchars,
                            &outlcpinfo->previoussuffix);
        */
#endif
      }
    }
    if (outlcpinfo != NULL && outlcpinfo->assideeffect)
    {
      if (bucketspec.specialsinbucket > 0)
      {
        unsigned long suffixvalue
          = gt_suffixsortspace_get(suffixsortspace,
                                   0,
                                   bucketspec.left
                                     + bucketspec.nonspecialsinbucket);
        minprefixindex = bucketends(outlcpinfo,
                                    &outlcpinfo->previoussuffix,
                                    /* first special element in bucket */
                                    suffixvalue,
                                    outlcpinfo->minchanged,
                                    bucketspec.specialsinbucket,
                                    code,
                                    bcktab);
        /* there is at least one special element: this is the last element
           in the bucket, and thus the previoussuffix for the next round */
        outlcpinfo->previoussuffix.defined = true;
        outlcpinfo->previoussuffix.code = code;
        outlcpinfo->previoussuffix.prefixindex = minprefixindex;
#ifdef SKDEBUG
        outlcpinfo->previoussuffix.startpos
          = gt_suffixsortspace_get(suffixsortspace,
                                   0,
                                   bucketspec.left
                                     + bucketspec.nonspecialsinbucket +
                                       bucketspec.specialsinbucket - 1);
        /*
         consistencyofsuffix(__LINE__,
                             encseq,readmode,bcktab,numofchars,
                             &outlcpinfo->previoussuffix);
        */
#endif
      } else
      {
        if (bucketspec.nonspecialsinbucket > 0)
        {
          /* if there is at least one element in the bucket, then the last
             one becomes the next previous suffix */
          outlcpinfo->previoussuffix.defined = true;
          outlcpinfo->previoussuffix.code = code;
          outlcpinfo->previoussuffix.prefixindex = prefixlength;
#ifdef SKDEBUG
          outlcpinfo->previoussuffix.startpos
            = gt_suffixsortspace_get(suffixsortspace,
                                     0,
                                     bucketspec.left
                                       + bucketspec.nonspecialsinbucket-1);
          /*
          consistencyofsuffix(__LINE__,
                              encseq,readmode,bcktab,numofchars,
                              &outlcpinfo->previoussuffix);
          */
#endif
        }
      }
      if (bucketspec.nonspecialsinbucket + bucketspec.specialsinbucket == 0)
      {
        outlcpinfo->previousbucketwasempty = true;
      } else
      {
        outlcpinfo->previousbucketwasempty = false;
      }
    }
  }
  wrapBentsedgresources(&bsr,
                        partwidth,
                        outlcpinfo == NULL ? NULL : &outlcpinfo->lcpsubtab,
                        outlcpinfo == NULL ? NULL : outlcpinfo->outfplcptab,
                        outlcpinfo == NULL ? NULL : outlcpinfo->outfpllvtab,
                        logger);
}

/*
  The following function is from match/sfx-diffcov.c, but we do not want to
  include match/sfx-diffcov.h here, as this would result into
  cyclic dependencies. So we directly put the forward declaration here */

void dc_setsuffixsortspace(void *voiddcov,GtSuffixsortspace *sssp);

/*
   The following function is used for sorting the sample making up the
   difference cover and for sorting with the difference cover.
*/

void gt_sortbucketofsuffixes(bool setdcovsuffixsortspace,
                             GtSuffixsortspace *suffixsortspace,
                             unsigned long numberofsuffixes,
                             GtBucketspec2 *bucketspec2,
                             const GtEncseq *encseq,
                             GtReadmode readmode,
                             GtCodetype mincode,
                             GtCodetype maxcode,
                             const Bcktab *bcktab,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             const Sfxstrategy *sfxstrategy,
                             void *voiddcov,
                             Dc_processunsortedrange dc_processunsortedrange,
                             GtLogger *logger)
{
  Bentsedgresources bsr;
  Bucketspecification bucketspec;
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtCodetype code;

  if (setdcovsuffixsortspace)
  {
    dc_setsuffixsortspace(voiddcov,suffixsortspace);
  }
  initBentsedgresources(&bsr,
                        suffixsortspace,
                        encseq,
                        readmode,
                        NULL, /* bcktab unused */
                        0,    /* mincode unused */
                        0,    /* maxcode unused */
                        0,    /* partwidth unused */
                        numofchars,
                        prefixlength,
                        NULL,  /* outlcpinfo unused */
                        sfxstrategy);
  bsr.voiddcov = voiddcov;
  bsr.dc_processunsortedrange = dc_processunsortedrange;
  for (code = mincode; code <= maxcode; code++)
  {
    if (bucketspec2 != NULL)
    {
      if (gt_copysort_checkhardwork(bucketspec2,code))
      {
        rightchar = (unsigned int) (code % numofchars);
      } else
      {
        continue;
      }
    }
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                         bcktab,
                                         code,
                                         maxcode,
                                         numberofsuffixes,
                                         rightchar,
                                         numofchars);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      /*fprintf(stderr,"set bucketleftidx = %lu\n",bsr.sssp->bucketleftidx);*/
      gt_suffixsortspace_bucketleftidx_set(bsr.sssp,bucketspec.left);
      bentleysedgewick(&bsr,
                       bucketspec.nonspecialsinbucket,
                       (unsigned long) prefixlength);
      gt_suffixsortspace_bucketleftidx_set(bsr.sssp,0);
    }
  }
  wrapBentsedgresources(&bsr,
                        0, /* partwidth value unused because lcptab == NULL */
                        NULL,
                        NULL,
                        NULL,
                        logger);
}
