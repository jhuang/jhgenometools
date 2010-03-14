/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core/dataalign.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/minmax.h"
#include "core/seqiterator.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/symboldef.h"
#include "core/unused_api.h"
#include "match/seqpos-def.h"
#include "match/verbose-def.h"
#include "match/intcode-def.h"
#include "match/encseq-def.h"

#include "match/eis-encidxseq.h"
#include "match/eis-sa-common.h"
#include "match/eis-sequencemultiread.h"
#include "match/eis-suffixerator-interface.h"

struct sfxInterface
{
  struct SASeqSrc baseClass;
  Readmode readmode;
  unsigned int prefixlength, numofparts;
  const Sfxstrategy *sfxstrategy;
  const GtAlphabet *alpha;
  const Encodedsequence *encseq;
  struct seqStats *stats;
  Sfxiterator *sfi;
  DefinedSeqpos rot0Pos;
  bool specialsuffixes;
  /* data relevant to holding portions of the suffix array */
  Seqpos lastGeneratedLen, lastGeneratedStart;
  const Seqpos *lastGeneratedSufTabSegment;
};

static SeqDataTranslator
SfxIRequest2XltorFunc(sfxInterface *sfxi,
                      enum sfxDataRequest rtype)
{
  SeqDataTranslator tr = { { NULL }, NULL };
  switch (rtype)
  {
    union saXltorState readState;
    struct saTaggedXltorState *stateStore;
  case SFX_REQUEST_BWTTAB:
    readState.encSeqTr.readmode = sfxi->readmode;
    readState.encSeqTr.encseq = sfxi->encseq;
    stateStore = addSuffixarrayXltor(&sfxi->baseClass.xltorStates,
                                     rtype, readState);
    tr.state.ref = &stateStore->state.encSeqTr;
    tr.translateData = (seqDataTranslateFunc)translateSuftab2BWT;
    break;
  case SFX_REQUEST_SUFTAB:
    tr.state.elemSize = sizeof (Seqpos);
    break;
  case SFX_REQUEST_LCPTAB:
    readState.lcpState.readmode = sfxi->readmode;
    readState.lcpState.encseq = sfxi->encseq;
    readState.lcpState.lastSufIdx = -1;
    stateStore = addSuffixarrayXltor(&sfxi->baseClass.xltorStates,
                                     rtype, readState);
    tr.state.ref = &stateStore->state.lcpState;
    tr.translateData = (seqDataTranslateFunc)translateSuftab2LCP;
    break;
  default:
    fprintf(stderr, "error: unimplemented request!\n");
    abort();
    break;
  }
  return tr;
}

static inline sfxInterface *
SASS2SfxI(SASeqSrc *baseClass)
{
  return (sfxInterface *)((char *)baseClass
                          - offsetof(sfxInterface, baseClass));
}

extern struct SASeqSrc *
SfxI2SASS(sfxInterface *sfxi)
{
  return &sfxi->baseClass;
}

static inline const sfxInterface *
constSASS2SfxI(const SASeqSrc *baseClass)
{
  return (const sfxInterface *)((const char *)baseClass
                                - offsetof(sfxInterface, baseClass));
}

static SeqDataTranslator
SfxIBaseRequest2XltorFunc(SASeqSrc *baseClass,
                          enum sfxDataRequest rtype)
{
  return SfxIRequest2XltorFunc(SASS2SfxI(baseClass), rtype);
}

static DefinedSeqpos
SfxIBaseGetRot0Pos(const SASeqSrc *baseClass)
{
  return SfxIGetRot0Pos(constSASS2SfxI(baseClass));
}

static const struct seqStats *
SfxIBaseGetSeqStats(const SASeqSrc *baseClass)
{
  return SfxIGetSeqStats(constSASS2SfxI(baseClass));
}

static MRAEnc *
SfxIBaseNewMRAEnc(const SASeqSrc *baseClass)
{
  return SfxINewMRAEnc(constSASS2SfxI(baseClass));
}

static void
deleteSfxInterfaceBase(SASeqSrc *baseClass)
{
  deleteSfxInterface(SASS2SfxI(baseClass));
}

static size_t
SfxIGenerate(void *iface, void *backlogState,
             move2BacklogFunc move2Backlog, void *output, Seqpos generateStart,
             size_t len, SeqDataTranslator xltor);

extern sfxInterface *
newSfxInterface(Readmode readmode,
                unsigned int prefixlength,
                unsigned int numofparts,
                const Sfxstrategy *sfxstrategy,
                const Encodedsequence *encseq,
                Sfxprogress *sfxprogress,
                Seqpos length,
                Verboseinfo *verbosity,
                GtError *err)
{
  return newSfxInterfaceWithReaders(readmode,
                                    prefixlength,
                                    numofparts,
                                    sfxstrategy,
                                    0,
                                    NULL,
                                    NULL,
                                    encseq,
                                    sfxprogress,
                                    length,
                                    verbosity,
                                    err);
}

static struct seqStats *
newSeqStatsFromCharDist(const Encodedsequence *encseq,
                        const GtAlphabet *alpha, Seqpos len)
{
  struct seqStats *stats = NULL;
  unsigned i, numofchars;
  Seqpos regularSymsSum = 0;
  stats = gt_malloc(offsetAlign(sizeof (*stats), sizeof (Seqpos))
                    + (UINT8_MAX + 1) * sizeof (Seqpos));
  unsigned int numOfSeqs;

  numOfSeqs = getencseqnumofdbsequences(encseq);
  stats->sourceAlphaType = sourceUInt8;
  stats->symbolDistributionTable =
    (Seqpos *)((char *)stats + offsetAlign(sizeof (*stats), sizeof (Seqpos)));
  memset(stats->symbolDistributionTable, 0, sizeof (Seqpos) * (UINT8_MAX + 1));
  numofchars = gt_alphabet_num_of_chars(alpha);
  for (i = 0; i < numofchars; ++i)
  {
    stats->symbolDistributionTable[i]
      = (Seqpos) getencseqcharactercount(encseq,(GtUchar) i);
    regularSymsSum += stats->symbolDistributionTable[i];
  }
  stats->symbolDistributionTable[WILDCARD] = len - regularSymsSum - numOfSeqs;
  stats->symbolDistributionTable[SEPARATOR] += numOfSeqs;
  stats->symbolDistributionTable[UNDEFBWTCHAR] += 1;
  return stats;
}

static void
deleteSeqStats(struct seqStats *stats)
{
  gt_free(stats);
}

#define newSfxInterfaceWithReadersErrRet()        \
  do {                                            \
    if (sfxi->stats)                             \
      deleteSeqStats(sfxi->stats);               \
    if (sfxi) gt_free(sfxi);                    \
    sfxi = NULL;                                 \
  } while (0)

extern sfxInterface *
newSfxInterfaceWithReaders(Readmode readmode,
                           unsigned int prefixlength,
                           unsigned int numofparts,
                           const Sfxstrategy *sfxstrategy,
                           size_t numReaders,
                           enum sfxDataRequest readerRequests[],
                           SeqDataReader readers[],
                           const Encodedsequence *encseq,
                           Sfxprogress *sfxprogress,
                           Seqpos length,
                           Verboseinfo *verbosity, GtError *err)
{
  sfxInterface *sfxi = NULL;
  gt_error_check(err);

  sfxi = gt_calloc(1, sizeof (*sfxi));
  {
    RandomSeqAccessor origSeqAccess = { SfxIGetOrigSeq, sfxi };
    initSASeqSrc(&sfxi->baseClass, length, SfxIBaseRequest2XltorFunc, NULL,
                 SfxIBaseGetRot0Pos, SfxIBaseGetSeqStats,
                 origSeqAccess, deleteSfxInterfaceBase, SfxIBaseNewMRAEnc,
                 SfxIGenerate, sfxi);
  }
  sfxi->readmode = readmode;
  sfxi->encseq = encseq;
  sfxi->alpha = getencseqAlphabet(encseq);
  sfxi->stats = newSeqStatsFromCharDist(encseq,sfxi->alpha, length);
  if (!(sfxi->sfi = newSfxiterator(encseq,
                                   readmode,
                                   prefixlength,
                                   numofparts,
                                   NULL,
                                   sfxstrategy,
                                   sfxprogress,
                                   verbosity, err)))
    newSfxInterfaceWithReadersErrRet();
  sfxi->rot0Pos.defined = false;
  sfxi->specialsuffixes = false;

  sfxi->lastGeneratedStart = sfxi->lastGeneratedLen = 0;
  sfxi->lastGeneratedSufTabSegment = NULL;

  {
    size_t i;
    for (i = 0; i < numReaders; ++i)
    {
      readers[i] = SfxIRegisterReader(sfxi, readerRequests[i]);
      if (!readers[i].readData)
        newSfxInterfaceWithReadersErrRet();
    }
  }
  return sfxi;
}

extern const Sfxiterator *SfxInterface2Sfxiterator(const sfxInterface *sfxi)
{
  return sfxi->sfi;
}

extern void
deleteSfxInterface(sfxInterface *sfxi)
{
  destructSASeqSrc(&sfxi->baseClass);
  freeSfxiterator(&sfxi->sfi);
  deleteSeqStats(sfxi->stats);
  gt_free(sfxi);
}

extern const GtAlphabet *
SfxIGetAlphabet(const sfxInterface *si)
{
  return si->alpha;
}

extern MRAEnc *
SfxINewMRAEnc(const sfxInterface *si)
{
  MRAEnc *alphabet;
  gt_assert(si);
  alphabet = MRAEncGTAlphaNew(SfxIGetAlphabet(si));
  MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}

Seqpos
SfxIGetLength(const sfxInterface *si)
{
  gt_assert(si);
  return si->baseClass.seqLen;
}

extern const struct seqStats *
SfxIGetSeqStats(const sfxInterface *si)
{
  return si->stats;
}

extern DefinedSeqpos
SfxIGetRot0Pos(const struct sfxInterface *si)
{
  return si->rot0Pos;
}

extern const Encodedsequence *
SfxIGetEncSeq(const sfxInterface *si)
{
  return si->encseq;
}

extern Readmode
SfxIGetReadmode(const sfxInterface *si)
{
  return si->readmode;
}

extern SeqDataReader
SfxIRegisterReader(sfxInterface *sfxi, enum sfxDataRequest rtype)
{
  return seqReaderSetRegisterConsumer(
    &sfxi->baseClass.readerSet, rtype, SfxIRequest2XltorFunc(sfxi, rtype));
}

extern size_t
SfxIGetOrigSeq(const void *state, Symbol *dest, Seqpos pos, size_t len)
{
  const struct sfxInterface *sfxi;
  gt_assert(state);
  sfxi = state;
  return EncSeqGetSubSeq(sfxi->encseq, sfxi->readmode, pos, len, dest);
}

/** writes substring of suffix table to output, stores older data into
 * cache if necessary */
static size_t
SfxIGenerate(void *iface, void *backlogState,
             move2BacklogFunc move2Backlog, void *output, Seqpos generateStart,
             size_t len, SeqDataTranslator xltor)
{
  sfxInterface *sfxi = iface;
  size_t elemsLeft = len;
  gt_assert(sfxi && backlogState && move2Backlog && output);
  gt_assert(generateStart + len <= SfxIGetLength(sfxi));
  do
  {
    if (generateStart < sfxi->lastGeneratedStart + sfxi->lastGeneratedLen)
    {
      size_t copyLen = MIN(elemsLeft, sfxi->lastGeneratedStart
                           + sfxi->lastGeneratedLen - generateStart),
        charsWritten =
        SDRTranslate(xltor, output, sfxi->lastGeneratedSufTabSegment
                     + generateStart - sfxi->lastGeneratedStart, copyLen);
      generateStart += copyLen;
      elemsLeft -= copyLen;
      output = (char *)output + charsWritten;
    }
    /* 1. read next region of sequence by calling nextSfxIterator */
    if (elemsLeft)
    {
      move2Backlog(backlogState, sfxi->lastGeneratedSufTabSegment,
                   sfxi->lastGeneratedStart, sfxi->lastGeneratedLen);
      sfxi->lastGeneratedStart += sfxi->lastGeneratedLen;
      if ((sfxi->lastGeneratedSufTabSegment =
           nextSfxiterator(&sfxi->lastGeneratedLen, &sfxi->specialsuffixes,
                           sfxi->sfi)))
      {
        size_t pos;
        /* size_t because the current approach cannot generate more
         * than memory will hold anyway */
        size_t lastGeneratedLen = sfxi->lastGeneratedLen;
        const Seqpos *suftab = sfxi->lastGeneratedSufTabSegment;
        if (!sfxi->rot0Pos.defined)
          for (pos=0; pos < lastGeneratedLen; pos++)
          {
            if (suftab[pos] == 0)
            {
              sfxi->rot0Pos.defined = true;
              sfxi->rot0Pos.valueseqpos = sfxi->lastGeneratedStart + pos;
              break;
            }
          }
        /* uncomment this to reenable synchronous writing of tables */
/*if (sfxi->lastGeneratedSufTabSegment == NULL */
/*    || suftab2file(&sfxi->outfileinfo, sfxi->lastGeneratedSufTabSegment, */
/*                   sfxi->so.readmode, sfxi->lastGeneratedLen, err) != 0) */
/*       break; */
      }
      else
        break;
    }
    /* 5. if positions in region don't suffice go back to step 3. */
  } while (elemsLeft);
  return len - elemsLeft;
}
