/*
  Copyright (c) 2010 Jiabin Huang <jiabin.huang@studium.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef MAXMAT4_DFS_BITTAB_H
#define MAXMAT4_DFS_BITTAB_H

#include "core/stack-inlined.h"
#ifndef S_SPLINT_S
#include "match/eis-bwtseq.h"
#endif
#include "core/bittab_api.h"

typedef struct {
  unsigned long depth,
                  lower,
                  upper;
  GtBittab *prefixofsuffixbits;
  GtCodetype code;
} Maxmat4NodeBittab;

GT_STACK_DECLARESTRUCT(Maxmat4NodeBittab, 256UL);

#ifndef S_SPLINT_S
int gt_pck_bitparallelism_bittab(const GtUchar *query,
                                  unsigned long querylen,
                                  const BWTSeq *bwtSeq,
                                  const GtEncseq *encseq,
                                  const Mbtab **mbtab,
                                  unsigned int maxdepth,
                                  unsigned long totallength,
                                  unsigned long leastlength,
                                  const GtMatchmode matchmode,
                                  Processmatchfunction processmatch,
                                  Showspecinfo *showspecinfo,
                                  unsigned long bitlength,
                                  bool showbitparallelismfactor,
                                  bool showtime,
                                  GtProgressTimer *timer,
                                  GT_UNUSED GtLogger *logger,
                                  GtError *err);
#endif

#endif
