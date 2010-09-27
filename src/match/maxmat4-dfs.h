/*
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

#ifndef MAXMAT4_DFS_H
#define MAXMAT4_DFS_H

#include "core/stack-inlined.h"
#include "match/eis-voiditf.h"

typedef struct {
  unsigned long depth,
                prefixofsuffixbits,
                lower,
                upper;
} Maxmat4Node;

GT_STACK_DECLARESTRUCT(Maxmat4Node, 256UL);

int gt_pck_bitparallelism(const GtUchar *query,
                          unsigned long querylen,
                          const FMindex *index,
                          const GtEncseq *encseq,
                          unsigned long totallength,
                          unsigned long leastlength,
                          //Findmatchfunction findmatchfunction,
                          Processmatchfunction processmatch,
                          Showspecinfo *showspecinfo,
                          bool showtime,
                          GtProgressTimer *timer,
                          GtLogger *logger,
                          GT_UNUSED GtError *err);

#endif
