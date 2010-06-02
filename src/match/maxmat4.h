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

#ifndef MAXMAT4_H
#define MAXMAT4_H

#include "core/defined-types.h"
#include "core/encseq.h"
#include "match/matchmode_api.h"

int gt_findmum(const GtEncseq *encseq,
               const void *genericindex,
               unsigned long totallength,
               const GtAlphabet *alphabet,
               const GtStrArray *queryfilenames,
               GtMatchmode matchmode,
               Definedunsignedlong leastlength,
               bool nucleotidesonly,
               bool bothdirections,
               bool reversecomplement,
               bool showstring,
               bool showreversepositions,
               bool showsequencelengths,
               bool verbose,
               GtError *err);

#endif
