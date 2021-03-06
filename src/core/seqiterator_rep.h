/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef SEQITERATOR_REP_H
#define SEQITERATOR_REP_H

#include "core/seqiterator.h"

typedef void       (*GtSeqIteratorSetSymbolmapFunc)(GtSeqIterator*,
                                                    const GtUchar *symbolmap);
typedef void       (*GtSeqIteratorSetSequenceOutFunc)(GtSeqIterator*, bool);
typedef int        (*GtSeqIteratorNextFunc)(GtSeqIterator *seq_iterator,
                                            const GtUchar **sequence,
                                            unsigned long *len,
                                            char **description, GtError*);
typedef const unsigned long long*
                   (*GtSeqIteratorGetCurrCounterFunc)(GtSeqIterator*,
                                                      unsigned long long);
typedef void       (*GtSeqIteratorSetQualBufferFunc)(GtSeqIterator*,
                                                     const GtUchar **qualities);
typedef void       (*GtSeqIteratorFreeFunc)(GtSeqIterator*);

struct GtSeqIterator {
  const GtSeqIteratorClass *c_class;
};

GtSeqIterator*            gt_seqiterator_create(const GtSeqIteratorClass*);
void*                     gt_seqiterator_cast(const GtSeqIteratorClass*,
                                              GtSeqIterator*);
const GtSeqIteratorClass* gt_seqiterator_class_new(size_t size,
                         GtSeqIteratorSetSymbolmapFunc set_symbolmap,
                         GtSeqIteratorSetSequenceOutFunc set_seqout,
                         GtSeqIteratorNextFunc next_func,
                         GtSeqIteratorGetCurrCounterFunc get_curr_counter_func,
                         GtSeqIteratorSetQualBufferFunc set_qual_buffer_func,
                         GtSeqIteratorFreeFunc free_func);

#endif
