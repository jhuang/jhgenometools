/*
  Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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

#ifndef PCK_COUNT_NODES_H
#define PCK_COUNT_NODES_H

#include <stdbool.h>

#include "core/stack-inlined.h"

#include "match/eis-voiditf.h"

typedef struct {
  unsigned long leaves, branching, lower, upper;
  unsigned int parentOffset;
  bool visited, on_branch;
} Nodecount;

GT_STACK_DECLARESTRUCT(Nodecount, 256UL);

void gt_pck_count_nodes_dfs(const FMindex *index,
                        unsigned long totallength,
                        unsigned int numofchars);

#endif
