/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SA_FILTER_H
#define SA_FILTER_H

#include "core/option.h"
#include "gth/sa.h"

/* a filter for spliced alignments */
typedef struct GthSAFilter GthSAFilter;

GthSAFilter* gth_sa_filter_new(void);
void         gth_sa_filter_register_options(GtOptionParser*, GthSAFilter*,
                                        bool extended_options);
void         gth_sa_filter_set_min_alignmentscore(GthSAFilter*, double);
void         gth_sa_filter_set_max_alignmentscore(GthSAFilter*, double);
void         gth_sa_filter_set_min_coverage(GthSAFilter*, double);
void         gth_sa_filter_set_max_coverage(GthSAFilter*, double);
double       gth_sa_filter_get_min_alignmentscore(const GthSAFilter*);
double       gth_sa_filter_get_max_alignmentscore(const GthSAFilter*);
double       gth_sa_filter_get_min_coverage(const GthSAFilter*);
double       gth_sa_filter_get_max_coverage(const GthSAFilter*);
bool         gth_sa_filter_filter_sa(const GthSAFilter*, GthSA*);
void         gth_sa_filter_delete(GthSAFilter*);

#endif
