/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_REMAINSORT_H
#define SFX_REMAINSORT_H

#include "core/error_api.h"
#include "core/codetype.h"
#include "core/readmode.h"
#include "core/defined-types.h"

#include "bcktab.h"
#include "compressedtab.h"
#include "sfx-suffixgetset.h"

typedef struct Rmnsufinfo Rmnsufinfo;

Rmnsufinfo *gt_rmnsufinfo_new(GtSuffixsortspace *suffixsortspace,
                              int mmapfiledesc,
                              GtStr *mmapfilename,
                              const GtEncseq *encseq,
                              Bcktab *bcktab,
                              GtCodetype maxcode,
                              unsigned int numofchars,
                              unsigned int prefixlength,
                              GtReadmode readmode,
                              unsigned long partwidth,
                              bool hashexceptions,
                              bool absoluteinversesuftab);

void gt_rmnsufinfo_addunsortedrange(Rmnsufinfo *rmnsufinfo,
                                    unsigned long left,
                                    unsigned long right,
                                    unsigned long depth);

void gt_rmnsufinfo_bcktab2firstlevelintervals(Rmnsufinfo *rmnsufinfo );

Compressedtable *gt_rmnsufinfo_delete(Rmnsufinfo **rmnsufinfoptr,
                                      bool withlcptab);

#endif
