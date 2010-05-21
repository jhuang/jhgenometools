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

#ifndef MATCHMODE_API_H
#define MATCHMODE_API_H

#include "core/error_api.h"

typedef enum
{
  GT_MATCHMODE_MUM = 0,
  GT_MATCHMODE_MUMREFERENCE,
  GT_MATCHMODE_MAXMATCH
} GtMatchmode;

/* Returns the descriptive string for the matchmode <matchmode>. */
//const char* gt_matchmode_show(GtMatchmode matchmode);
/* Returns the <GtMatchmode> for the description <string>, which must be one
   of "mum","mumreference" or "maxmatch. If <string> does not equal any of them,
   -1 is returned and <err> is set accordingly. */
//int         gt_matchmode_parse(const char *string, GtError *err);

#endif
