/*
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/assert_api.h"
#include "core/chardef.h"
#include "maxmat4-initeqsvec.h"

void gt_maxmat4_initeqsvectorrev(GtBittab **C,
                      unsigned long eqslen,
                      unsigned long bitlen,
                      unsigned long leastlength,
                      const GtUchar *pattern,
                      unsigned long patternlength)
{
  const GtUchar *pptr;
  unsigned long i;
  long position;  /* have to use long, unsigned long doesn't work */
    
  for (i = 0; i < eqslen; i++) {
    gt_bittab_unset(C[i*leastlength]);
  }

  for (position = patternlength-1; position >= 0; position--) {
		pptr = pattern+position;
		gt_assert(*pptr != (GtUchar) SEPARATOR);
    if (*pptr != (GtUchar) WILDCARD)
    {
			gt_bittab_set_bit(C[(unsigned long) (*pptr)*leastlength], position + bitlen - patternlength);
    }   
  } 
}
