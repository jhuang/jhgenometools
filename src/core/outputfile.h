/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#include "core/option.h"

#define GT_FORCE_OPT_CSTR  "force"

typedef struct GtOutputFileInfo GtOutputFileInfo;

GtOutputFileInfo* gt_outputfileinfo_new(void);
void              gt_outputfile_register_options(GtOptionParser*,
                                                 GtFile **outfp,
                                                 GtOutputFileInfo*);
void              gt_outputfileinfo_delete(GtOutputFileInfo*);

/* Helper funtion for (rare) tools which do not use the full <GtOutputFileInfo>
   (usually if directories are involved). */
GtFile*           gt_outputfile_xopen_forcecheck(const char *path,
                                                 const char *mode, bool force,
                                                 GtError *err);

#endif
