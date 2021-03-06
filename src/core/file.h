/*
  Copyright (c) 2005-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdlib.h>
#include "core/file_api.h"

typedef enum {
  GT_FILE_MODE_UNCOMPRESSED,
  GT_FILE_MODE_GZIP,
  GT_FILE_MODE_BZIP2
} GtFileMode;

/* Returns GFM_GZIP if file with <path> ends with '.gz', GFM_BZIP2 if it ends
   with '.bz2', and GFM_UNCOMPRESSED otherwise. */
GtFileMode  gt_file_mode_determine(const char *path);

/* Returns ".gz" if <mode> is GFM_GZIP, ".bz2" if <mode> is GFM_BZIP2, and ""
   otherwise. */
const char* gt_file_mode_suffix(GtFileMode mode);

/* Returns the length of the ``basename'' of <path>. That is, the length of path
   without '.gz' or '.bz2' suffixes. */
size_t      gt_file_basename_length(const char *path);

/* Create a new GtFile object and open the underlying file handle, returns
   NULL and sets <err> if the file <path> could not be opened. */
GtFile*     gt_file_open(GtFileMode, const char *path, const char *mode,
                         GtError*);

/* Create a new GtFile object and open the underlying file handle, abort if
   the file <path> does not exist. The <file_mode> has to be given
   explicitly. */
GtFile*     gt_file_xopen_file_mode(GtFileMode file_mode, const char *path,
                                    const char *mode);

/* Create a new GtFile object and open the underlying file handle. Aborts if
   the file <path> could not be opened. The GtFileMode is determined
   automatically via gt_file_mode_determine(path). */
GtFile*     gt_file_xopen(const char *path, const char *mode);

/* Create a new GtFile object from a normal file pointer. */
GtFile*     gt_file_new_from_fileptr(FILE*);

GtFileMode  gt_file_mode(const GtFile*);

/* Return next character from <file> of EOF, if end-of-file is reached. */
int         gt_file_xfgetc(GtFile *file);

/* Unget character <c> to <file> (which obviously cannot be NULL).
   Can only be used once at a time. */
void        gt_file_unget_char(GtFile *file, char c);

/* printf(3) for generic files */
void        gt_file_xprintf(GtFile*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));

void        gt_file_xfputc(int c, GtFile*);

/* Read up to <nbytes> and store result in <buf>, returns bytes read. */
int         gt_file_xread(GtFile*, void *buf, size_t nbytes);

/* Write <nbytes> from <buf> to given generic file. */
void        gt_file_xwrite(GtFile*, void *buf, size_t nbytes);

/* Rewind the file. */
void        gt_file_xrewind(GtFile*);

/* Destroy the file handle object, but do not close the underlying handle. */
void        gt_file_delete_without_handle(GtFile*);

#endif
