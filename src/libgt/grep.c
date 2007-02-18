/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <sys/types.h>
#include <assert.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include "ensure.h"
#include "error.h"
#include "grep.h"

static void grep_error(int errcode, regex_t *matcher, Error *err)
{
  char sbuf[BUFSIZ], *buf;
  size_t bufsize;
  error_check(err);
  bufsize = regerror(errcode, matcher, NULL, 0);
  buf = malloc(bufsize);
  (void) regerror(errcode, matcher, buf ? buf : sbuf, buf ? bufsize : BUFSIZ);
  error_set(err, "grep(): %s", buf ? buf : sbuf);
  free(buf);
}

int grep(bool *match, const char *pattern, const char *line, Error *err)
{
  regex_t matcher;
  int rval, has_err = 0;
  error_check(err);
  assert(pattern && line);
  if ((rval = regcomp(&matcher, pattern, REG_EXTENDED | REG_NOSUB))) {
    grep_error(rval, &matcher, err);
    has_err = -1;
  }
  if (!has_err) {
    rval = regexec(&matcher, line, 0, NULL, 0);
    if (rval && rval != REG_NOMATCH) {
      grep_error(rval, &matcher, err);
      has_err = -1;
    }
  }
  regfree(&matcher);
  if (!has_err) {
    if (rval)
      *match = false;
    else
      *match = true;
  }
  return has_err;
}

int grep_unit_test(void)
{
  int has_err;
  bool match;

  has_err = grep(&match, "a", "a", NULL);
  ensure(!has_err);
  ensure(match);

  has_err = grep(&match, "b", "a", NULL);
  ensure(!has_err);
  ensure(!match);

  has_err =  grep(&match, "aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  ensure(!has_err);
  ensure(match);

  has_err = grep(&match, "aba", "wenbapzbpqSayhzzaBaZZqyghaAAahhaA", NULL);
  ensure(!has_err);
  ensure(!match);

  has_err = grep(&match, "^aba", "abawenbapzbpqSayhzzZZqyghaAAahhaA", NULL);
  ensure(!has_err);
  ensure(match);

  has_err = grep(&match, "^aba", "wenbapzbpqSayhzzabaZZqyghaAAahhaA", NULL);
  ensure(!has_err);
  ensure(!match);

  return EXIT_SUCCESS;
}
