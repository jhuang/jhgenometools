/*
  Copyright (c) 2004-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "gth/gthoutput.h"
#include "gth/gthsortags.h"
#include "gth/ags.h"
#include "gth/pgl.h"

static int compareAGSs(const void *dataA, const void *dataB)
{
  GthAGS *agsA = *(GthAGS**) dataA;
  GthAGS *agsB = *(GthAGS**) dataB;
  gt_assert(agsA->overallscore != UNDEFHIGHPRECPROBTYPE);
  gt_assert(agsB->overallscore != UNDEFHIGHPRECPROBTYPE);
  if (agsA->overallscore < agsB->overallscore)
    return 1;
  if (agsA->overallscore > agsB->overallscore)
    return -1;
  return 0;
}

static void determineAGSscore(GthAGS *ags, double sortagswf)
{
  HIGHPRECPROBTYPE average_exon_score              = 0.0,
                   average_splice_site_probability = 0.0;
  unsigned long i, numofexons = gt_array_size(ags->exons);
  Splicesiteprob *splicesiteprob;
  gt_assert(numofexons > 0);

  if (numofexons == 1) {
    /* if the AGS contains only one exon, the exonscore equals the overall score
    */
    ags->overallscore = ((ExoninfoAGS*) gt_array_get_first(ags->exons))
                        ->exonscore;
  }
  else {
    /* compute weighted mean of the average exon score and the average splice
       site probability */
    for (i = 0; i < numofexons; i++) {
      /* sum them up */
      average_exon_score += ((ExoninfoAGS*) gt_array_get(ags->exons, i))
                            ->exonscore;
      if (i > 0) {
        splicesiteprob = gt_array_get(ags->splicesiteprobs, i-1);
        average_splice_site_probability += splicesiteprob->donorsiteprob;
        average_splice_site_probability += splicesiteprob->acceptorsiteprob;
      }
    }
    /* calc average */
    average_exon_score              /= numofexons;
    average_splice_site_probability /= (2 * (numofexons - 1));

    /* calc weighted mean */
    ags->overallscore = (sortagswf * average_exon_score +
                         average_splice_site_probability) /
                        (sortagswf + 1.0);

  }
}

static void sortAGSsinsinglePGL(GthPGL *pgl, double sortagswf)
{
  unsigned long i;

  /* determine AGS scores */
  for (i = 0; i < gt_array_size(pgl->assemblies); i++)
    determineAGSscore(*(GthAGS**) gt_array_get(pgl->assemblies, i), sortagswf);

  /* sort AGSs according to the scores */
  qsort(gt_array_get_space(pgl->assemblies), gt_array_size(pgl->assemblies),
        sizeof (GthAGS*), compareAGSs);
}

void sortAGSs(GtArray *pgls, double sortagswf)
{
  unsigned long i;
  for (i = 0; i < gt_array_size(pgls); i++)
    sortAGSsinsinglePGL(*(GthPGL**) gt_array_get(pgls, i), sortagswf);
}