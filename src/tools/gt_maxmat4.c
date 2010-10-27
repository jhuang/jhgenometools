/*
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

#include "match/maxmat4.h"
#include "match/matchmode_api.h"
#include "tools/gt_maxmat4.h"
#include "core/progress_timer_api.h"
#include "core/showtime.h"

#include "match/esa-map.h"
#include "match/eis-voiditf.h"

typedef struct {
  bool mum,
       mumreference,
       maxmatch,
       bothdirections,                       /* computebothdirections */
       reversecomplement,                    /* onlyreversecomplementmatches */
       showstring,                           /* showstring */
       showreversepositions,                 /* showreversepositions */
       fourcolumn,                           /* is option -F on? */
       showsequencelengths,                  /* showsequencelengths */
       prebwt,
       bitparallelism,
       bitparallelismwithbittab,
       showbitparallelismfactor,
       showtime;
  Definedunsignedlong leastlength,           /* leastlength */
                      bitlength;             /* num_of_bits */
} GtMaxmat4Arguments;

static void* gt_maxmat4_arguments_new(void)
{
  GtMaxmat4Arguments *arguments;
  arguments = gt_malloc(sizeof (GtMaxmat4Arguments));
  return arguments;
}

static void gt_maxmat4_arguments_delete(void *tool_arguments)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  if (!arguments)
  {
    return;
  }
  gt_free(arguments);
}

static GtOptionParser* gt_maxmat4_option_parser_new(void *tool_arguments)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  GtOption *option_mum, *option_mumreference, *option_maxmatch, *option_l,
           *option_b, *option_r, *option_s, *option_c, *option_F, *option_L,
           *option_p, *option_bit, *option_bittab, *option_bitlen, *option_f,
           *option_showtime;
  GtOptionParser *op;

  gt_assert(arguments);

  op = gt_option_parser_new(
      "[options] reference-indexname query-file-names",
      "Find and output (to stdout) the positions and length of all "
      "sufficiently long (unique) maximal matches of a substring in "
      "reference-indexname (in format of packed index) and query-file-names. "
      "Match only the characters a, c, g, or t, "
      "they can be in upper or in lower case.");
  gt_option_parser_set_mailaddress(op,"<gt-users@genometools.org>");

  /* -mum */
  option_mum = gt_option_new_bool("mum",
      "compute MUMs, i.e. maximal matches that are unique "
      "in both sequences ",
      &arguments->mum, false);
  gt_option_parser_add_option(op, option_mum);

  /* -mumreference */
  option_mumreference = gt_option_new_bool("mumreference",
      "compute MUM-candidates, i.e. maximal matches that are unique "
      "in the subject-sequence but not necessarily in the query-sequence "
      "If none of options [mum, mumreference and maxmatch] are used, "
      "then the program will default to mumreference mode ",
      &arguments->mumreference, false);
  gt_option_parser_add_option(op, option_mumreference);

  /* -maxmatch */
  option_maxmatch = gt_option_new_bool("maxmatch",
      "compute all maximal matches, regardless of their uniqueness ",
      &arguments->maxmatch, false);
  gt_option_parser_add_option(op, option_maxmatch);

  /* -l for leastlength */
  option_l = gt_option_new_ulong_min("l",
      "set the minimum length of a match "
      "if not set, the default value is 20",
      &arguments->leastlength.valueunsignedlong,
      (unsigned long)20,(unsigned long) 1);
  gt_option_parser_add_option(op, option_l);

  /* -b for computebothdirections */
  option_b = gt_option_new_bool("b",
      "compute forward and reverse complement matches",
      &arguments->bothdirections, false);
  gt_option_parser_add_option(op, option_b);

  /* -r for onlyreversecomplementmatches */
  option_r = gt_option_new_bool("r",
      "only compute reverse complement matches",
      &arguments->reversecomplement, false);
  gt_option_parser_add_option(op, option_r);

  /* -s for showstring */
  option_s = gt_option_new_bool("s",
      "show the matching substrings",
      &arguments->showstring, false);
  gt_option_parser_add_option(op, option_s);

  /* -c for showreversepositions */
  option_c = gt_option_new_bool("c",
      "report the query-position of a reverse complement match "
      "relative to the original query sequence",
      &arguments->showreversepositions, false);
  gt_option_parser_add_option(op, option_c);

  /* -F for fourcolumn */
  option_F = gt_option_new_bool("F",
      "force 4 column output format that prepends every match line with "
      "the reference sequence identifier",
      &arguments->fourcolumn, false);
  gt_option_parser_add_option(op, option_F);

  /* -L for showsequencelengths */
  option_L = gt_option_new_bool("L",
      "show the length of the query sequence on the header line",
      &arguments->showsequencelengths, false);
  gt_option_parser_add_option(op, option_L);

  /* -p for prebwt */
  option_p = gt_option_new_bool("p",
      "use precomputed bwt-bounds files (file type: pbt), which can "
      "be precomputed by command 'gt prebwt'. If no pbt-file exists, it "
      "will be computed automatically in normal process",
      &arguments->prebwt, false);
  gt_option_parser_add_option(op, option_p);

  /* -bit for bitparallelism */
  option_bit = gt_option_new_bool("bit",
      "use bit parallelism",
      &arguments->bitparallelism, false);
  gt_option_parser_add_option(op, option_bit);

  /* -bittab for bitparallelism with bittab */
  option_bittab = gt_option_new_bool("bittab",
      "use bit parallelism with bittab",
      &arguments->bitparallelismwithbittab, false);
  gt_option_parser_add_option(op, option_bittab);

  /* -bitlen for bittab length */
  option_bitlen = gt_option_new_ulong_min("bitlen",
      "set custom length of bittab for option bittab "
      "if not set, the default value is 32",
      &arguments->bitlength.valueunsignedlong,
      (unsigned long)32,(unsigned long) 1);
  gt_option_parser_add_option(op, option_bitlen);

  /* -f for showbitparallelismfactor */
  option_f = gt_option_new_bool("f",
      "show bit parallelism factor for every query, \
      that is， 1/average offset",
      &arguments->showbitparallelismfactor, false);
  gt_option_parser_add_option(op, option_f);

  /* -showtime for showtime */
  option_showtime = gt_option_new_bool("showtime",
      "show the progress time",
      &arguments->showtime, false);
  gt_option_parser_add_option(op, option_showtime);

  /* option implications */
  gt_option_imply_either_2(option_c, option_b, option_r);
  gt_option_imply_either_2(option_f, option_bit, option_bittab);
  gt_option_imply(option_bitlen, option_bittab);

  /* option exclusions */
  gt_option_exclude(option_mum, option_mumreference);
  gt_option_exclude(option_mum, option_maxmatch);
  gt_option_exclude(option_mumreference, option_maxmatch);
  gt_option_exclude(option_b, option_r);
  gt_option_exclude(option_bit, option_bittab);

  /* set minimal arugments */
  gt_option_parser_set_min_args(op, (unsigned int)2);

  return op;
}

static int gt_maxmat4_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments,
                                      GtError *err)
{
  GT_UNUSED GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

static int gt_maxmat4_runner(GT_UNUSED int argc,
                             GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args,
                             void *tool_arguments,
                             GtError *err)
{
  GtMaxmat4Arguments *arguments = (GtMaxmat4Arguments*)tool_arguments;
  int arg = parsed_args;

  Suffixarray suffixarray;
  GtLogger *logger = NULL;
  bool haserr = false;
  const GtAlphabet *alphabet = NULL;
  unsigned long totallength;
  unsigned int mappedbits;
  int idx;

  void *packedindex = NULL;
  const Mbtab **mbtab = NULL;
  unsigned int maxdepth = 0;

  GtStr *referencefile = NULL;
  GtStrArray *queryfiles = NULL;

  GtProgressTimer *maxmat4progress = NULL;
  if (gt_showtime_enabled() || arguments->showtime)
  {
    maxmat4progress = \
    gt_progress_timer_new("finding maximal matches of some minimum length "
                       "between a reference sequence and a query-sequence");
  }

  /* init the referencefile and queryfiles */
  referencefile = gt_str_new_cstr(argv[arg]);
  queryfiles = gt_str_array_new();

  for (idx=arg+1; idx<argc; idx++)
  {
    gt_str_array_add_cstr(queryfiles, argv[idx]);
  }

  /* load packed index in memory */
  gt_error_check(err);
  gt_assert(tool_arguments);
  logger = gt_logger_new(false, GT_LOGGER_DEFLT_PREFIX, stdout);
  mappedbits = SARR_ESQTAB|SARR_SSPTAB|SARR_DESTAB|SARR_SDSTAB;
  /*
  if (arguments->prebwt)
  {
     mappedbits |= SARR_PBTTAB;
  }
  */

  /*
   * map suffixarray from referencefile, a referencefile contains maybe many
   * sequence, is this merged in a suffixarray?
   */
  if (gt_mapsuffixarray(&suffixarray,
                        mappedbits,
                        gt_str_get(referencefile),
                        logger,
                        err) != 0)
  {
    gt_assert(!haserr);
    /* gt_error_set(err, "error in reading reference file"); */
    haserr = true;
    totallength = 0;
  } else
  {
    alphabet = gt_encseq_alphabet(suffixarray.encseq);
    totallength = gt_encseq_total_length(suffixarray.encseq);
  }

  if (!haserr)
  {
    /* load packed index from suffxiarray */
    packedindex = gt_loadvoidBWTSeqForSA(gt_str_get(referencefile),
                                       arguments->prebwt,  /* bool withpckbwt */
                                       err);
    if (packedindex == NULL)
    {
      haserr = true;
    }

    if (arguments->prebwt)
    {
      mbtab = gt_bwtseq2mbtab(packedindex);
      maxdepth = gt_bwtseq2maxdepth(packedindex);
    }
  }

  if (!haserr)
  {
    /* set match mode */
    GtMatchmode matchmode = GT_MATCHMODE_MUMREFERENCE;
    if (arguments->mum)
    {
      matchmode = GT_MATCHMODE_MUM;
    }
    else if (arguments->mumreference)
    {
      matchmode = GT_MATCHMODE_MUMREFERENCE;
    }
    else if (arguments->maxmatch)
    {
      matchmode = GT_MATCHMODE_MAXMATCH;
    }

    if (gt_findmum(suffixarray.encseq,
                   packedindex,
                   /* if option prebwt is false -> mbtab is NULL and maxdepth
                    * is 0, but it doesn't matter because the both parameters
                    * won't be used later in this case */
                   mbtab,
                   maxdepth,
                   totallength,
                   alphabet,
                   queryfiles,
                   matchmode,
                   arguments->leastlength,
                   arguments->bothdirections,
                   arguments->reversecomplement,
                   arguments->showstring,
                   arguments->showreversepositions,
                   arguments->fourcolumn,
                   arguments->showsequencelengths,
                   arguments->prebwt,
                   arguments->bitparallelism,
                   arguments->bitparallelismwithbittab,
                   arguments->bitlength,
                   arguments->showbitparallelismfactor,
                   arguments->showtime,
                   maxmat4progress,
                   logger,
                   err) != 0)
    {
      haserr = true;
    }
  }

  if (packedindex != NULL)
  {
    gt_deletevoidBWTSeq(packedindex);
  }
  gt_freesuffixarray(&suffixarray);

  gt_logger_delete(logger);
  gt_str_delete(referencefile);
  gt_str_array_delete(queryfiles);

  if (maxmat4progress != NULL)
  {
    gt_progress_timer_start_new_state(maxmat4progress,NULL,stdout);
    gt_progress_timer_delete(maxmat4progress);
  }
  return haserr ? -1 : 0;
}

GtTool* gt_maxmat4(void)
{
  return gt_tool_new(gt_maxmat4_arguments_new,
                     gt_maxmat4_arguments_delete,
                     gt_maxmat4_option_parser_new,
                     gt_maxmat4_arguments_check,
                     gt_maxmat4_runner);
}
