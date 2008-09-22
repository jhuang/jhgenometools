/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/option.h"
#include "core/versionfunc.h"
#include "core/warning.h"
#include "extended/add_introns_stream.h"
#include "extended/gff3_in_stream.h"
#include "extended/splice_site_info_stream.h"
#include "extended/gtdatahelp.h"
#include "extended/seqid2file.h"
#include "tools/gt_splicesiteinfo.h"

typedef struct {
  GtStr *seqfile,
      *regionmapping;
  bool addintrons;
} SpliceSiteInfoArguments;

static OPrval parse_options(int *parsed_args,
                            SpliceSiteInfoArguments *arguments, int argc,
                            const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option;
  OPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("[option ...] [GFF3_file ...]", "Show information "
                         "about splice sites given in GFF3 files.");

  /* -seqfile and -regionmapping */
  seqid2file_options(op, arguments->seqfile, arguments->regionmapping);

  /* -addintrons */
  option = gt_option_new_bool("addintrons",
                           "add intron features between existing "
                           "exon features\n(before computing the information "
                           "to be shown)", &arguments->addintrons, false);
  gt_option_parser_add_option(op, option);

  /* parse */
  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_splicesiteinfo(int argc, const char **argv, GtError *err)
{
  GtNodeStream *gff3_in_stream = NULL,
               *add_introns_stream = NULL,
               *splice_site_info_stream = NULL;
  GtGenomeNode *gn;
  SpliceSiteInfoArguments arguments;
  GtRegionMapping *regionmapping;
  int parsed_args, had_err = 0;
  gt_error_check(err);

  /* option parsing */
  arguments.seqfile = gt_str_new();
  arguments.regionmapping = gt_str_new();
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      gt_str_delete(arguments.regionmapping);
      gt_str_delete(arguments.seqfile);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      gt_str_delete(arguments.regionmapping);
      gt_str_delete(arguments.seqfile);
      return 0;
  }

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                    argv + parsed_args,
                                                    false, false);

    /* create region mapping */
    regionmapping = seqid2file_regionmapping_new(arguments.seqfile,
                                                 arguments.regionmapping, err);
    if (!regionmapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create addintrons stream (if necessary) */
    if (arguments.addintrons)
      add_introns_stream = gt_add_introns_stream_new(gff3_in_stream);

    /* create extract feature stream */
    splice_site_info_stream = splice_site_info_stream_new(arguments.addintrons
                                                          ? add_introns_stream
                                                          : gff3_in_stream,
                                                          regionmapping);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = gt_node_stream_next(splice_site_info_stream, &gn,
                                               err)) && gn) {
      gt_genome_node_rec_delete(gn);
    }
  }

  if (!had_err) {
    if (!splice_site_info_stream_show(splice_site_info_stream)) {
      gt_warning("input file(s) contained no intron, use option -addintrons to "
                 "add introns automatically");
    }
  }

  /* free */
  gt_node_stream_delete(splice_site_info_stream);
  gt_node_stream_delete(add_introns_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_str_delete(arguments.regionmapping);
  gt_str_delete(arguments.seqfile);

  return had_err;
}
