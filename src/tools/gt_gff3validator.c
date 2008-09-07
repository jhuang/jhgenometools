/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/ma.h"
#include "extended/gff3_in_stream.h"
#include "extended/feature_type_factory_obo.h"
#include "tools/gt_gff3validator.h"

typedef struct {
  Str *typecheck;
} GFF3ValidatorArguments;

static void* gt_gff3validator_arguments_new(void)
{
  GFF3ValidatorArguments *arguments = ma_calloc(1, sizeof *arguments);
  arguments->typecheck = str_new();
  return arguments;
}

static void gt_gff3validator_arguments_delete(void *tool_arguments)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  if (!arguments) return;
  str_delete(arguments->typecheck);
  ma_free(arguments);
}

static OptionParser* gt_gff3validator_option_parser_new(void *tool_arguments)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *option;
  assert(arguments);

  /* init */
  op = option_parser_new("[option ...] [GFF3_file ...]",
                         "Strictly validate given GFF3 files.");

  /* -typecheck */
  option = option_new_filename("typecheck", "check GFF3 types against \"id\" "
                               "and \"name\" tags in given OBO file",
                               arguments->typecheck);
  option_parser_add_option(op, option);

  return op;
}

static int gt_gff3validator_runner(int argc, const char **argv, int parsed_args,
                                   void *tool_arguments, GT_Error *err)
{
  GFF3ValidatorArguments *arguments = tool_arguments;
  FeatureTypeFactory *ftf = NULL;
  GenomeStream *gff3_in_stream;
  GT_GenomeNode *gn;
  int had_err = 0;

  gt_error_check(err);
  assert(arguments);

  /* create a gff3 input stream */
  gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                               argv + parsed_args,
                                               false, true);

  /* set different type checker if necessary */
  if (str_length(arguments->typecheck)) {
    if (!(ftf = feature_type_factory_obo_new(str_get(arguments->typecheck),
                                             err))) {
        had_err = -1;
    }
    if (!had_err)
      gff3_in_stream_set_feature_type_factory(gff3_in_stream, ftf);
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err) {
    while (!(had_err = genome_stream_next_tree(gff3_in_stream, &gn, err)) &&
           gn) {
      gt_genome_node_rec_delete(gn);
    }
  }

  if (!had_err)
    printf("input is valid GFF3\n");

  /* free */
  genome_stream_delete(gff3_in_stream);
  feature_type_factory_delete(ftf);

  return had_err;
}

Tool* gt_gff3validator(void)
{
  return tool_new(gt_gff3validator_arguments_new,
                  gt_gff3validator_arguments_delete,
                  gt_gff3validator_option_parser_new,
                  NULL,
                  gt_gff3validator_runner);
}
