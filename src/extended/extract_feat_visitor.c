/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/codon_iterator_simple.h"
#include "core/fasta.h"
#include "core/symbol.h"
#include "core/translator.h"
#include "extended/extract_feat_sequence.h"
#include "extended/extract_feat_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_rep.h"

struct GtExtractFeatVisitor {
  const GtNodeVisitor parent_instance;
  const char *type;
  bool join,
       translate;
  unsigned long fastaseq_counter,
                width;
  GtRegionMapping *region_mapping;
  GtFile *outfp;
};

#define gt_extract_feat_visitor_cast(GV)\
        gt_node_visitor_cast(gt_extract_feat_visitor_class(), GV)

static void extract_feat_visitor_free(GtNodeVisitor *nv)
{
  GtExtractFeatVisitor *extract_feat_visitor = gt_extract_feat_visitor_cast(nv);
  gt_assert(extract_feat_visitor);
  gt_region_mapping_delete(extract_feat_visitor->region_mapping);
}

static void construct_description(GtStr *description, const char *type,
                                  unsigned long counter, bool join,
                                  bool translate)
{
  gt_assert(!gt_str_length(description));
  gt_str_append_cstr(description, type);
  gt_str_append_char(description, '_');
  gt_str_append_ulong(description, counter);
  if (join)
    gt_str_append_cstr(description, " (joined)");
  if (translate)
    gt_str_append_cstr(description, " (translated)");
}

static void show_entry(GtStr *description, GtStr *sequence, bool translate,
                       unsigned long width, GtFile *outfp)
{
  if (translate) {
    char translated;
    int rval;
    unsigned int frame;
    GtStr *protein = gt_str_new();

    GtCodonIterator *ci = gt_codon_iterator_simple_new(gt_str_get(sequence),
                                                       gt_str_length(sequence),
                                                       NULL);
    GtTranslator* tr = gt_translator_new(ci);
    rval = gt_translator_next(tr, &translated, &frame, NULL);
    while (!rval && translated) {
      gt_str_append_char(protein, translated);
      rval = gt_translator_next(tr, &translated, &frame, NULL);
    }
    gt_fasta_show_entry(gt_str_get(description), gt_str_get(protein),
                        gt_str_length(protein), width, outfp);
    gt_str_delete(protein);
    gt_translator_delete(tr);
    gt_codon_iterator_delete(ci);
  }
  else {
    gt_fasta_show_entry(gt_str_get(description), gt_str_get(sequence),
                        gt_str_length(sequence), width, outfp);
  }
}

static int extract_feat_visitor_feature_node(GtNodeVisitor *nv,
                                             GtFeatureNode *fn, GtError *err)
{
  GtExtractFeatVisitor *efv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  GtStr *description,
      *sequence;
  int had_err = 0;
  gt_error_check(err);
  efv = gt_extract_feat_visitor_cast(nv);
  gt_assert(efv->region_mapping);
  fni = gt_feature_node_iterator_new(fn);
  description = gt_str_new();
  sequence = gt_str_new();
  while (!had_err && (child = gt_feature_node_iterator_next(fni))) {
    if (gt_extract_feat_sequence(sequence, (GtGenomeNode*) child, efv->type,
                                 efv->join, efv->region_mapping, err)) {
      had_err = -1;
    }

    if (!had_err && gt_str_length(sequence)) {
      efv->fastaseq_counter++;
      construct_description(description, efv->type, efv->fastaseq_counter,
                            efv->join, efv->translate);
      show_entry(description, sequence, efv->translate, efv->width, efv->outfp);
      gt_str_reset(description);
      gt_str_reset(sequence);
    }
  }
  gt_str_delete(sequence);
  gt_str_delete(description);
  gt_feature_node_iterator_delete(fni);
  return had_err;
}

const GtNodeVisitorClass* gt_extract_feat_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtExtractFeatVisitor),
                                    extract_feat_visitor_free,
                                    NULL,
                                    extract_feat_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_extract_feat_visitor_new(GtRegionMapping *rm,
                                           const char *type, bool join,
                                           bool translate, unsigned long width,
                                           GtFile *outfp)
{
  GtNodeVisitor *nv;
  GtExtractFeatVisitor *efv;
  gt_assert(rm);
  nv = gt_node_visitor_create(gt_extract_feat_visitor_class());
  efv= gt_extract_feat_visitor_cast(nv);
  efv->type = gt_symbol(type);
  efv->join = join;
  efv->translate = translate;
  efv->fastaseq_counter = 0;
  efv->region_mapping = rm;
  efv->width = width;
  efv->outfp = outfp;
  return nv;
}
