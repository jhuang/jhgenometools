/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/hashmap.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/merge_feature_visitor.h"
#include "extended/node_visitor_rep.h"

struct GtMergeFeatureVisitor {
  const GtNodeVisitor parent_instance;
  GtGenomeNode *current_tree;
  GtHashmap *hm; /* type -> previous node */
  GtArray *nodes_to_remove;
};

#define gt_merge_feature_visitor_cast(GV)\
        gt_node_visitor_cast(gt_merge_feature_visitor_class(), GV)

static void merge_feature_visitor_free(GtNodeVisitor *nv)
{
  GtMergeFeatureVisitor *merge_feature_visitor =
    gt_merge_feature_visitor_cast(nv);
  gt_assert(merge_feature_visitor);
  gt_hashmap_delete(merge_feature_visitor->hm);
  gt_array_delete(merge_feature_visitor->nodes_to_remove);
}

static int mergefeat_in_children(GtGenomeNode *gn, void *data,
                                 GT_UNUSED GtError *err)
{
  GtMergeFeatureVisitor *v = (GtMergeFeatureVisitor*) data;
  GtFeatureNode *previous_feature, *current_feature;
  GtRange previous_range, current_range;
  gt_error_check(err);
  current_feature = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(current_feature);
  if ((previous_feature =
        gt_hashmap_get(v->hm, gt_feature_node_get_type(current_feature)))) {
    /* previous feature found -> check if merging is necessary */
    gt_assert(gt_feature_node_get_type(previous_feature) ==
              gt_feature_node_get_type(current_feature));
    previous_range = gt_genome_node_get_range((GtGenomeNode*)
                                              previous_feature);
    current_range = gt_genome_node_get_range((GtGenomeNode*) current_feature);
    /* sorted */
    gt_assert(gt_range_compare(&previous_range, &current_range) <= 0);
    if (previous_range.end + 1 == current_range.start) {
      /* merge nodes */
      gt_feature_node_set_end(previous_feature, current_range.end);
      /* XXX: compute average score ? */
      gt_feature_node_unset_score(previous_feature);
      gt_assert(!gt_genome_node_number_of_children((GtGenomeNode*)
                                                current_feature));
      gt_array_add(v->nodes_to_remove, current_feature);
    }
    /* remove previous feature */
    gt_hashmap_remove(v->hm, gt_feature_node_get_type(previous_feature));
  }
  /* add current feature */
  gt_hashmap_add(v->hm, (char*) gt_feature_node_get_type(current_feature),
                 current_feature);
  return 0;
}

static int mergefeat_if_necessary(GtGenomeNode *gn, void *data, GtError *err)
{
  GtMergeFeatureVisitor *v = (GtMergeFeatureVisitor*) data;
  GtFeatureNode *fn;
  gt_error_check(err);
  fn = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(fn);
  v->current_tree = gn;
  gt_hashmap_reset(v->hm);
  return gt_genome_node_traverse_direct_children(gn, v, mergefeat_in_children,
                                                 err);
}

static int merge_feature_visitor_feature_node(GtNodeVisitor *nv,
                                          GtFeatureNode *fn, GtError *err)
{
  GtMergeFeatureVisitor *v;
  GtGenomeNode *leaf;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  v = gt_merge_feature_visitor_cast(nv);
  gt_array_reset(v->nodes_to_remove);
  had_err = gt_genome_node_traverse_children((GtGenomeNode*) fn, v,
                                          mergefeat_if_necessary, false, err);
  if (!had_err) {
    for (i = 0; i < gt_array_size(v->nodes_to_remove); i++) {
      leaf = *(GtGenomeNode**) gt_array_get(v->nodes_to_remove, i);
      gt_genome_node_remove_leaf((GtGenomeNode*) fn, leaf);
      gt_genome_node_delete(gt_feature_node_cast(leaf));
    }
  }
  return had_err;
}

const GtNodeVisitorClass* gt_merge_feature_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtMergeFeatureVisitor),
                                    merge_feature_visitor_free,
                                    NULL,
                                    merge_feature_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_merge_feature_visitor_new(void)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_merge_feature_visitor_class());
  GtMergeFeatureVisitor *merge_feature_visitor =
    gt_merge_feature_visitor_cast(nv);
  merge_feature_visitor->hm = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  merge_feature_visitor->nodes_to_remove = gt_array_new(sizeof (GtGenomeNode*));
  return nv;
}
