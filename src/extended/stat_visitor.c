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
#include "core/disc_distri.h"
#include "core/unused_api.h"
#include "extended/node_visitor_rep.h"
#include "extended/stat_visitor.h"

struct GtStatVisitor {
  const GtNodeVisitor parent_instance;
  unsigned long number_of_sequence_regions,
                number_of_multi_features,
                number_of_genes,
                number_of_protein_coding_genes,
                number_of_mRNAs,
                number_of_exons,
                number_of_CDSs,
                number_of_LTR_retrotransposons,
                exon_number_for_distri,
                cds_length_for_distri;
  unsigned long long total_length_of_sequence_regions;
  GtDiscDistri *gene_length_distribution,
               *gene_score_distribution,
               *exon_length_distribution,
               *exon_number_distribution,
               *intron_length_distribution,
               *cds_length_distribution;
};

#define stat_visitor_cast(GV)\
        gt_node_visitor_cast(gt_stat_visitor_class(), GV)

static void stat_visitor_free(GtNodeVisitor *nv)
{
  GtStatVisitor *sv = stat_visitor_cast(nv);
  gt_disc_distri_delete(sv->cds_length_distribution);
  gt_disc_distri_delete(sv->intron_length_distribution);
  gt_disc_distri_delete(sv->exon_number_distribution);
  gt_disc_distri_delete(sv->exon_length_distribution);
  gt_disc_distri_delete(sv->gene_score_distribution);
  gt_disc_distri_delete(sv->gene_length_distribution);
}

static int add_exon_or_cds_number(GtGenomeNode *gn, void *data,
                                  GT_UNUSED GtError *err)
{
  GtStatVisitor *sv = (GtStatVisitor*) data;
  GtFeatureNode *fn = (GtFeatureNode*) gn;
  gt_error_check(err);
  gt_assert(sv && fn);
  if (gt_feature_node_has_type(fn, gt_ft_exon))
    sv->exon_number_for_distri++;
  else if (gt_feature_node_has_type(fn, gt_ft_CDS)) {
    GtRange range = gt_genome_node_get_range(gn);
    sv->cds_length_for_distri += gt_range_length(&range);
  }
  return 0;
}

static void compute_type_statistics(GtFeatureNode *fn, GtStatVisitor *sv)
{
  GtRange range;
  gt_assert(fn && sv);
  if (gt_feature_node_has_type(fn, gt_ft_gene)) {
    sv->number_of_genes++;
    if (gt_feature_node_has_CDS(fn))
      sv->number_of_protein_coding_genes++;
    if (sv->gene_length_distribution) {
      range = gt_genome_node_get_range((GtGenomeNode*) fn);
      gt_disc_distri_add(sv->gene_length_distribution, gt_range_length(&range));
    }
    if (sv->gene_score_distribution) {
      gt_disc_distri_add(sv->gene_score_distribution,
                         gt_feature_node_get_score(fn) * 100.0);
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_mRNA)) {
    sv->number_of_mRNAs++;
  }
  else if (gt_feature_node_has_type(fn, gt_ft_exon)) {
    sv->number_of_exons++;
    if (sv->exon_length_distribution) {
      range = gt_genome_node_get_range((GtGenomeNode*) fn);
      gt_disc_distri_add(sv->exon_length_distribution,
                         gt_range_length(&range));
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_CDS)) {
    sv->number_of_CDSs++;
  }
  else if (gt_feature_node_has_type(fn, gt_ft_intron)) {
    if (sv->intron_length_distribution) {
      range = gt_genome_node_get_range((GtGenomeNode*) fn);
      gt_disc_distri_add(sv->intron_length_distribution,
                         gt_range_length(&range));
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_LTR_retrotransposon)) {
    sv->number_of_LTR_retrotransposons++;
  }
}

static int compute_statistics(GtGenomeNode *gn, void *data, GtError *err)
{
  GtStatVisitor *sv;
  GtFeatureNode *fn;
  int rval;
  gt_error_check(err);
  gt_assert(data);
  sv = (GtStatVisitor*) data;
  fn = (GtFeatureNode*) gn;
  if (gt_feature_node_is_multi(fn) &&
      gt_feature_node_get_multi_representative(fn) == fn) {
    sv->number_of_multi_features++;
  }
  compute_type_statistics(fn, sv);
  if (sv->exon_number_distribution || sv->cds_length_distribution) {
    sv->exon_number_for_distri = 0;
    sv->cds_length_for_distri = 0;
    rval = gt_genome_node_traverse_direct_children(gn, sv,
                                                   add_exon_or_cds_number, err);
    gt_assert(!rval); /* add_exon_or_cds_number() is sane */
    if (sv->exon_number_distribution && sv->exon_number_for_distri) {
      gt_disc_distri_add(sv->exon_number_distribution,
                         sv->exon_number_for_distri);
    }
    if (sv->cds_length_distribution && sv->cds_length_for_distri) {
      gt_disc_distri_add(sv->cds_length_distribution,
                         sv->cds_length_for_distri);
    }
  }
  return 0;
}

static int stat_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                     GtError *err)
{
  GtStatVisitor *sv;
  gt_error_check(err);
  sv = stat_visitor_cast(nv);
  return gt_genome_node_traverse_children((GtGenomeNode*) fn, sv,
                                          compute_statistics, false, err);
}

static int stat_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                    GT_UNUSED GtError *err)
{
  GtStatVisitor *sv;
  GtRange range;
  gt_error_check(err);
  sv = stat_visitor_cast(nv);
  sv->number_of_sequence_regions++;
  range = gt_genome_node_get_range((GtGenomeNode*) rn);
  sv->total_length_of_sequence_regions += gt_range_length(&range);
  return 0;
}

const GtNodeVisitorClass* gt_stat_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtStatVisitor),
                                    stat_visitor_free,
                                    NULL,
                                    stat_visitor_feature_node,
                                    stat_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_stat_visitor_new(bool gene_length_distri,
                                   bool gene_score_distri,
                                   bool exon_length_distri,
                                   bool exon_number_distri,
                                   bool intron_length_distri,
                                   bool cds_length_distri)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_stat_visitor_class());
  GtStatVisitor *sv = stat_visitor_cast(nv);
  if (gene_length_distri)
    sv->gene_length_distribution = gt_disc_distri_new();
  if (gene_score_distri)
    sv->gene_score_distribution = gt_disc_distri_new();
  if (exon_length_distri)
    sv->exon_length_distribution = gt_disc_distri_new();
  if (exon_number_distri)
    sv->exon_number_distribution = gt_disc_distri_new();
  if (intron_length_distri)
    sv->intron_length_distribution = gt_disc_distri_new();
  if (cds_length_distri)
    sv->cds_length_distribution = gt_disc_distri_new();
  return nv;
}

void gt_stat_visitor_show_stats(GtNodeVisitor *nv, GtFile *outfp)
{
  GtStatVisitor *sv = stat_visitor_cast(nv);
  if (sv->number_of_sequence_regions) {
    gt_file_xprintf(outfp, "sequence regions: %lu (total length: %llu)\n",
                    sv->number_of_sequence_regions,
                    sv->total_length_of_sequence_regions);
  }
  if (sv->number_of_multi_features) {
    gt_file_xprintf(outfp, "multi-features: %lu\n",
                    sv->number_of_multi_features);
  }
  if (sv->number_of_genes)
    gt_file_xprintf(outfp, "genes: %lu\n", sv->number_of_genes);
  if (sv->number_of_protein_coding_genes) {
    gt_file_xprintf(outfp, "protein-coding genes: %lu\n",
                    sv->number_of_protein_coding_genes);
  }
  if (sv->number_of_mRNAs)
    gt_file_xprintf(outfp, "mRNAs: %lu\n", sv->number_of_mRNAs);
  if (sv->number_of_exons)
    gt_file_xprintf(outfp, "exons: %lu\n", sv->number_of_exons);
  if (sv->number_of_CDSs)
    gt_file_xprintf(outfp, "CDSs: %lu\n", sv->number_of_CDSs);
  if (sv->number_of_LTR_retrotransposons) {
    gt_file_xprintf(outfp, "LTR_retrotransposons: %lu\n",
                    sv->number_of_LTR_retrotransposons);
  }
  if (sv->gene_length_distribution) {
    gt_file_xprintf(outfp, "gene length distribution:\n");
    gt_disc_distri_show(sv->gene_length_distribution, outfp);
  }
  if (sv->gene_score_distribution) {
    gt_file_xprintf(outfp, "gene score distribution:\n");
    gt_disc_distri_show(sv->gene_score_distribution, outfp);
  }
  if (sv->exon_length_distribution) {
    gt_file_xprintf(outfp, "exon length distribution:\n");
    gt_disc_distri_show(sv->exon_length_distribution, outfp);
  }
  if (sv->exon_number_distribution) {
    gt_file_xprintf(outfp, "exon number distribution:\n");
    gt_disc_distri_show(sv->exon_number_distribution, outfp);
  }
  if (sv->intron_length_distribution) {
    gt_file_xprintf(outfp, "intron length distribution:\n");
    gt_disc_distri_show(sv->intron_length_distribution, outfp);
  }
  if (sv->cds_length_distribution) {
    gt_file_xprintf(outfp, "CDS length distribution:\n");
    gt_disc_distri_show(sv->cds_length_distribution, outfp);
  }
}
