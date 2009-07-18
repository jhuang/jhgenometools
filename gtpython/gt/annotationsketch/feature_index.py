#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

from gt.dlload import gtlib
from gt.core.array import Array
from gt.core.error import Error, gterror
from gt.core.gtrange import Range
from gt.core.str_array import StrArray
from gt.extended.feature_node import FeatureNode


class FeatureIndex:

    def __init__(self, *args):
        raise NotImplementedError, \
            'Please call the constructor of a FeatureIndex implementation.'

    def __del__(self):
        try:
            gtlib.gt_feature_index_delete(self.fi)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, FeatureIndex):
            raise TypeError, "argument must be a FeatureIndex"
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def get_features_for_seqid(self, seqid):
        rval = gtlib.gt_feature_index_get_features_for_seqid(self.fi,
                seqid)
        if rval:
            a = Array(rval, True)
            result = []
            for i in range(a.size()):
                fptr = gtlib.gt_genome_node_ref(a.get(i))
                result.append(FeatureNode.create_from_ptr(fptr))
            return result
        else:
            return None

    def add_gff3file(self, filename):
        err = Error()
        rval = gtlib.gt_feature_index_add_gff3file(self.fi, filename,
                err)
        if rval != 0:
            gterror(err)

    def get_first_seqid(self):
        return gtlib.gt_feature_index_get_first_seqid(self.fi)

    def get_seqids(self):
        result = []
        stra = StrArray(gtlib.gt_feature_index_get_seqids(self.fi))
        for i in range(stra.size()):
            result.append(stra.get(i))
        return result

    def get_range_for_seqid(self, seqid):
        from ctypes import byref
        if gtlib.gt_feature_index_has_seqid(self.fi, seqid) == 0:
            gterror("feature_index does not contain seqid")
        range = Range()
        gtlib.gt_feature_index_get_range_for_seqid(self.fi, byref(range),
                seqid)
        return range

    def get_features_for_range(self, start, end, seqid):
        from ctypes import byref
        a = Array.create()
        err = Error()
        rng = Range(start, end)
        rval = gtlib.gt_feature_index_get_features_for_range(self.fi, a,
                seqid, byref(rng), err)
        if rval != 0:
            gterror(err)
        result = []
        for i in range(a.size()):
            fptr = gtlib.gt_genome_node_ref(a.get(i))
            result.append(FeatureNode.create_from_ptr(fptr))
        return result

    def register(cls, gtlib):
        from ctypes import c_char_p, c_void_p, c_int, POINTER
        gtlib.gt_feature_index_get_features_for_seqid.restype = c_void_p
        gtlib.gt_feature_index_add_gff3file.argtypes = [c_void_p,
                c_char_p, Error]
        gtlib.gt_feature_index_get_first_seqid.restype = c_char_p
        gtlib.gt_feature_index_get_seqids.restype = c_void_p
        gtlib.gt_feature_index_has_seqid.argtypes = [c_void_p, c_char_p]
        gtlib.gt_feature_index_get_range_for_seqid.argtypes = [c_void_p,
                POINTER(Range), c_char_p]
        gtlib.gt_feature_index_get_features_for_range.argtypes = [c_void_p,
                Array, c_char_p, POINTER(Range), Error]

    register = classmethod(register)


class FeatureIndexMemory(FeatureIndex):

    def __init__(self):
        self.fi = gtlib.gt_feature_index_memory_new()
        self._as_parameter_ = self.fi

    def from_param(cls, obj):
        if not isinstance(obj, FeatureIndexMemory):
            raise TypeError, "argument must be a FeatureIndexMemory"
        return obj._as_parameter_

    from_param = classmethod(from_param)

