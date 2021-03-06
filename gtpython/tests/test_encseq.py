#!/usr/bin/python
# -*- coding: utf-8 -*-

from gt.core.encseq import *
from gt.core.error import GTError
from gt.core import readmode
import unittest
import tempfile
import sys, os
import string

class EncodedsequenceTest(unittest.TestCase):

    def setUp(self):
        self.dnaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.dseq1 = "agtccagctgtcagctagcgggcccgatgatatttt"
        self.dseq2 = "gtgctgtac"
        self.dnaseqfile.write(">seq1\n"+self.dseq1+"\n")
        self.dnaseqfile.write(">seq2\n"+self.dseq2+"\n")
        self.aaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.aaseq1 = "MVHFTAEEKAAVTSLWSKMNVEEAGGEALG"
        self.aaseq2 = "KMNAVE"
        self.aaseqfile.write(">seq1\n"+self.aaseq1+"\n")
        self.aaseqfile.write(">seq2\n"+self.aaseq2+"\n")
        self.dnaseqfile.close()
        self.aaseqfile.close()
        self.idxsuffixes = ['esq','des','ssp','sds','al1']

    def tearDown(self):
        os.unlink(self.dnaseqfile.name)
        os.unlink(self.aaseqfile.name)

    def create_es(self, indexname):
        ee = EncseqEncoder()
        return ee.encode([self.dnaseqfile.name], indexname)                

    def create_es_protein(self, indexname):
        ee = EncseqEncoder()
        return ee.encode([self.aaseqfile.name], indexname)
	        
    def create_mem(self):
        a = Alphabet.create_dna()
        eb = EncseqBuilder(a)
        eb.enable_description_support()
        eb.enable_multiseq_support()
        eb.add_string(self.dseq1, 'seq1')
        eb.add_string(self.dseq2, 'seq2')
        return eb.build() 

    def create_mem_protein(self):
        a = Alphabet.create_protein()
        eb = EncseqBuilder(a)
        eb.enable_description_support()
        eb.enable_multiseq_support()
        eb.add_string(self.aaseq1, 'seq1')
        eb.add_string(self.aaseq2, 'seq2')
        return eb.build() 
        
    def delete_idx(self, indexname):
        for suf in self.idxsuffixes:
            os.unlink(indexname+"."+suf)

    def test_create_new(self):
        val = self.create_es("foo")
        for suf in self.idxsuffixes:
            self.assertTrue(os.path.isfile("foo."+suf))
        self.delete_idx("foo")

    def test_create_mapped(self):
        self.create_es("foo_mapped")
        el = EncseqLoader() 
        es = el.load("foo_mapped")
        self.assertNotEqual(es, None)
        self.delete_idx("foo_mapped")

    def test_map_fail(self):
        el = EncseqLoader()
        self.assertRaises(IOError, el.load, "foo_fail")

    def test_dna(self):
        self.create_es("foo")
        el = EncseqLoader()
        es = el.load("foo")
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(es, self.dseq1, self.dseq2)
        self.run_test_num_seqs(es)
        self.run_test_num_files(es)
        self.run_test_seq_length(es)
        self.run_test_seq_startpos(es)
        self.run_test_seq_substr_encoded(es, self.dseq1, self.dseq2)
        self.run_test_seq_substr_plain(es, self.dseq1, self.dseq2)
        self.run_test_seq_substr_sequential(es, self.dseq1, self.dseq2)
        self.run_test_total_length(es)
        self.delete_idx("foo")
        es = self.create_mem()
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(es, self.dseq1, self.dseq2)
        self.run_test_num_seqs(es)
        self.run_test_num_files_mem(es)
        self.run_test_seq_length(es)
        self.run_test_seq_substr_encoded(es, self.dseq1, self.dseq2)
        self.run_test_seq_substr_plain(es, self.dseq1, self.dseq2)
        self.run_test_seq_substr_sequential(es, self.dseq1, self.dseq2)
        self.run_test_total_length(es)

    def test_protein(self):
        self.create_es_protein("foo")
        el = EncseqLoader()
        es = el.load("foo")
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(es, self.aaseq1, self.aaseq2)
        self.run_test_num_seqs(es)
        self.run_test_num_files(es)
        self.run_test_seq_startpos_protein(es)
        self.run_test_seq_length_protein(es)
        self.run_test_file_length_protein(es)
        self.run_test_seq_substr_encoded(es, self.aaseq1, self.aaseq2)
        self.run_test_seq_substr_plain(es, self.aaseq1, self.aaseq2)
        self.run_test_seq_substr_sequential(es, self.aaseq1, self.aaseq2)
        self.run_test_total_length_protein(es)
        self.delete_idx("foo")
        es = self.create_mem_protein()
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(es, self.aaseq1, self.aaseq2)
        self.run_test_num_seqs(es)
        self.run_test_num_files_mem(es)
        self.run_test_seq_length_protein(es)
        self.run_test_seq_substr_encoded(es, self.aaseq1, self.aaseq2)
        self.run_test_seq_substr_plain(es, self.aaseq1, self.aaseq2)
        self.run_test_seq_substr_sequential(es, self.aaseq1, self.aaseq2)
        self.run_test_total_length_protein(es)

    def run_test_num_seqs(self, es):
        self.assertEquals(es.num_of_sequences(), 2)

    def run_test_num_files(self, es):
        self.assertEquals(es.num_of_files(), 1)

    def run_test_num_files_mem(self, es):
        self.assertEquals(es.num_of_files(), 1)

    def run_test_descriptions(self, es):
        self.assertRaises(GTError, es.description, 2)
        self.assertEquals(es.description(0), "seq1")
        self.assertEquals(es.description(1), "seq2")

    def run_test_total_length(self, es):
        self.assertEquals(es.total_length(), 46)

    def run_test_total_length_protein(self, es):
        self.assertEquals(es.total_length(), 37)

    def run_test_get_encoded_char(self, es, seq1, seq2):
        a = es.alphabet()
        for i, c in enumerate(seq1):
            encchar = es.get_encoded_char(i, readmode.FORWARD)
            self.assertEquals(a.decode(encchar), c)
        for i, c in enumerate(seq2[::-1]):
            encchar = es.get_encoded_char(i, readmode.REVERSE)
            self.assertEquals(a.decode(encchar), c)

    def run_test_seq_startpos(self, es):
        self.assertEquals(es.seq_startpos(0), 0)
        self.assertEquals(es.seq_startpos(1), 37)

    def run_test_seq_startpos_protein(self, es):
        self.assertEquals(es.seq_startpos(0), 0)
        self.assertEquals(es.seq_startpos(1), 31)

    def run_test_seq_length(self, es):
        self.assertEquals(es.seq_length(0), 36)
        self.assertEquals(es.seq_length(1), 9)

    def run_test_file_length(self, es):
        self.assertEquals(es.file_effective_length(0), 46)

    def run_test_seq_length_protein(self, es):
        self.assertEquals(es.seq_length(0), 30)
        self.assertEquals(es.seq_length(1), 6)

    def run_test_file_length_protein(self, es):
        self.assertEquals(es.file_effective_length(0), 37)
        
    def run_test_seq_substr_encoded(self, es, seq1, seq2):
        start = 3
        end = 13
        res = es.seq_encoded(0, start, end)
        a = es.alphabet()
        for i in range(start, end):
            self.assertEquals(a.decode(res[i-start]), seq1[i])
        start = 0
        end = 5
        res = es.seq_encoded(1, start, end)
        for i in range(start, end):
            self.assertEquals(a.decode(res[i-start]), seq2[i])

    def run_test_seq_substr_plain(self, es, seq1, seq2):
        start = 3
        end = 13
        self.assertEquals(es.seq_plain(0, start, end), seq1[start:end+1])
        start = 0
        end = 5
        self.assertEquals(es.seq_plain(1, start, end), seq2[start:end+1])

    def run_test_seq_substr_sequential(self, es, seq1, seq2):
        start = 3
        end = 13
        er = es.create_reader(readmode.FORWARD, start)
        a = es.alphabet()
        for i in range(start, end):
            self.assertEquals(a.decode(er.next_encoded_char()), seq1[i])
        start = es.seq_startpos(1)
        end = start + 5
        er = es.create_reader(readmode.FORWARD, start)
        for i in range(start-start, end-start):
            self.assertEquals(a.decode(er.next_encoded_char()), seq2[i])

if __name__ == "__main__":
    unittest.main()
