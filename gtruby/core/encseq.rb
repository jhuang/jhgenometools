#
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
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

require 'gtdlload'
require 'core/alphabet'
require 'core/error'
require 'core/readmode'
require 'core/str'
require 'dl/struct'
require 'dl'

module GT
  extend DL::Importable
  gtdlload "libgenometools"
  extern "GtEncseqEncoder* gt_encseq_encoder_new()"
  extern "int gt_encseq_encoder_use_representation(GtEncseqEncoder*, GtStr*,
                                                   GtError*)"
  extern "int gt_encseq_encoder_use_symbolmap_file(GtEncseqEncoder*, GtStr*,
                                                   GtError*)"
  extern "void gt_encseq_encoder_enable_description_support(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_disable_description_support(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_enable_multiseq_support(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_disable_multiseq_support(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_create_tis_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_do_not_create_tis_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_create_des_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_do_not_create_des_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_create_ssp_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_do_not_create_ssp_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_create_sds_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_do_not_create_sds_tab(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_set_input_dna(GtEncseqEncoder*)"
  extern "void gt_encseq_encoder_set_input_protein(GtEncseqEncoder*)"
  extern "int  gt_encseq_encoder_encode(GtEncseqEncoder*, GtStrArray*,
                                        GtStr*, GtError*)"
  extern "void gt_encseq_encoder_delete(GtEncseqEncoder*)"

  extern "GtEncseqLoader* gt_encseq_loader_new()"
  extern "void gt_encseq_loader_require_description_support(GtEncseqLoader*)"
  extern "void gt_encseq_loader_drop_description_support(GtEncseqLoader*)"
  extern "void gt_encseq_loader_require_multiseq_support(GtEncseqLoader*)"
  extern "void gt_encseq_loader_drop_multiseq_support(GtEncseqLoader*)"
  extern "void gt_encseq_loader_require_tis_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_do_not_require_tis_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_require_des_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_do_not_require_des_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_require_ssp_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_do_not_require_ssp_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_require_sds_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_do_not_require_sds_tab(GtEncseqLoader*)"
  extern "void gt_encseq_loader_enable_range_iterator(GtEncseqLoader*)"
  extern "void gt_encseq_loader_disable_range_iterator(GtEncseqLoader*)"
  extern "GtEncseq* gt_encseq_loader_load(GtEncseqLoader*, GtStr*, GtError*)"
  extern "void gt_encseq_loader_delete(GtEncseqLoader*)"

  extern "GtEncseqBuilder* gt_encseq_builder_new(GtAlphabet*)"
  extern "void gt_encseq_builder_enable_range_iterator(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_disable_range_iterator(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_enable_description_support(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_disable_description_support(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_enable_multiseq_support(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_disable_multiseq_support(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_create_des_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_do_not_create_des_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_create_ssp_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_do_not_create_ssp_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_create_sds_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_do_not_create_sds_tab(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_add_cstr(GtEncseqBuilder*,const char*,
                                          unsigned long, const char*)"
  extern "void gt_encseq_builder_add_encoded(GtEncseqBuilder*, const GtUchar*,
                                             unsigned long, const char*)"
  extern "GtEncseq* gt_encseq_builder_build(GtEncseqBuilder*, GtError*)"
  extern "void gt_encseq_builder_reset(GtEncseqBuilder*)"
  extern "void gt_encseq_builder_delete(GtEncseqBuilder*)"

  typealias "GtUchar", "unsigned char"
  typealias "GtUlong", "unsigned long"
  typealias "GtReadmode", "int"
  extern "unsigned long gt_encseq_total_length(const GtEncseq*)"
  extern "unsigned long gt_encseq_num_of_sequences(const GtEncseq*)"
  extern "GtUchar gt_encseq_get_encoded_char(const GtEncseq*, unsigned long,
                                             GtReadmode)"
  extern "void gt_encseq_extract_substring(const GtEncseq*, GtUchar*,
                                           unsigned long, unsigned long)"
  extern "unsigned long gt_encseq_seqlength(const GtEncseq*, unsigned long)"
  extern "unsigned long gt_encseq_seqstartpos(const GtEncseq *, unsigned long)"
  extern "void* gt_encseq_description(const GtEncseq*, unsigned long*,
                                      unsigned long)"
  extern "GtEncseqReader* gt_encseq_create_reader_with_readmode(const GtEncseq*,
                                                                GtReadmode,
                                                                unsigned long)"
  extern "GtAlphabet* gt_encseq_alphabet(const GtEncseq*)"
  extern "const GtStrArray* gt_encseq_filenames(const GtEncseq*)"
  extern "void gt_encseq_delete(GtEncseq*)"

  extern "GtUchar gt_encseq_reader_next_encoded_char(GtEncseqReader*)"
  extern "void gt_encseq_reader_delete(GtEncseqReader*)"

  class EncseqEncoder
    def initialize(eenc = nil, own = true)
      if eenc.nil? then
        @eenc = GT.gt_encseq_encoder_new()
      else
        @eenc = eenc 
      end
      if own then
        @eenc.free = GT::symbol("gt_encseq_encoder_delete", "0P")
      end
    end

    def use_representation(sat)
      strsat = Str.new(sat)
      err = Error.new
      if (GT.gt_encseq_encoder_use_representation(@eenc, strsat, err) != 0)
        GT.gterror(err)
      end
    end

    def use_symbolmap_file(smap)
      strsmap = Str.new(smap)
      err = Error.new
      if (GT.gt_encseq_encoder_use_symbolmap_file(@eenc, strsmap, err) != 0)
        GT.gterror(err)
      end
    end

    def enable_description_support
      GT.gt_encseq_encoder_enable_description_support(@eenc)
    end

    def disable_description_support
      GT.gt_encseq_encoder_disable_description_support(@eenc)
    end

    def enable_multiseq_support
      GT.gt_encseq_encoder_enable_multiseq_support(@eenc)
    end

    def disable_multiseq_support
      GT.gt_encseq_encoder_disable_multiseq_support(@eenc)
    end

    def create_tis_tab
      GT.gt_encseq_encoder_create_tis_tab(@eenc)
    end

    def do_not_create_tis_tab
      GT.gt_encseq_encoder_do_not_create_tis_tab(@eenc)
    end

    def create_des_tab
      GT.gt_encseq_encoder_create_des_tab(@eenc)
    end

    def do_not_create_des_tab
      GT.gt_encseq_encoder_do_not_create_des_tab(@eenc)
    end

    def create_ssp_tab
      GT.gt_encseq_encoder_create_ssp_tab(@eenc)
    end

    def do_not_create_ssp_tab
      GT.gt_encseq_encoder_do_not_create_ssp_tab(@eenc)
    end

    def create_sds_tab
      GT.gt_encseq_encoder_create_sds_tab(@eenc)
    end

    def do_not_create_sds_tab
      GT.gt_encseq_encoder_do_not_create_sds_tab(@eenc)
    end

    def set_input_dna
      GT.gt_encseq_encoder_set_input_dna(@eenc)
    end

    def set_input_protein
      GT.gt_encseq_encoder_set_input_protein(@eenc)
    end

    def encode(files, indexname)
      files = files.to_a
      if files.length == 0 then
        GT.gterror("list of input sequence files must be non-empty")
      end
      sa = StrArray.new
      files.each do |f|
        fn = f.to_s
        if (!File.exists?(fn))
          GT.gterror("file not found: #{fn}")
        end
        sa.add(fn)
      end
      idxname = Str.new(indexname.to_s)
      err = Error.new
      rval = GT.gt_encseq_encoder_encode(@eenc, sa.to_ptr, idxname.to_ptr, \
                                         err.to_ptr)
      if rval != 0 then
        GT.gterror(err)
      end
    end
    
    def to_ptr
      @eenc
    end
  end

  class EncseqLoader
    def initialize(eldr = nil, own = true)
      if eldr.nil? then
        @eldr = GT.gt_encseq_loader_new()
      else
        @eldr = eldr 
      end
      if own then
        @eldr.free = GT::symbol("gt_encseq_loader_delete", "0P")
      end
      @destab = true
      @ssptab = true
      @sdstab = true
    end

    def require_description_support
      @destab = true
      @sdstab = true
      GT.gt_encseq_loader_require_description_support(@eldr)
    end

    def drop_description_support
      @destab = false
      @sdstab = false
      GT.gt_encseq_loader_drop_description_support(@eldr)
    end

    def require_multiseq_support
      @ssptab = true
      GT.gt_encseq_loader_require_multiseq_support(@eldr)
    end

    def drop_multiseq_support
      @ssptab = false
      GT.gt_encseq_loader_drop_multiseq_support(@eldr)
    end

    def require_tis_tab
      GT.gt_encseq_loader_require_tis_tab(@eldr)
    end

    def do_not_require_tis_tab
      GT.gt_encseq_loader_do_not_require_tis_tab(@eldr)
    end

    def require_des_tab
      @destab = true
      GT.gt_encseq_loader_require_des_tab(@eldr)
    end

    def do_not_require_des_tab
      @destab = false
      GT.gt_encseq_loader_do_not_require_des_tab(@eldr)
    end

    def require_ssp_tab
      @ssptab = true
      GT.gt_encseq_loader_require_ssp_tab(@eldr)
    end

    def do_not_require_ssp_tab
       @ssptab = false
      GT.gt_encseq_loader_do_not_require_ssp_tab(@eldr)
    end

    def require_sds_tab
      @sdstab = true
      GT.gt_encseq_loader_require_sds_tab(@eldr)
    end

    def do_not_require_sds_tab
      @sdstab = false
      GT.gt_encseq_loader_do_not_require_sds_tab(@eldr)
    end

    def enable_range_iterator
      GT.gt_encseq_loader_enable_range_iterator(@eldr)
    end
    
    def disable_range_iterator
      GT.gt_encseq_loader_disable_range_iterator(@eldr)
    end

    def load(indexname)
      if !File.exists?("#{indexname}.esq") then
        GT.gterror("file not found: #{indexname}.esq")
      end
      if @destab and !File.exists?("#{indexname}.des") then
        GT.gterror("file not found: #{indexname}.des")
      end
      if @ssptab and !File.exists?("#{indexname}.ssp") then
        GT.gterror("file not found: #{indexname}.ssp")
      end
      if @sdstab and !File.exists?("#{indexname}.sds") then
        GT.gterror("file not found: #{indexname}.sds")
      end
      idxname = Str.new(indexname.to_s)
      err = Error.new
      rval = GT.gt_encseq_loader_load(@eldr, idxname.to_ptr, err.to_ptr)
      if rval == GT::NULL then
        GT.gterror(err)
      end
      Encseq.new(rval, true)
    end
    
    def to_ptr
      @eldr
    end
  end

  class EncseqBuilder

    def self.create(alphabet)
      if !alphabet.is_a?(Alphabet) then
        GT.gterror("argument must be an Alphabet")
      end
      EncseqBuilder.new(GT.gt_encseq_builder_new(alphabet.to_ptr))
    end
  
    def initialize(ebld, own = true)
      @ebld = ebld
      if own then
        @ebld.free = GT::symbol("gt_encseq_builder_delete", "0P")
      end
    end

    def enable_description_support
      GT.gt_encseq_builder_enable_description_support(@ebld)
    end

    def disable_description_support
      GT.gt_encseq_builder_disable_description_support(@ebld)
    end

    def enable_multiseq_support
      GT.gt_encseq_builder_enable_multiseq_support(@ebld)
    end

    def disable_multiseq_support
      GT.gt_encseq_builder_disable_multiseq_support(@ebld)
    end

    def create_tis_tab
      GT.gt_encseq_builder_create_tis_tab(@ebld)
    end

    def do_not_create_tis_tab
      GT.gt_encseq_builder_do_not_create_tis_tab(@ebld)
    end

    def create_des_tab
      GT.gt_encseq_builder_create_des_tab(@ebld)
    end

    def do_not_create_des_tab
      GT.gt_encseq_builder_do_not_create_des_tab(@ebld)
    end

    def create_ssp_tab
      GT.gt_encseq_builder_create_ssp_tab(@ebld)
    end

    def do_not_create_ssp_tab
      GT.gt_encseq_builder_do_not_create_ssp_tab(@ebld)
    end

    def create_sds_tab
      GT.gt_encseq_builder_create_sds_tab(@ebld)
    end

    def do_not_create_sds_tab
      GT.gt_encseq_builder_do_not_create_sds_tab(@ebld)
    end

    def enable_range_iterator
      GT.gt_encseq_builder_enable_range_iterator(@eldr)
    end
    
    def disable_range_iterator
      GT.gt_encseq_builder_disable_range_iterator(@eldr)
    end

    def add_string(string, desc = '')
      str = string.to_s
      GT.gt_encseq_builder_add_cstr(@ebld, str, str.length, desc)
    end

    def build
      err = Error.new
      rval = GT.gt_encseq_builder_build(@ebld, err.to_ptr)
      if rval == GT::NULL then
        GT.gterror(err)
      end
      Encseq.new(rval, true)
    end
    
    def to_ptr
      @ebld
    end
  end

  UlongParam = struct [
    "GtUlong val"
  ]

  SeqInfo = struct [
    "GtUlong startpos",
    "GtUlong length"
  ]

  class Encseq
    def initialize(encseq, own = true)
      @encseq = encseq
      if own then
        @encseq.free = GT::symbol("gt_encseq_delete", "0P")
      end
    end

    def to_ptr
      @encseq
    end

    def num_of_sequences
      return GT.gt_encseq_num_of_sequences(@encseq)
    end

    def description(num)
      num = num.to_i
      if num >= self.num_of_sequences then
        GT.gterror("invalid sequence number #{num}")
      end
      len = UlongParam.malloc
      charptr = GT.gt_encseq_description(@encseq, len, num)
      if len.val > 0 and charptr != GT::NULL then
        charptr.to_s(len.val)
      else
        ""
      end
    end

    def get_encoded_char(pos, readmode = GT::READMODE_FORWARD)
      if readmode < 0 or readmode > 3 then
          GT.gterror("invalid readmode!")
      end
      return GT.gt_encseq_get_encoded_char(@encseq, pos, readmode)
    end
   
    def total_length
      GT.gt_encseq_total_length(@encseq)
    end

    def alphabet
      aptr = GT.gt_encseq_alphabet(@encseq)
      if aptr == GT::NULL then
        nil
      else
        Alphabet.new(aptr, false)
      end
    end

    def seq_length(num)
      if num >= self.num_of_sequences then
        GT.gterror("invalid sequence number #{num}")
      end
      GT.gt_encseq_seqlength(@encseq, num)
    end

    def seq_startpos(num)
      if num >= self.num_of_sequences then
        GT.gterror("invalid sequence number #{num}")
      end
      GT.gt_encseq_seqstartpos(@encseq, num)
    end

    def create_reader(readmode, startpos)
      if readmode < 0 or readmode > 3 then
          GT.gterror("invalid readmode!")
      end
      if startpos < 0 or startpos >= total_length then
        gterror("invalid start position: #{startpos} " + \
                "(allowed: 0-#{total_length-1})")
      end
      er_ptr = GT.gt_encseq_create_reader_with_readmode(@encseq, readmode, \
                                                        startpos)
      return EncseqReader.new(er_ptr)
    end

    def seq_encoded(num, start, stop)
      if num >= self.num_of_sequences then
        GT.gterror("invalid sequence number #{num}")
      end
      seqlength = self.seq_length(num)
      seqstartpos = self.seq_startpos(num)
      if start < 0 or stop >= seqlength then
        gterror("invalid coordinates: #{start}-#{stop} " + \
                "(allowed: #{0}-#{seqlength-1})")
      end
      buf = DL.malloc(DL::sizeof('C') * (stop-start+1))
      GT.gt_encseq_extract_substring(@encseq, buf, seqstartpos+start, \
                                     seqstartpos+stop)
      buf.to_a('C')
    end

    def seq_plain(num, start, stop)
      self.alphabet.decode_seq(self.seq_encoded(num, start, stop))
    end
  end

   class EncseqReader
    def initialize(er, own = true)
      @er = er
      if own then
        @er.free = GT::symbol("gt_encseq_reader_delete", "0P")
      end
    end

    def next_encoded_char
      return GT.gt_encseq_reader_next_encoded_char(@er)
    end

    def to_ptr
      @er
    end
  end
end