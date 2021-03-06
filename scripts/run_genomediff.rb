#!/usr/bin/ruby
#
# Copyright (c) 2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
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
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'fileutils'
require 'genomediff'
require 'optparse'
require 'ostruct'

GT_DIR = File.join(File.dirname(__FILE__), "..")
ENV['GTDIR'] = GT_DIR

options = OpenStruct.new
options.esa = true
options.redN = false
options.parts = 1
options.name = "esa"

opt= OptionParser.new do |opt|
  opt.banner = "USAGE: #$0 [options] <files>\n
    The Program calls the tool Genomediff after preparing the data.\n
    Files for which the N-stretches are to be removed, have to be FASTA.
    the indices will be named esa or pck (packed), this can be changed by
    the --name flag"
  opt.separator ""
  opt.separator "Options"
  opt.on("--fm", "Use FM-Index rather then ESA") do
    options.esa = false
    if options.name == "esa"
      options.name = "pck"
    end
  end
  opt.on("--reduceN", "reduce stretches of Wildcards in the sequences to one",
         "single N, only works with FASTA!") do |val|
    options.redN = val
  end
  opt.on("-p", "--parts PARTS", Integer, "number of PARTS to build index,",
         "reduces peak memory during index construction") do |val|
    options.parts = val
  end
  opt.on("--name NAME", String, "the baseNAME of the index files") do |val|
    options.name = val
  end 
  opt.on_tail("-h", "--help", "prints this help and exits") do
    puts opt
    exit 0
  end
end

opt.parse!(ARGV)

files = ""
ARGV.each do |file|
  unless File.exist?(file)
    puts "non existing file " + file
  else
    if options.redN
      file = Genomediff.reduceN(file)
    end
    files += Genomediff.reverse_and_concat(file, false) + " "
  end
end


if options.esa
  puts "***ESA***"
  puts Genomediff.esa_index(files, options.parts, options.name)
  puts Genomediff.esa_genomediff(options.name,"")
else
  puts "***FM-INDEX***"
  puts Genomediff.pck_index(files, 8, options.parts, options.name)
  puts Genomediff.pck_genomediff(options.name,"")
end
