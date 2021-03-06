def outoptionsnobck
  return "-tis -suf -des -sds -ssp -lcp -bwt"
end

def outoptions
  return outoptionsnobck + " -bck"
end

def trials()
  return "-scantrials 10 -multicharcmptrials 1000"
end

def checksfx(parts,withsmap,cmp,doubling,filelist,alldirs=true)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  if alldirs 
    dirlist=["fwd","rev","cpl ","rcl "]
  else
    dirlist=["fwd ","rev"]
  end
  dirlist.each do |dirarg|
    extra=""
    if cmp
      extra=extra + " -cmpcharbychar"
      if dirarg == "fwd" and doubling
        extra=extra + " -maxdepth"
      end
    end
    run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl " +
             "-algbds 10 31 80 #{extra} #{outoptions} " +
             "-indexname esa -dir " + dirarg + " -db " + filearg
    if dirarg == "cpl" or dirarg = "rcl"
      run_test "#{$bin}gt dev sfxmap #{trials()} #{outoptions} -v " +
               "-esa esa",
               :maxtime => 600
    else
      if dirarg == "fwd"
        dirarg_rev = "rev"
      else
        dirarg_rev = "fwd"
      end
      run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl " +
               "-algbds 10 31 80 #{extra} #{outoptions} " +
               "-indexname esa-#{dirarg_rev} -dir " + dirarg_rev + 
               " -db " + filearg
      run_test "#{$bin}gt packedindex mkindex -indexname pck -dir " + dirarg +
               " -db " + filearg
      run_test "#{$bin}gt dev sfxmap #{trials()} #{outoptions} -v " +
               "-esa esa -pck pck -cmpsuf",
               :maxtime => 600
      run_test "#{$bin}gt dev sfxmap #{trials()} #{outoptions} -v " +
               "-esa esa-#{dirarg_rev} -pck pck -cmplcp",
               :maxtime => 600
    end
  end
end

def checkdc(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -v -pl -dc 64 -suf -ssp -tis " +
           "-indexname sfx -db " + filearg
  run_test "#{$bin}gt dev sfxmap #{trials()} -suf -tis -ssp -v -esa sfx",
           :maxtime => 600
  run_test "#{$bin}gt suffixerator -v -pl -parts 3 -dc 64 -suf -tis " +
           "-indexname sfx3 -db " + filearg
  run "diff sfx3.suf sfx.suf"
end

def flattenfilelist(filelist)
  s=""
  filelist.each do |f|
    s += "#{$testdata}#{f} "
  end
  return s
end

def checkbwt(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -pl #{outoptions} -indexname sfx -db " +
           flattenfilelist(filelist)
end

allfiles = []
all_fastafiles = ["Atinsert.fna",
                  "Duplicate.fna",
                  "Random-Small.fna",
                  "Random.fna",
                  "Copysorttest.fna",
                  "Random159.fna",
                  "Random160.fna",
                  "RandomN.fna",
                  "TTT-small.fna",
                  "trna_glutamine.fna",
                  "atC99826.fna"]

allfiles += all_fastafiles
allfiles += (all_genbankfiles = all_fastafiles.collect{ |f|
                                                        f.gsub(".fna",".gbk")
                                                      })
allfiles += (all_emblfiles = all_fastafiles.collect{ |f|
                                                     f.gsub(".fna",".embl")
                                                   })

allmultifiles = []
all_multifastafiles = ["Atinsert.fna",
                       "Duplicate.fna",
                       "Random159.fna",
                       "Random160.fna"]

allmultifiles += all_multifastafiles
allmultifiles += (all_multigenbankfiles = all_multifastafiles.collect{ |f|
                                                         f.gsub(".fna",".gbk")
                                                                     })
allmultifiles += (all_multiemblfiles = all_multifastafiles.collect{ |f|
                                                         f.gsub(".fna",".embl")
                                                                  })

all_fastqfiles = ["fastq_long.fastq",
                  "test10_multiline.fastq",
                  "test1.fastq",
                  "test5_tricky.fastq"]

allmultifiles += all_fastqfiles
allfiles += all_fastqfiles

alldir = ["fwd","cpl","rev","rcl"]

{"FASTA" => all_fastafiles,
 "EMBL" => all_emblfiles,
 "GenBank" => all_genbankfiles,
 "FastQ" => all_fastqfiles}.each do |k,filelist|
    Name "gt suffixerator (#{k})"
    Keywords "gt_suffixerator tis"
    Test do
    run_test "#{$bin}gt suffixerator -tis -ssp -indexname sfx -db " +
             flattenfilelist(filelist)
    run_test "#{$bin}gt dev sfxmap -ssp -tis -esa sfx"
    filelist.each do |filename|
      run_test "#{$bin}gt suffixerator -tis -ssp -indexname sfx -db " +
               "#{$testdata}#{filename}"
      run_test "#{$bin}gt dev sfxmap -ssp -tis -esa sfx"
    end
  end
end

allfiles.each do |filename|
  Name "gt suffixerator #{filename}"
  Keywords "gt_suffixerator tis"
  Test do
    ["direct", "bit", "uchar", "ushort", "uint32"].each do |sat|
      run_test "#{$bin}gt suffixerator -tis -indexname sfx -sat #{sat} " +
               "-db #{$testdata}#{filename}"
    end
  end
end

Name "gt suffixerator file of reads of equal length"
Keywords "gt_suffixerator reads"
Test do
  run_test "#{$bin}/gt suffixerator -des -tis -ssp -dna " +
           "-db #{$testdata}U89959_genomic.fas -indexname u8idx"
  run_test "#{$bin}/gt simreads -coverage 4 -len 100 -force -o u8.reads u8idx"
  run_test "#{$bin}/gt suffixerator -v #{outoptionsnobck} -dna -db u8.reads"
  run "grep -q '# init character encoding (eqlen ' #{$last_stdout}"
  run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa u8.reads", \
           :maxtime => 200
  run_test "#{$bin}/gt suffixerator -v #{outoptionsnobck} -dir rev -dna -db u8.reads"
  run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa u8.reads", \
           :maxtime => 200
end

all_fastafiles.each do |filename|
  Name "gt suffixerator -dc 64 -parts 1+3 #{filename}"
  Keywords "gt_suffixerator dc"
  Test do
    checkdc([filename])
  end
end

Name "gt suffixerator -dc 64 -parts 1+3 all-fastafiles"
Keywords "gt_suffixerator dc"
Test do
  checkdc(all_fastafiles)
end

alldir.each do |dir|
  {"FASTA" => all_fastafiles,
   "EMBL" => all_emblfiles,
   "GenBank" => all_genbankfiles,
   "FastQ" => all_fastqfiles}.each do |k,filelist|
    Name "gt suffixerator #{dir} (#{k})"
    Keywords "gt_suffixerator"
    Test do
       run_test "#{$bin}gt suffixerator -dir #{dir} -tis -suf -bwt -lcp " +
                "-indexname sfx -pl -db " +
                flattenfilelist(filelist)
       run_test "#{$bin}gt suffixerator -storespecialcodes -dir #{dir} -tis " +
                "-suf -lcp -indexname sfx -pl -db " +
                flattenfilelist(filelist)
       run_test "#{$bin}gt suffixerator -tis -bwt -lcp -pl -ii sfx"
    end
  end
end

faillist = ["-indexname sfx -db /nothing",
            "-indexname /nothing/sfx -db #{$testdata}TTT-small.fna",
            "-smap /nothing -db #{$testdata}TTT-small.fna",
            "-dna -db #{$testdata}sw100K1.fsa",
            "-protein -dir cpl -db #{$testdata}sw100K1.fsa",
            "-dna -db #{$testdata}Random.fna RandomN.fna",
            "-dna -suf -pl 10 -db #{$testdata}Random.fna",
            "-dna -tis -sat plain -db #{$testdata}TTT-small.fna"]

faillist.each do |failcommand|
  Name "gt suffixerator failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis " + failcommand,:retval => 1
  end
end

allmultifiles.each do |filename|
  Name "gt suffixerator sfxmap-failure #{filename}"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -dna -indexname localidx " +
             "-db #{$testdata}#{filename}"
    run_test "#{$bin}gt suffixerator -suf -lcp -pl -dir rev -ii localidx"
    run_test "#{$bin}gt dev sfxmap -tis -des -esa localidx",
             :retval => 1
    # In short read files with equal read lengths, ssptabs need not
    # be built explicitly.
    # Thus these tests only fail for non-equal length files.
    if !all_fastqfiles.include?(filename) then
      run_test "#{$bin}gt dev sfxmap -tis -ssp -esa localidx",
               :retval => 1
      run_test "#{$bin}gt dev sfxmap -ssp -esa localidx",
               :retval => 1
    end
    run_test "#{$bin}gt dev sfxmap -des -esa localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -tis -bck -esa localidx",
             :retval => 1
  end
end

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(all_fastafiles)
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    extra=""
    if withsmap == 1
      extra="-protein"
      extraname="protein"
    elsif withsmap == 2
      extra="-smap TransProt11"
      extraname="TransProt11"
    end
    if parts == 1
     doubling=true
    else
     doubling=false
    end
    Name "gt suffixerator+map protein filelist #{extraname} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,extra,true,doubling,
               ["sw100K1.fsa","sw100K2.fsa"],false)
    end
  end
end

0.upto(2) do |cmpval|
  1.upto(2) do |parts|
    [0,2].each do |withsmap|
      extra=""
      if withsmap == 1
	extra="-dna"
	extraname=" dna"
      elsif withsmap == 2
	extra="-smap TransDNA"
	extraname=" trans"
      end
      doublingname=""
      if cmpval == 0
	cmp=false
	doubling=false
      elsif cmpval == 1
	cmp=true
	doubling=false
      else
	cmp=true
	if parts == 1
	  doubling=true
	  doublingname=" doubling "
	else
	  doubling=false
	end
      end
      all_fastafiles.each do |filename|
	Name "gt suffixerator+map #{filename}#{extraname} #{parts} parts " +
             "#{doubling}"
	Keywords "gt_suffixerator"
	Test do
	  checksfx(parts,extra,cmp,doubling,[filename])
	end
      end
      filelist=["RandomN.fna","Random.fna","Atinsert.fna"]
      Name "gt suffixerator+map dna filelist#{extraname} " +
	     "#{parts} parts #{doubling}"
      Keywords "gt_suffixerator"
      Test do
	checksfx(parts,extra,cmp,doubling,filelist)
      end
    end
  end
end

def checkmapped(args)
  Name "gt suffixerator checkmapped"
  Keywords "gt_suffixerator gttestdata"
  Test do
    run_test "#{$bin}gt suffixerator #{outoptions} -algbds 3 34 90 " +
             "-indexname sfxidx #{args}",
             :maxtime => 1200
    run_test "#{$bin}gt dev sfxmap #{outoptions} #{trials()} -v -esa sfxidx",
             :maxtime => 2400
    run_test "#{$bin}gt dev sfxmap #{outoptionsnobck} -stream -v -esa sfxidx",
             :maxtime => 2400
  end
end

def grumbach()
  return "#{$gttestdata}DNA-mix/Grumbach.fna/"
end

if $gttestdata then
  checkmapped("-db " +
              "#{$gttestdata}Iowa/at100K1 " +
              "#{grumbach()}Wildcards.fna " +
              "#{grumbach()}chntxx.fna " +
              "#{grumbach()}hs5hcmvcg.fna " +
              "#{grumbach()}humdystrop.fna " +
              "#{grumbach()}humghcsa.fna " +
              "#{grumbach()}humhdabcd.fna " +
              "#{grumbach()}humhprtb.fna " +
              "#{grumbach()}mipacga.fna " +
              "#{grumbach()}mpocpcg.fna " +
              "#{grumbach()}ychrIII.fna " +
              "-parts 3 -pl")

  checkmapped("-parts 1 -pl -db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB -parts 3 -pl")

  checkmapped("-parts 2 -pl -smap TransDNA -db  #{$gttestdata}Iowa/at100K1")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K -parts 1 -pl -smap " +
              "TransProt11")
end

SATS = ["direct", "bytecompress", "eqlen", "bit", "uchar", "ushort", "uint32"]

EQLENDNAFILE = {:filename => "#{$testdata}test1.fasta",
                :desc => "equal length DNA",
                :msgs => {
                  "bytecompress" => "cannot use bytecompress on DNA sequences"}}
DNAFILE   = {:filename => "#{$testdata}at1MB",
             :desc => "non-equal length DNA",
             :msgs => {
                "bytecompress" => "cannot use bytecompress on DNA sequences",
                "eqlen" => "all sequences are of equal length and no " + \
                "sequence contains"}}
EQLENAAFILE = {:filename => "#{$testdata}trembl-eqlen.faa",
                :desc => "equal length AA",
                :msgs => {
                  "eqlen" => "as the sequence is not DNA",
                  "bit" => "as the sequence is not DNA",
                  "uchar" => "as the sequence is not DNA",
                  "ushort" => "as the sequence is not DNA",
                  "uint32" => "as the sequence is not DNA"}}
AAFILE    = {:filename => "#{$testdata}trembl.faa",
                :desc => "non-equal length AA",
                :msgs => {
                  "eqlen" => "as the sequence is not DNA",
                  "bit" => "as the sequence is not DNA",
                  "uchar" => "as the sequence is not DNA",
                  "ushort" => "as the sequence is not DNA",
                  "uint32" => "as the sequence is not DNA"}}

SATTESTFILES = [EQLENDNAFILE, DNAFILE, EQLENAAFILE, AAFILE]

SATTESTFILES.each do |file|
  SATS.each do |sat|
    Name "gt suffixerator sat #{sat} -> #{file[:desc]}"
    Keywords "gt_suffixerator sats"
    Test do
      if !file[:msgs][sat].nil? then
        retval = 1
      else
        retval = 0
      end
      run_test "#{$bin}/gt suffixerator -sat #{sat} -v -suf -lcp -des -sds " + \
               "-ssp -tis -db #{file[:filename]} -indexname myidx", \
               :retval => retval
      if !file[:msgs][sat].nil? then
        grep($last_stderr, /#{file[:msgs][sat]}/)
      end
      run_test "#{$bin}/gt dev sfxmap -suf -lcp -des -sds -ssp -esa myidx", \
               :retval => retval
    end
  end
end  
