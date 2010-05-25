#allfiles = ["Atinsert.fna",
            #"Duplicate.fna",
            #"Random-Small.fna",
            #"Random.fna",
            #"Random159.fna",
            #"Random160.fna",
            #"RandomN.fna",
            #"TTT-small.fna",
            #"trna_glutamine.fna"]

repfindtestfiles=["Duplicate.fna",
                  "Wildcards.fna",
                  "hs5hcmvcg.fna",
                  "humhbb.fna",
                  "mpomtcg.fna",
                  "at1MB",
                  "ychrIII.fna"]
#repfindtestfiles=["Duplicate.fna",
                  #"at1MB"]
def addfilepath(filename)
  if filename == 'Duplicate.fna' or filename == 'at1MB'
    return "#{$testdata}#{filename}"
  else
    return "#{$gttestdata}DNA-mix/Grumbach.fna/#{filename}"
  end
end

def determinemumreferenceminlength(reffile)
  if reffile == 'at1MB'
    return 22
  else
    return 14
  end
end

def determinemaxmatchminlength(reffile)
  if reffile == 'at1MB'
    return 22
  else
    return 14
  end
end


#def checktagerator(queryfile,ms)
  #run "#{$bin}gt shredder -minlength 12 -maxlength 15 #{queryfile} | " +
      #"#{$bin}gt seqfilter -minlength 12 - | " +
      #"sed -e \'s/^>.*/>/\' > patternfile"
  #if File.size("patternfile") > 0
    #run_test("#{$bin}gt tagerator -rw -cmp -e 0 -esa sfx -q patternfile",
             #:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -e 1 -esa sfx -q patternfile " +
             #"-withwildcards",:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -e 2 -esa sfx -q patternfile " +
             #"-withwildcards",:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -esa sfx -q patternfile " +
             #" -maxocc 10",
             #:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -e 0 -pck pck -q patternfile",
             #:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -e 1 -pck pck -q patternfile",
             #:maxtime => 100)
    #run_test("#{$bin}gt tagerator -rw -cmp -e 2 -pck pck -q patternfile",
             #:maxtime => 200)
    #run_test("#{$bin}gt tagerator -rw -cmp -pck pck -q patternfile " +
             #"-maxocc 10",
             #:maxtime => 100)
  #end
#end







#def determineminlength(reffile)
  #if reffile == 'Duplicate.fna'
    #return 6
  #elsif reffile == 'Wildcards.fna'
    #return 6
  #elsif reffile == 'hs5hcmvcg.fna'
    #return 10
  #elsif reffile == 'humhbb.fna'
    #return 10
  #elsif reffile == 'mpomtcg.fna'
    #return 10
  #elsif reffile == 'at1MB'
    #return 18
  #elsif reffile == 'ychrIII.fna'  
    #return 10
  #else
    #return 14
  #end
#end




## 自身检查，在我这里没有
#def checkrepfind(reffile)
  #reffilepath=addfilepath(reffile)
  #run_test("#{$bin}gt suffixerator -algbds 3 40 120 -db " +
           #"#{reffilepath} -indexname sfxidx -dna -suf -tis -lcp -ssp -pl",
           #:maxtime => 320)
  #minlength = determineminlength(reffile)
  #run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx", :maxtime => 320)
  #resultfile="#{$gttestdata}repfind-result/#{reffile}.result"
  #run "cmp -s #{$last_stdout} #{resultfile}"
  #run_test("#{$bin}gt repfind -l #{minlength} -r -ii sfxidx", :maxtime => 320)
  #resultfile="#{$gttestdata}repfind-result/#{reffile}-r.result"
  #run "cmp -s #{$last_stdout} #{resultfile}"
#end



def checkmaxmat4withmummer(reffile,queryfile,matchmode)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)

  pckname=reffile + "-pck"                                                 
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{reffilepath} -indexname #{pckname} -sprank -dna -ssp -des -sds -pl"   #-sprank -dna -tis -ssp -des -sds -pl
  run_test("#{$bin}gt dev maxmat4 -#{matchmode} -b -l #{determinemaxmatchminlength(reffile)} -L -s -c #{pckname} #{queryfilepath}", :maxtime => 320)
  #run "grep -v '^>' #{$last_stdout} | sort"  
  if (matchmode=="mumreference") 
      run "sed -e '/^>/d' #{$last_stdout}"
  elsif (matchmode=="maxmatch") 
      run "sed -e '/^>/d' #{$last_stdout} | sort"
  end             
  run "diff -w #{$last_stdout} #{$gttestdata}maxmat4-result/#{reffile}-#{queryfile}_#{matchmode}.result"
end


def generatemummerresults(reffile, queryfile, matchmode)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)
  run "~/Desktop/maxmat3/mm3src/maxmat3.x -#{matchmode} -b -l #{determinemaxmatchminlength(reffile)} -L -s -c #{reffilepath} #{queryfilepath}"
  #run "grep -v '^>' #{$last_stdout} | sort" 
  if (matchmode=="mumreference") 
      run "sed -e '/^>/d' #{$last_stdout}"
  elsif (matchmode=="maxmatch") 
      run "sed -e '/^>/d' #{$last_stdout} | sort"  
  end  
  run "cp #{$last_stdout} #{$gttestdata}maxmat4-result/#{reffile}-#{queryfile}_#{matchmode}.result"
end



if $gttestdata then
  repfindtestfiles.each do |reffile|  
    repfindtestfiles.each do |queryfile|
      if reffile != queryfile
        Name "generate benchmark-file #{reffile}-#{queryfile}_mumreference.result with MUMmer"
				Keywords "generate_mummer_mumreference_benchmark"
				Test do
					generatemummerresults(reffile,queryfile,"mumreference")
				end
        Name "generate benchmark-file #{reffile}-#{queryfile}_maxmatch.result with MUMmer"
				Keywords "generate_mummer_maxmatch_benchmark"
				Test do
					generatemummerresults(reffile,queryfile,"maxmatch")
				end

        Name "check maxmat4 -mumreference #{reffile} versus #{queryfile} comparing with MUMmer"
        Keywords "check_maxmat4_mumreference_with_mummer"
        Test do
          checkmaxmat4withmummer(reffile,queryfile,"mumreference")
        end
        Name "check maxmat4 -maxmatch #{reffile} versus #{queryfile} comparing with MUMmer"
        Keywords "check_maxmat4_maxmatch_with_mummer"
        Test do
          checkmaxmat4withmummer(reffile,queryfile,"maxmatch")
        end       
      end
    end
  end
end

#Name "gt maxmat4 small"
#Keywords "gt_maxmat4_repfind"   # 这里对应jhtest.sh
#Test do
  #puts $bin               # ./bin/
  #puts $testdata          # ./testdata/

  #run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           #"-indexname sfx -dna -tis -suf -lcp -ssp -pl"
  #run_test "#{$bin}gt repfind -l 8 -ii sfx"                # once for normal
  #run "grep -v '^#' #{$last_stdout}"     # from stdout_2 invert-match save o stdout_3    
  ##puts $last_stdout       
  #run "diff -w #{$last_stdout} #{$testdata}repfind-8-Atinsert.txt"  # w Ignore all white space
  #run_test "#{$bin}gt repfind -scan -l 8 -ii sfx"
  #run "grep -v '^#' #{$last_stdout}"
  #run "diff -w #{$last_stdout} #{$testdata}repfind-8-Atinsert.txt"  # once for -scan
  #run_test "#{$bin}gt repfind -samples 40 -l 6 -ii sfx",:maxtime => 600  # test time not exceed the limited time
#end

#allfiles.each do |reffile|
  #allfiles.each do |queryfile|
    #if queryfile != reffile
      #Name "gt greedyfwdmat #{reffile} #{queryfile}"
      #Keywords "gt_greedyfwdmat small"
      #Test do
        #createandcheckgreedyfwdmat("#{$testdata}/#{reffile}",
                                   #"#{$testdata}/#{queryfile}")
        #checktagerator("#{$testdata}/#{reffile}",
                       #"#{$testdata}/#{queryfile}")
        #run "rm -f sfx.* fmi.* pck.*"
      #end
    #end
  #end
#end

#allfiles.each do |reffile|
  #allfiles.each do |queryfile|
    #if queryfile != reffile
      #Name "gt idxlocali #{reffile} #{queryfile}"
      #Keywords "gt_idxlocali"
      #Test do
        #run("#{$bin}gt packedindex mkindex -ssp -tis -indexname pck -db " +
            #"#{$testdata}/#{reffile} -sprank -dna -pl -bsize 10 " +
            #"-locfreq 32 -dir rev",
            #:maxtime => 100)
        #run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck " +
                 #"-q #{$testdata}/#{queryfile}",
                 #:maxtime => 100)
        #run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck -online " +
                 #"-q #{$testdata}/#{queryfile}",
                 #:maxtime => 100)
        #run_test "#{$bin}gt suffixerator -indexname sfx -ssp -tis -suf -dna " +
                 #"-v -db #{$testdata}/#{reffile}"
        #run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx " +
                 #"-q #{$testdata}/#{queryfile}",
                 #:maxtime => 100)
        #run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx -online " +
                 #"-q #{$testdata}/#{queryfile}",
                 #:maxtime => 100)
      #end
    #end
  #end
#end

#allfiles.each do |reffile|
  #Name "gt packedindex #{reffile}"
  #Keywords "gt_packedindex small"
  #Test do
    #run_test("#{$bin}gt packedindex mkindex -tis -ssp -indexname pck " +
             #"-sprank -db #{$testdata}/#{reffile} -dna -pl -bsize 10 " +
             #" -locfreq 32 -dir rev",
             #:maxtime => 1200)
  #end
#end





def checkmaxmat4withrepfind(reffile,queryfile,minlength)
  #reffilepath=addfilepath(reffile)
  #queryfilepath=addfilepath(queryfile)
  reffilepath="#{$testdata}#{reffile}"
  queryfilepath="#{$testdata}#{queryfile}"
  
  # generate repfind result
  idxname=reffile + "-idx"
  run_test "#{$bin}gt suffixerator -algbds 3 40 120 -db " +
           "#{reffilepath} -indexname #{idxname} -dna -suf -tis -lcp -ssp -pl"
  run_test("#{$bin}gt repfind -v -l #{minlength} -ii #{idxname} -q #{queryfilepath}",
           :maxtime => 320)
  run "sort #{$last_stdout}"
  run "mv #{$last_stdout} repfind.result"
    
  # generate maxmat4 result 
  pckname=reffile + "-pck"                                                 
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{reffilepath} -indexname #{pckname} -sprank -dna -ssp -des -sds -pl"  
  #run "#{$bin}gt suffixerator -indexname sfx -tis -suf -ssp -dna -v " +
           #"-db #{reffile}"
  #run("#{$bin}gt packedindex mkindex -tis -ssp -indexname pck -db #{reffile} " +
      #"-sprank -dna -pl -bsize 10 -locfreq 32 -dir rev", :maxtime => 100)
      
  #run "#{$bin}gt prebwt -maxdepth 4 -pck pck"
  #run "#{$bin}gt matstat -verify(uniquesub) -output querypos -min 1 -max 20 -query #{queryfile} -pck pck"
  run_test("#{$bin}gt dev maxmat4 -b -l #{minlength} -L -s -c #{pckname} #{queryfilepath}", :maxtime => 320)
  run "sed -e '/^>/d' #{$last_stdout} | sort"            
  run "mv #{$last_stdout} maxmat4.result"
  
  # compare both results
  #run "cmp -s repfind.result maxmat4.result"
  run "diff -w repfind.result maxmat4.result"
  
  # remove both results
  #run "rm -f sfx.* fmi.* pck.*"
end

Name "check maxmat4 at1MB versus U89959_genomic.fas comparing with repfind"
Keywords "check_maxmat4_with_repfind"
Test do
  checkmaxmat4withrepfind("at1MB","U89959_genomic.fas",15)
end
