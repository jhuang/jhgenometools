#allfiles = ["Atinsert.fna",
            #"Duplicate.fna",
            #"Random-Small.fna",
            #"Random.fna",
            #"Random159.fna",
            #"Random160.fna",
            #"RandomN.fna",
            #"TTT-small.fna",
            #"trna_glutamine.fna"]
allfiles = ["Atinsert.fna",
            "Random-Small.fna",
            "Random159.fna",
            "RandomN.fna",
            "trna_glutamine.fna"]


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


# format repfind output to maxmat4 output 
def formatrepfindoutput(filename)
	if (File.exist?(filename))
		begin 
			file = File.new(filename, "r")    
		rescue
			STDERR.print "Could not open file #{filename}!\n"
			exit 1
		end
	else
		STDERR.print "File #{filename} does not exist!\n"
		exit 1
	end

	data = file.readlines() 
	output = File.new(filename, "w")
	data.each do |line|
		line.strip!
		line_data = line.split(%r{\s+})
			
		output.printf("%d %d %s\n",line_data[2].to_i+1, line_data[6].to_i+1, line_data[0])
	end
	output.close
end


def formatmaxmat4output(filename)
	if (File.exist?(filename))
		begin 
			file = File.new(filename, "r")    
		rescue
			STDERR.print "Could not open file #{filename}!\n"
			exit 1
		end
	else
		STDERR.print "File #{filename} does not exist!\n"
		exit 1
	end

	data = file.readlines() 
	output = File.new(filename, "w")
	data.each do |line|
		line.strip!
		line_data = line.split(%r{\s+})
			
	  if (line_data.length==4) then		
		  output.printf("%d %d %d\n",line_data[1].to_i, line_data[2].to_i, line_data[3].to_i)
		elsif (line_data.length==3) then	
		  output.printf("%d %d %d\n",line_data[0].to_i, line_data[1].to_i, line_data[2].to_i)
		end	
	end
	output.close
end


def checkmaxmat4withrepfind(reffile,queryfile,minlength)
  #reffilepath=addfilepath(reffile)
  #queryfilepath=addfilepath(queryfile)
  reffilepath="#{$testdata}#{reffile}"
  queryfilepath="#{$testdata}#{queryfile}"
  #reffile="peakseq.fa"
  #queryfile="query.fa"
  #reffilepath="/home/jiabin/fasta/peakseq.fa"
  #queryfilepath="/home/jiabin/fasta/query.fa"
  
  # generate repfind result
  idxname=reffile + "-idx"
  run_test "#{$bin}gt suffixerator -algbds 3 40 120 -db " +
           "#{reffilepath} -indexname #{idxname} -dna -suf -tis -lcp -ssp -pl"
  run_test("#{$bin}gt repfind -l #{minlength} -ii #{idxname} -q #{queryfilepath}",
           :maxtime => 320)
  formatrepfindoutput("#{$last_stdout}")
  run "sed -e '/^\s*$/d' #{$last_stdout} | sort"
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
  run_test("#{$bin}gt dev maxmat4 -maxmatch -l #{minlength} -L #{pckname} #{queryfilepath}", :maxtime => 320)
  formatmaxmat4output("#{$last_stdout}")
  run "sed -e '/^>/d' #{$last_stdout} | sort"            
  run "mv #{$last_stdout} maxmat4.result"
  
  # compare both results
  #run "cmp -s repfind.result maxmat4.result"
  run "diff -w repfind.result maxmat4.result"
  
  # remove both results
  #run "rm -f repfind.result maxmat4.result"
end

#Name "check maxmat4 at1MB versus U89959_genomic.fas comparing with repfind"
#Keywords "check_maxmat4_with_repfind"
#Test do
  #checkmaxmat4withrepfind("at1MB","U89959_genomic.fas",15)
#end


allfiles.each do |reffile|
  allfiles.each do |queryfile|
    [12,15,18,21,24].each do |minlength| 
			if queryfile != reffile
				Name "check maxmat4 -l #{minlength} #{reffile} versus #{queryfile} comparing with repfind"
				Keywords "check_maxmat4_with_repfind"
				Test do
					checkmaxmat4withrepfind("at1MB","U89959_genomic.fas",minlength)
				end
      end
    end
  end
end
