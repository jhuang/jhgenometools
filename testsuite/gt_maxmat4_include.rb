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


def checkmaxmat4withmummer(reffile,queryfile,matchmode)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)

  pckname=reffile + "-pck"                                                 
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{reffilepath} -indexname #{pckname} -sprank -dna -ssp -des -sds -pl"   #-sprank -dna -tis -ssp -des -sds -pl
  run_test("#{$bin}gt dev maxmat4 -#{matchmode} -b -l #{determinemaxmatchminlength(reffile)} -L -s -c #{pckname} #{queryfilepath}", :maxtime => 320)
  #run "grep -v '^>' #{$last_stdout} | sort"  
  if (matchmode=="mumreference" || matchmode=="mum") 
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
  if (matchmode=="mumreference" || matchmode=="mum") 
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
        Name "store mummer -mumref #{reffile}/#{queryfile} result"
				Keywords "generate_mummer_mumreference_benchmark"
				Test do
					generatemummerresults(reffile,queryfile,"mumreference")
				end
        Name "store mummer -maxmatch #{reffile}/#{queryfile} result"
				Keywords "generate_mummer_maxmatch_benchmark"
				Test do
					generatemummerresults(reffile,queryfile,"maxmatch")
				end
        Name "store mummer -mum #{reffile}/#{queryfile} result"
				Keywords "generate_mummer_mum_benchmark"
				Test do
					generatemummerresults(reffile,queryfile,"mum")
				end
				
        Name "maxmat4 -mumref #{reffile}/#{queryfile} against mummer"
        Keywords "check_maxmat4_mumreference_with_mummer"
        Test do
          checkmaxmat4withmummer(reffile,queryfile,"mumreference")
        end
        Name "maxmat4 -maxmatch #{reffile}/#{queryfile} against mummer"
        Keywords "check_maxmat4_maxmatch_with_mummer"
        Test do
          checkmaxmat4withmummer(reffile,queryfile,"maxmatch")
        end   
        Name "maxmat4 -mum #{reffile}/#{queryfile} against mummer"
        Keywords "check_maxmat4_mum_with_mummer"
        Test do
          checkmaxmat4withmummer(reffile,queryfile,"mum")
        end      
      end
    end
  end
end

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
  reffilepath="#{$testdata}#{reffile}"
  queryfilepath="#{$testdata}#{queryfile}"
  
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
  run_test("#{$bin}gt dev maxmat4 -maxmatch -l #{minlength} -L #{pckname} #{queryfilepath}", :maxtime => 320)
  formatmaxmat4output("#{$last_stdout}")
  run "sed -e '/^>/d' #{$last_stdout} | sort"            
  run "mv #{$last_stdout} maxmat4.result"
  
  # compare both results
  run "diff -w repfind.result maxmat4.result"
end


allfiles.each do |reffile|
  allfiles.each do |queryfile|
    [12,15,18,21,24].each do |minlength| 
			if queryfile != reffile
				Name "maxmat4 -l #{minlength} #{reffile}/#{queryfile} (repfind)"
				Keywords "check_maxmat4_with_repfind"
				Test do
					#checkmaxmat4withrepfind("at1MB","U89959_genomic.fas",minlength)
					checkmaxmat4withrepfind(reffile,queryfile,minlength)
				end
      end
    end
  end
end
