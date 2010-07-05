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
  if (matchmode=="mumreference") 
			run_test("#{$bin}gt dev maxmat4 -#{matchmode} -b -l #{determinemaxmatchminlength(reffile)} -L -s -c #{pckname} #{queryfilepath}", :maxtime => 32000)
			#run "grep -v '^>' #{$last_stdout} | sort" 
      run "sed -e '/^>/d' #{$last_stdout}"
  elsif (matchmode=="mum")
			run_test("#{$bin}gt dev maxmat4 -#{matchmode} -r -l #{determinemaxmatchminlength(reffile)} -s -c #{pckname} #{queryfilepath}", :maxtime => 32000)
      run "sed -e '/^>/d' #{$last_stdout}"
  elsif (matchmode=="maxmatch") 
  		run_test("#{$bin}gt dev maxmat4 -#{matchmode} -l #{determinemaxmatchminlength(reffile)} -L -s #{pckname} #{queryfilepath}", :maxtime => 32000)
      run "sed -e '/^>/d' #{$last_stdout} | sort"
  end             
  run "diff -w #{$last_stdout} #{$gttestdata}maxmat4-result/#{reffile}-#{queryfile}_#{matchmode}.result"
end


def generatemummerresults(reffile, queryfile, matchmode)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)
  if (matchmode=="mumreference") 
			run "~/maxmat3/mm3src/maxmat3.x -#{matchmode} -b -l #{determinemaxmatchminlength(reffile)} -L -s -c -n #{reffilepath} #{queryfilepath}"
			#run "grep -v '^>' #{$last_stdout} | sort" 
      run "sed -e '/^>/d' #{$last_stdout}"       
  elsif (matchmode=="mum") 
 			run "~/maxmat3/mm3src/maxmat3.x -#{matchmode} -r -l #{determinemaxmatchminlength(reffile)} -s -c -n #{reffilepath} #{queryfilepath}"
      run "sed -e '/^>/d' #{$last_stdout}"          
  elsif (matchmode=="maxmatch") 
   	  run "~/maxmat3/mm3src/maxmat3.x -#{matchmode} -l #{determinemaxmatchminlength(reffile)} -L -s -n #{reffilepath} #{queryfilepath}"
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
           :maxtime => 32000)
  formatrepfindoutput("#{$last_stdout}")
  run "sed -e '/^\s*$/d' #{$last_stdout} | sort"
  run "mv #{$last_stdout} repfind.result"
  
  # generate maxmat4 result 
  pckname=reffile + "-pck"                                                 
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{reffilepath} -indexname #{pckname} -sprank -dna -ssp -des -sds -pl"      
  run_test("#{$bin}gt dev maxmat4 -maxmatch -l #{minlength} -L #{pckname} #{queryfilepath}", :maxtime => 32000)
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


Name "gt dev maxmat4 -help"
Keywords "gt_maxmat4_inputerror"
Test do
  run_test "#{$bin}gt dev maxmat4 -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt dev maxmat4 -version"
Keywords "gt_maxmat4_inputerror"
Test do
  run_test "#{$bin}gt dev maxmat4 -version"
  grep $last_stdout, "Compile flags"
end

Name "gt dev maxmat4 -noop"
Keywords "gt_maxmat4_inputerror"
Test do
  run_test "#{$bin}gt dev maxmat4 -noop", :retval => 1
  grep $last_stderr, "unknown option"
end

Name "gt dev maxmat4 referencefile (packed index) not exist"
Keywords "gt_maxmat4_inputerror"
Test do                 
	run_test "#{$bin}gt dev maxmat4 -maxmatch /nothing/pck #{$testdata}at1MB", :retval => 1                               
  grep $last_stderr, "No such file or directory"
end

Name "gt dev maxmat4 referencefile (packed index) empty"
Keywords "gt_maxmat4_inputerror"
Test do                 
	run_test "#{$bin}gt dev maxmat4 -maxmatch #{$testdata}empty_file #{$testdata}at1MB", :retval => 1                           
  grep $last_stderr, "No such file or directory"
end

Name "gt dev maxmat4 referencefile (packed index) incomplete"
Keywords "gt_maxmat4_inputerror"
Test do       
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -pl", :maxtime => 32000           
	run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck #{$testdata}at1MB", :retval => 1, :maxtime => 32000                         
  grep $last_stderr, "No such file or directory"
end

Name "gt dev maxmat4 referencefile corrupt"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}corrupt.fas -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :retval => 1   
  grep $last_stderr, "cannot guess file type of file"
end

Name "gt dev maxmat4 referencefile protein"
Keywords "gt_maxmat4_inputerror"
Test do       
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}sw100K1.fsa -indexname sw100K1pck -sprank -dna -ssp -des -sds -pl", :retval => 1                       
  grep $last_stderr, "illegal character"
end

Name "gt dev maxmat4 queryfile not exist"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck /nothing", :retval => 1, :maxtime => 32000   
  grep $last_stderr, "No such file or directory"
end

Name "gt dev maxmat4 queryfile empty"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck #{$testdata}empty_file", :retval => 1, :maxtime => 32000   
  grep $last_stderr, "unknown file contents"
end

Name "gt dev maxmat4 queryfile corrupt"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck #{$testdata}corrupt.fas", :retval => 1, :maxtime => 32000  
  grep $last_stderr, "cannot guess file type of file"
end

Name "gt dev maxmat4 queryfile protein"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck #{$testdata}sw100K1.fsa", :retval => 1, :maxtime => 32000   
  grep $last_stderr, "illegal character"
end

Name "gt dev maxmat4 more than one queryfile"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch at1MBpck #{$testdata}Random.fna #{$testdata}RandomN.fna", :maxtime => 32000   
end

Name "gt dev maxmat4 set option -s twice"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch -s -s at1MBpck #{$testdata}at1MB", :retval => 1, :maxtime => 32000 
  grep $last_stderr, "option \"s\" already set"
end

Name "gt dev maxmat4 set option -c without -b or -r"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch -c at1MBpck #{$testdata}at1MB", :retval => 1, :maxtime => 32000  
  grep $last_stderr, "option \"-c\" requires option \"-b\" or \"-r\""
end

Name "gt dev maxmat4 -b -r exclude each other"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch -b -r at1MBpck #{$testdata}at1MB", :retval => 1, :maxtime => 32000  
  grep $last_stderr, "option \"-b\" and option \"-r\" exclude each other"
end

Name "gt dev maxmat4 -maxmatch -mum exclude each other"
Keywords "gt_maxmat4_inputerror"
Test do                                               
  run_test "#{$bin}gt packedindex mkindex -bsize 10 -locfreq 8 -dir rev -db " +
           "#{$testdata}at1MB -indexname at1MBpck -sprank -dna -ssp -des -sds -pl", :maxtime => 32000 
  run_test "#{$bin}gt dev maxmat4 -maxmatch -mum at1MBpck #{$testdata}at1MB", :retval => 1, :maxtime => 32000  
  grep $last_stderr, "option \"-mum\" and option \"-maxmatch\" exclude each other"
end
