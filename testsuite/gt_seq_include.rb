Name "gt seq nonexistent file"
Keywords "gt_seq"
Test do
  run_test("#{$bin}gt seq #{$testdata}nonexistent_file", :retval => 1)
end

1.upto(7) do |i|
  Name "gt seq fail #{i}"
  Keywords "gt_seq"
  Test do
    run_test("#{$bin}gt seq -recreate #{$testdata}gt_bioseq_fail_#{i}.fas",
             :retval => 1)
  end
end

1.upto(7) do |i|
  Name "gt seq fail #{i} (stdin)"
  Keywords "gt_seq"
  Test do
    run_test("cat #{$testdata}gt_bioseq_fail_#{i}.fas | #{$bin}gt seq -recreate -", :retval => 1)
  end
end

Name "gt seq test 1"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -recreate #{$testdata}gt_bioseq_succ_1.fas"
  if not File.exists?("#{$testdata}gt_bioseq_succ_1.fas.gt_bsi") then
    raise TestFailed, "file \"#{$testdata}gt_bioseq_succ_1.fas.gt_bsi\" does not exist"
  end
  if not File.exists?("#{$testdata}gt_bioseq_succ_1.fas.gt_bsr") then
    raise TestFailed, "file \"#{$testdata}gt_bioseq_succ_1.fas.gt_bsr\" does not exist"
  end
  old_bsi_mtime = File.mtime("#{$testdata}gt_bioseq_succ_1.fas.gt_bsi")
  old_bsr_mtime = File.mtime("#{$testdata}gt_bioseq_succ_1.fas.gt_bsr")
  sleep(1)
  run_test "#{$bin}gt seq #{$testdata}gt_bioseq_succ_1.fas"
  new_bsi_mtime = File.mtime("#{$testdata}gt_bioseq_succ_1.fas.gt_bsi")
  new_bsr_mtime = File.mtime("#{$testdata}gt_bioseq_succ_1.fas.gt_bsr")
  # make sure the index file have not been recreated
  if (old_bsi_mtime != new_bsi_mtime) or (old_bsr_mtime != new_bsr_mtime) then
    raise TestFailed, "index files have been recreated"
  end
end

Name "gt seq test 1 (stdin)"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_1.fas | #{$memcheck} #{$bin}gt seq -"
end

Name "gt seq test 2"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -recreate #{$testdata}gt_bioseq_succ_2.fas"
  if not File.exists?("#{$testdata}gt_bioseq_succ_2.fas.gt_bsi") then
    raise TestFailed, "file \"#{$testdata}gt_bioseq_succ_2.fas.gt_bsi\" does not exist"
  end
  if not File.exists?("#{$testdata}gt_bioseq_succ_2.fas.gt_bsr") then
    raise TestFailed, "file \"#{$testdata}gt_bioseq_succ_2.fas.gt_bsr\" does not exist"
  end
  old_bsi_mtime = File.mtime("#{$testdata}gt_bioseq_succ_2.fas.gt_bsi")
  old_bsr_mtime = File.mtime("#{$testdata}gt_bioseq_succ_2.fas.gt_bsr")
  sleep(1)
  run_test "#{$bin}gt seq #{$testdata}gt_bioseq_succ_2.fas"
  new_bsi_mtime = File.mtime("#{$testdata}gt_bioseq_succ_2.fas.gt_bsi")
  new_bsr_mtime = File.mtime("#{$testdata}gt_bioseq_succ_2.fas.gt_bsr")
  # make sure the index file have not been recreated
  if (old_bsi_mtime != new_bsi_mtime) or (old_bsr_mtime != new_bsr_mtime) then
    raise TestFailed, "index files have been recreated"
  end
end

Name "gt seq test 2 (stdin)"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_2.fas | #{$memcheck} #{$bin}gt seq -"
end

Name "gt seq test 3"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -recreate -showfasta -width 70 #{$testdata}gt_bioseq_succ_3.fas"
  run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.fas"
end

Name "gt seq test 3 (stdin)"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_3.fas | #{$memcheck} #{$bin}gt seq -recreate -showfasta -width 70 -"
  run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.fas"
end

1.upto(3) do |i|
  Name "gt seq test 3 out #{i}"
  Keywords "gt_seq"
  Test do
    run_test "#{$bin}gt seq -showseqnum #{i} -width 70 #{$testdata}gt_bioseq_succ_3.fas"
    run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.out#{i}"
  end
end

1.upto(3) do |i|
  Name "gt seq test 3 out #{i} (stdin)"
  Keywords "gt_seq"
  Test do
    run "cat #{$testdata}gt_bioseq_succ_3.fas | #{$memcheck} #{$bin}gt seq -showseqnum #{i} -width 70 -"
    run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.out#{i}"
  end
end

Name "gt seq test 3 out 4 fail"
Keywords "gt_seq"
Test do
  run_test("#{$bin}gt seq -showseqnum 4 #{$testdata}gt_bioseq_succ_3.fas",
           :retval => 1)
end

Name "gt seq test 3 stat"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -stat #{$testdata}gt_bioseq_succ_3.fas"
end

Name "gt seq test 3 stat (stdin)"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_3.fas | #{$memcheck} #{$bin}gt seq -stat -"
end

Name "gt seq test multiple sequence files"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -recreate #{$testdata}gt_bioseq_succ_1.fas #{$testdata}gt_bioseq_succ_2.fas"
end

Name "gt seq test multiple sequence files (incl. stdin)"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_2.fas | #{$memcheck} #{$bin}gt seq -recreate #{$testdata}gt_bioseq_succ_1.fas -"
end

Name "gt seq -gc-content"
Keywords "gt_seq"
Test do
  run "cat #{$testdata}gt_bioseq_succ_3.fas | #{$memcheck} #{$bin}gt seq -gc-content -"
  run "diff #{$last_stdout} #{$testdata}gt_bioseq_succ_3.gc"
end

Name "gt seq -seqlengthdistri"
Keywords "gt_seq"
Test do
  run_test "#{$bin}gt seq -seqlengthdistri #{$testdata}sw100K1.fsa"
  run "diff #{$last_stdout} #{$testdata}gt_bioseq_seqlengthdistri.out"
end

Name "gt seq -showseqnum 1.5"
Keywords "gt_seq"
Test do
  run_test("#{$bin}gt seq -showseqnum 1.5 #{$testdata}sw100K1.fsa",
           :retval => 1)
end
