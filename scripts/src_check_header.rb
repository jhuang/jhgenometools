#!/usr/bin/ruby -w

require 'fileutils'

had_err = false

ARGV.each do |file|
  if !File.exist?(file)
    puts "file #{file} does not exist!"
    had_err = true
  end
  if !had_err
    File.open("src/header_test.c", 'w') do |cfile|
      cfile.printf("#include \"%s\"\nint main()\n{" +
                   "  return 0;\n}", file.sub(/.*src\//, ''))
    end
    err = `make obj/src/header_test.o`
    if $? != 0
      puts "missing inclusions in #{file}!"
      puts err
      had_err = true
    end
    FileUtils.rm_f(["obj/src/header_test.o",
                    "src/header_test.c"])
  end
end
exit 0 if !had_err
exit 1
