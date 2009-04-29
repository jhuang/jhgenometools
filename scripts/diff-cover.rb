#!/usr/bin/env ruby

def copy(v,r)
  return Array.new(r) {v}
end

def partialsums(bseq)
  sum = 0
  l = [0]
  bseq.each do |b|
    sum += b
    l = l.push(sum)
  end
  return l
end

def checkcover(v,psum,i)
  psum.each do |j|
    psum.each do |k|
      if (k-j) % v == i
        return true
      end
    end
  end
  return false
end

def checkallcover(v,psum)
  0.upto(v-1) do |i|
    if not checkcover(v,psum,i)
      STDERR.puts "cannot find difference elements for #{i} in #{v}-cover"
      exit(1)
    end
  end
end

def diff_cover(logmod)
  dc = case logmod 
    when 0 then [0]
    when 1 then [0,1]
    when 2 then [0,1,2]
    when 3 then [0,1,2,4]
    when 4 then [0,1,2,5,8]
    when 5 then [0,1,2,3,7,11,19]
    when 6 then [0,1,2,5,14,16,34,42,59]
    when 7 then [0,1,3,7,17,40,55,64,75,85,104,109,117]
    when 8 then [0,1,3,7,12,20,30,44,65,80,89,96,114,
                122,128,150,196,197,201,219]
    else
      v = 2**logmod
      r = 0
      while 24*r*r+36*r+13 < v
        r += 1
      end
      bseq = copy(1,r) + 
             copy(r+1,1) + 
             copy(2*r+1,r) + 
             copy(4*r+3,2*r+1) + 
             copy(2*r+2,r+1) + 
             copy(1,r)
      psum = partialsums(bseq)
      checkallcover(v,psum)
      partialsums(bseq)
  end
  return dc
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <maxcover>"
  exit(1)
end

maxcover=ARGV[0].to_i

puts "static Diffvalue differencecovertab[] ="
puts "{"
lentab = []
0.upto(maxcover) do |logmod|
  dc = diff_cover(logmod)
  lentab = lentab.push(dc.length)
  print "  /* #{2**logmod} */ " + dc.join("UL,")
  if logmod == maxcover
    puts "UL"
  else
    puts "UL,"
  end
end
puts "};"
puts ""

puts "static unsigned int differencecoversizes[]"
puts "  = {" + lentab.join("U,") + "U};"
