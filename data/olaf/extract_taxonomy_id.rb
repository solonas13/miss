#!/usr/bin/env ruby

require 'rubygems'

dirname = "genbank_entries"
Dir.entries(dirname).select{|f| not f =~ /^\./}.each do |f|
  line = File.open(File.join dirname, f).readlines.select{|l| l =~ /ORGA/}.first
  if line =~ /ORGANISM=([0-9]+)/
    puts "#{f} #{$1}" 
  end
end
