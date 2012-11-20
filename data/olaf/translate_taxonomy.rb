#!/usr/bin/env ruby

require 'rubygems'
require 'rnewick'

# Makes a new newick string translating the numbers (ncbi) into the names (genbank) according to "taxalist.txt"

taxalist = "taxalist.txt"
taxonomy = "taxonomy.out"
usage = "#{$0} #{taxonomy} #{taxalist}"
raise usage unless ARGV.size == 2
# parse the dic of names
genbank_names = {}
File.open(taxalist).readlines.each do |line|
  genbank, ncbi = line.split
  genbank_names[ncbi] = genbank
end
# translate names
t = NewickFile.new(taxonomy).newickStrings.first
t.str.gsub!(/([0-9]+)/) do |n| # assume taxa names are numbers and we have a  simple topology! 
  genbank_names[n]
end
puts t.str
