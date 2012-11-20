#!/usr/bin/env ruby

require 'rubygems'
require 'rphylip'


phy = Phylip.new "mt-co1_Rotifera.phylip"
phy.names.each do |genbank_id|
  puts genbank_id 
  url = "http://www.ncbi.nlm.nih.gov/nuccore/#{genbank_id}"
  Dir.chdir("genbank_entries") do
    system "wget #{url}"
  end
end
