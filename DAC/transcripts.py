#!/usr/bin/env python
import gff
import sys
from collections import defaultdict

reader =  gff.Reader(open(sys.argv[1]))

#Summing the lengths of transcripts should be something like:


transcript_lengths = defaultdict(int)
for record in reader:
	if record.type == "exon" or record.type == "tRNAscan":
		transcript_id = record.attributes["transcript_id"][0]
		exon_length = record.end - record.start + 1
		transcript_lengths[transcript_id] += exon_length

for tnx in transcript_lengths.keys():
	print "%s\t%s" % (tnx, transcript_lengths[tnx])

