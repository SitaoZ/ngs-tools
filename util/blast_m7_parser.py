#!/usr/bin/python 

import re,sys


"""
Author : Zhu Sitao
Date : 2018-3-30
Dest : convert blast's result from -m7 (xml) to -m8 (tabular, without comment lines)
 
"""

if len(sys.argv) != 3:
	print ("usage: python blast_m7_parser.py <InFile,blast.xml> <OutFile,blast.tabular>")
	sys.exit(1)

InFile,OutFile = sys.argv[1],sys.argv[2]

OutPath = open(OutFile,'w')
header = "Query_id\tSubject_id\tIdentity\tAlign_length\tMiss_match\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE_value\tScore\tSubject_annotation\n"
OutPath.writelines(header)
#print (header)
# pattern description
queryPatt = re.compile(r'<(Iteration_query-def)>(\S+).*<\/\1')
hitIdPatt = re.compile(r'<(Hit_id)>(.*)<\/\1')
hitDefPatt = re.compile(r'<(Hit_def)>([^\s]*)\s+(.*)<\/\1')
bitsPatt = re.compile(r'<(Hsp_bit-score)>(.*)<\/\1')
evaluePatt = re.compile(r'<(Hsp_evalue)>(.*)<\/\1')
qFromPatt = re.compile(r'<(Hsp_query-from)>(.*)<\/\1')
qToPatt = re.compile(r'<(Hsp_query-to)>(.*)<\/\1')
hFromPatt = re.compile(r'<(Hsp_hit-from)>(.*)<\/\1')
hToPatt = re.compile(r'<(Hsp_hit-to)>(.*)<\/\1')
framePatt = re.compile(r'<(Hsp_query-frame)>(.*)<\/\1') ## Determine the direction of alignment
identityPatt = re.compile(r'<(Hsp_identity)>(.*)<\/\1')
gapPatt = re.compile(r'<(Hsp_gaps)>(.*)<\/\1') ## maybe not exist ;so init gap=0 firstly
lengthPatt = re.compile(r'<(Hsp_align-len)>(.*)<\/\1')
endPatt = re.compile(r'<\/Hsp>') ## match on hit end

with open(InFile,'r') as F:
	# init gap=0
	gap = 0
	for line in F:
		line = line.strip()
		if queryPatt.match(line):
			query = queryPatt.search(line).group(2)
		elif hitIdPatt.match(line):
			hitId = hitIdPatt.match(line).group(2)
		elif hitDefPatt.match(line):
			hitDef = hitDefPatt.match(line).group(2)
			annot = hitDefPatt.match(line).group(3) if hitDefPatt.match(line).group(3) and len(hitDefPatt.match(line).group(3)) != 0 else "none"
		elif bitsPatt.match(line):
			bits = bitsPatt.match(line).group(2)
		elif evaluePatt.match(line):
			evalue = evaluePatt.match(line).group(2)
		elif qFromPatt.match(line):
			qFrom = qFromPatt.match(line).group(2)
		elif qToPatt.match(line):
			qTo = qToPatt.match(line).group(2)
		elif hFromPatt.match(line):
			hFrom = hFromPatt.match(line).group(2)
		elif hToPatt.match(line):
			hTo = hToPatt.match(line).group(2)
		elif framePatt.match(line):
			frame = framePatt.match(line).group(2)
		elif identityPatt.match(line):
			identity = identityPatt.match(line).group(2)
		elif gapPatt.match(line):
			gap = gapPatt.match(line).group(2)
		elif lengthPatt.match(line):
			length = lengthPatt.match(line).group(2)
		elif endPatt.match(line):
			percent = str("%.2f" % (int(identity)/int(length)*100))
			hit = hitId if 'gi' in hitId else hitDef # nr database hitID start with "gi"
			start_end_format = qFrom+"\t"+qTo if int(frame)>0 else qTo+"\t"+qFrom
			misMatch = str(int(length)-int(identity))
## Query_id,Subject_id,Identity,Align_length,Miss_match,Gap,Query_start,Query_end,Subject_start,Subject_end,E_value,Score,Subject_annotation
			oneHitOutput = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(query,hit,percent,length,misMatch,gap,start_end_format,hFrom,hTo,evalue,bits,annot)
			OutPath.writelines(oneHitOutput)
OutPath.close()
