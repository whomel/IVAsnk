import os, sys

#script to complete the (incompletes) headers from a `bedtools maskfasta` output.
#python fix_bedtool_id.py {input.blastFasta} {output.blastFastaMasked}

#dictionary 
with open (sys.argv[1],'r') as fi:
	recdict={}
	for line in fi:
		if line[0] ==">":
			recdict[line[1:].split(" ")[0]]=line

with open (sys.argv[2],'r') as fi:
	res=""
	for line in fi:
		if line[0] == ">":
			res+=recdict[line[1:-1]]
		else:
			res+=line

print(res)

