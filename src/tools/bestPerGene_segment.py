#Difference with v1. (1) Replace M2,M1,NEP etc. (2) one line for "Segment x" replace. (3) Include one result per each HxNx subtype.

import os, sys

segments=["(HA)", "(MP)", "(NA)", "(NP)", "(NS)", "(PA)", "(PB1)", "(PB2)"]

def ideclade(s):
	#print(s)
	return "(H"+s.split("(H")[1].split(")")[0]

genes={}
with open (sys.argv[1],'r') as fi:
	for line in fi:
		#if "complete" in line and "partial" not in line:
		line=line.replace("RNA","")
		line=line.replace(":","")
		line=line.replace(","," ")
		line=line.replace("M2","(MP)").replace("M1","(MP)").replace("NEP","(NS)").replace("NS1","(NS)").replace("BM2","(MP)").replace(" HA ","(HA)").replace(" MP ","(MP)").replace(" NA ","(NA)").replace(" NP ","(NP)").replace(" NS ","(NS)").replace(" PA ","(PA)").replace(" PB1 ","(PB1)").replace(" PB2 ","(PB2)")
                #In case not segment is present, it replaces "segment X" substring
		if not any(x in line for x in segments):
			line=line.replace("segment 1","(PB2)").replace("segment 2","(PB1)").replace("segment 3","(PA)").replace("segment 4","(HA)").replace("segment 5","(NP)").replace("segment 6","(NA)").replace("segment 7","(MP)").replace("segment 8","(NS)").replace("seg 1","(PB2)").replace("seg 2","(PB1)").replace("seg 3","(PA)").replace("seg 4","(HA)").replace("seg 5","(NP)").replace("seg 6","(NA)").replace("seg 7","(MP)").replace("seg 8","(NS)")
		#Group results by segment
		if "B virus" not in line:
			subtype=ideclade(line)
		else:
			subtype="B virus"
		for segment in segments:
			#Only check subtype for HA and NA
			if (segment == "(HA)" or segment == "(NA)") and "B virus" not in line:
				if segment in line.split("\t")[9].split(subtype)[1]:
					if segment+"_"+subtype in genes.keys():
						genes[segment+"_"+subtype]+=line
					else:
						genes[segment+"_"+subtype]=line
					break
			else:
				if segment in line.split("\t")[9].split(subtype)[1]:
                               		if segment in genes.keys():
                                       		genes[segment]+=line
                               		else:
                                      		genes[segment]=line
                               		break
#print (genes["(MP)"])
#print ("_-_-_-_-_")
#Take the longest per each segment and per each subtype!
res=""
#contig.00017    gi|1732444596|gb|LC497148|Influenza     0.0     1853    1732444596      LC497148        N/A     N/A     N/A     A virus (A/duck/Vietnam/HU9-521/2018(H6N6)) viral cRNA, segment 7, complete sequence    1006    1     $
#contig.00019    gi|1025607882|gb|LC148831|Influenza     0.0     1454    1025607882      LC148831        N/A     N/A     N/A     A virus (A/duck/Hokkaido/W9/2015(H1N1)) viral cRNA, segment: 7, complete sequence

for k in genes.keys():
	maxi=0
	maxiline=""
	for line in genes[k].split("\n"):
		if line!='':
			if int(line.split("\t")[10])>maxi:
				maxi=int(line.split("\t")[10])
				maxiline=line
	res+=maxiline+"\n"

print(res)

