import os, sys

segments=["(HA)", "(MP)", "(NA)", "(NP)", "(NS)", "(PA)", "(PB1)", "(PB2)","(M2)","(M1)","(NEP)","(NS1)","(BM2)"," NEP "," HA ", " MP ", " NA ", " NP ", " NS ", " PA ", " PB1 ", " PB2 "," M2 "," M1 "," NEP "," NS1 "," BM2 "]

genes={}
with open (sys.argv[1],'r') as fi:
	for line in fi:
		#print (line)
		line=line.replace("RNA","")
		line=line.replace(":","")
		line=line.replace(","," ")
                #In case not segment is present, it replaces segment substring
		if not any(x in line for x in segments):
			line=line.replace("segment 1","(PB2)")
			line=line.replace("segment 2","(PB1)")
			line=line.replace("segment 3","(PA)")
			line=line.replace("segment 4","(HA)")
			line=line.replace("segment 5","(NP)")
			line=line.replace("segment 6","(NA)")
			line=line.replace("segment 7","(MP)")
			line=line.replace("segment 8","(NS)")
			line=line.replace("seg 1","(PB2)")
			line=line.replace("seg 2","(PB1)")
			line=line.replace("seg 3","(PA)")
			line=line.replace("seg 4","(HA)")
			line=line.replace("seg 5","(NP)")
			line=line.replace("seg 6","(NA)")
			line=line.replace("seg 7","(MP)")
			line=line.replace("seg 8","(NS)")
		#print(line)
		for segment in segments:
			if segment in line:
				if segment in genes.keys():
					genes[segment]+=line	
				else:
					genes[segment]=line
				break

#print (genes["(NEP)"])
#print ("_-_-_-_-_")
res=""
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
