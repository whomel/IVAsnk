import os,sys

#Filter vcf file. Discard all the entries where ALT number of reads in DP4 is lower than the REF.
#It is the same as `bcftools view -e '(DP4[2]+DP4[3])<(DP4[0]+DP4[1])' {output.vcf}.gz -o {output.vcfilt}`
#But, moreover, it checks if the ALT freq is > 20% AND minimum depth of 200x. In that case, it keeps the variation.

#python filterVCF.py vcf.file freq > filteredVCF.file


def minor(s,freq, line):
	ref1=int(s.split(",")[0])
	ref2=int(s.split(",")[1])
	alt1=int(s.split(",")[2])
	alt2=int(s.split(",")[3])
        #if ALT > 20% and depth>200x and it's NOT an indel, it keeps the mutation
	if (alt1+alt2)>(freq*(alt1+alt2+ref1+ref2)/100) and (alt1+alt2+ref1+ref2)>200 and 'INDEL' not in line:
		return False
	#if it's an INDEL, keep it only if it's MAJOR.
	elif (alt1+alt2)>(freq*(alt1+alt2+ref1+ref2)/100) and (alt1+alt2+ref1+ref2)>200 and 'INDEL' in line:
		if (alt1+alt2)>(ref1+ref2):
                	return False
		else:
			return True
        #if depth < 200x, take major (so discard the mutation if ALT<REF, if not bcftools take the ALT, even if it's not the major)
	elif (alt1+alt2+ref1+ref2)<200 and (alt1+alt2)<(ref1+ref2):
		return True
	elif (alt1+alt2+ref1+ref2)<200 and (alt1+alt2)>(ref1+ref2):
		return False
	else:
		return True

res=""
with open (sys.argv[1],'r') as fi:
	for line in fi:
		if line[0]!="#":
			#DP4=38,18,144,33;
			dp4=line.split("DP4=")[1].split(";")[0]
			if not minor(dp4, int(sys.argv[2]), line):
				res+=line
		else:
			if line!="\n" and line!="":
				res+=line.replace('\n\n','\n')

print(res)
