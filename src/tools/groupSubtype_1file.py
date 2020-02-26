#Same as groupSubtype.py but instead a fasta folder, only one fasta.

import os,sys
import subprocess
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

#python groupSubtype.py /home/yi-mo/pipelines/miguel/paired/out_cutadapt_iseq/iSeqRun7_2019Avian2_miguelRun/consensus/post/

segN={"PB2":"1","PB1":"2","PA":"3","HA":"4","NP":"5","NA":"6","MP":"7","NS":"8"}

filename=sys.argv[1].replace("/rename","")
#for filename in glob.glob(sys.argv[1]+"/*.fa"):
with open (filename,'r') as fi:
    res=""
    fluType=""
    subtypeHA = ""
    subtypeNA = ""
    for line in fi:
        if "B VIRUS" in line:
            fluType="B"
        elif "A VIRUS" in line:
            #if "H" in line.split("))")[0].split("(")[-1]:
                #subtype = line.split("))")[0].split("(")[-1]
            fluType="A"
            lineclean=line.upper()
            lineclean=lineclean.replace(":",'').replace("SEGMENT 4","(HA)").replace("SEG 4","(HA)")
            lineclean=lineclean.replace("SEGMENT 6","(NA)").replace("SEG 6","(NA)")
            if "(HA)" in lineclean or "(NA)" in lineclean:
                if "(HA)" in lineclean:
                    if "H" in line.split("))")[0].split("(")[-1]:
                        subtypeHA = line.split("))")[0].split("(")[-1].split("N")[0]
                elif "(NA)" in lineclean:
                    if "H" in line.split("))")[0].split("(")[-1]:
                        subtypeNA = "N"+line.split("))")[0].split("(")[-1].split("N")[1]
    if fluType=="B":
        bname = os.path.basename(filename)
        #Create folder
        command = "mkdir -p "+ os.path.dirname(filename) +"/rename/type/FLUB/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #Move file
        command = "cp "+ os.path.dirname(filename) +"/rename/"+bname+" " + os.path.dirname(filename) +"/rename/type/FLUB/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #Rename "_HA" to ".4"
        for index,record in enumerate(SeqIO.parse(os.path.dirname(filename) +"/rename/type/FLUB/"+bname, 'fasta')):
            res+=">"+record.id.split("_")[0]+"."+segN[record.id.split("_")[1]]+"\n"+record.seq+"\n"
        with open (os.path.dirname(filename) +"/rename/type/FLUB/"+bname,"w") as fo:
            fo.write(str(res))
    elif fluType=="A":
        bname = os.path.basename(filename)
        #Create folder
        command = "mkdir -p "+ os.path.dirname(filename) +"/rename/type/FLUA/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #Create folder
        command = "mkdir "+ os.path.dirname(filename) +"/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #Move file
        command = "cp "+ os.path.dirname(filename) +"/rename/"+bname+" "+ os.path.dirname(filename) +"/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"
        process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #Rename "_HA" to ".4"
        for index,record in enumerate(SeqIO.parse(os.path.dirname(filename) +"/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"+bname, 'fasta')):
            res+=">"+record.id.split("_")[0]+"."+segN[record.id.split("_")[1]]+"\n"+record.seq+"\n"
        with open (os.path.dirname(filename) +"/rename/type/FLUA/"+subtypeHA+subtypeNA+"/"+bname,"w") as fo:
            fo.write(str(res))
