import biomodule
from Bio import SeqIO
import os
import sys

sequeunces = {}
sequence2Assemble = open("toAssemble.fasta","w")

#Collects spades contigs and elonged repeats
for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    locus = str(seq_record.id)
    if not locus in sequeunces:
        sequeunces[locus] = str(seq_record.seq)



infile = open(sys.argv[2])
for a in range(9):
    line = infile.readline().rstrip()
coord = []
field = line.split("\t")
for item in field:
    if not item == '':
        coord.append(item)

locus1 = coord[0]
coord1_start = int(coord[3])
coord1_end = int(coord[4])
startingScaffold = sequeunces[locus1]

step =1 
halfSteps = open("halpfSteps.fasta","w")
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    coord = []
    field = line.split("\t")
    for item in field:
        if not item == '':
            coord.append(item)
    locus2 = coord[0]
    coord2_start = int(coord[3])
    coord2_end = int(coord[4])

    
    
    tempAssembly = open("tempAssembly.fasta",'w')
    if coord2_start > coord2_end:
        tempAssembly.write(">PartialGenome"+"\n"+startingScaffold+"\n"+">"+locus2+"\n"+biomodule.reverseComplement(sequeunces[locus2])+"\n") 
    else:
        tempAssembly.write(">PartialGenome"+"\n"+startingScaffold+"\n"+">"+locus2+"\n"+sequeunces[locus2]+"\n")
    tempAssembly.close()

    os.system("phrap tempAssembly.fasta")

    numSeq = 0
    for seq_record in SeqIO.parse("tempAssembly.fasta.contigs","fasta"):
        startingScaffold = str(seq_record.seq)
        numSeq += 1
        halfSteps.write(">"+str(step)+"\n"+startingScaffold+"\n")
        if numSeq>1:
            print "More than one seq"
            exit()
    if numSeq == 0:
        print "no assembly"
        print locus1,locus2
        exit()
    
    step += 1
    coord1_start = coord2_start
    coord1_end = coord2_end




