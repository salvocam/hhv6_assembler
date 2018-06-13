import biomodule
from Bio import SeqIO
import os
import sys


def fuseSequences(s1,s2):
    a = -5
    while True:
        pos = s2.find(s1[a:])
        if not pos == -1:
            a = a -1
        else:
            break
        if abs(a) == len(s1):
            break
    fusedSequence = s1[:a+1]+s2[s2.find(s1[a+1:]):]
    return fusedSequence


sequences = {}
sequence2Assemble = open("toAssemble.fasta","w")

#Collects spades contigs and elonged repeats
for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    locus = str(seq_record.id)
    if not locus in sequences:
        sequences[locus] = str(seq_record.seq)



infile = open(sys.argv[2])

line = infile.readline().rstrip()
while not "#contig" in line:
    line = infile.readline().rstrip()
    print line

line = infile.readline().rstrip()
coord = []
field = line.split("\t")
for item in field:
    if not item == '':
        coord.append(item)

locus1 = coord[0]
refCoord_start1 = int(coord[1])
refCoord_end1 = int(coord[2])
coord1_start = int(coord[3])
coord1_end = int(coord[4])
if coord1_start > coord1_end:
    startingScaffold = biomodule.reverseComplement(sequences[locus1])
else:
    startingScaffold = sequences[locus1]


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
    refCoord_start2 = int(coord[1])
    refCoord_end2 = int(coord[2])
    coord2_start = int(coord[3])
    coord2_end = int(coord[4])



    if coord2_start > coord2_end:
        s2 = biomodule.reverseComplement(sequences[locus2])
    else:
        s2 = sequences[locus2]
    
    if refCoord_start2 > refCoord_end1:
        print "Entra qui!",refCoord_start2,refCoord_end1
        for a in range(int(refCoord_start2)-int(refCoord_end1)+1):
            startingScaffold += "N"
        startingScaffold += s2
    else:
        startingScaffold = fuseSequences(startingScaffold,s2)


    halfSteps.write(">"+str(step)+"\n"+startingScaffold+"\n")

    step += 1
    coord1_start = coord2_start
    coord1_end = coord2_end
    refCoord_start1 = refCoord_start2
    refCoord_end1 = refCoord_end2


finalScaffold = open("finalScaffold.fasta","w")
finalScaffold.write(">FinalScaffold\n"+startingScaffold)





