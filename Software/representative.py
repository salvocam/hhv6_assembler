from Bio import SeqIO
import numpy as nm
import sys
sequences = {}
maxCov = 1
maxLength = 1



covDist = {}
maxLen = 0
for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    locus = str(seq_record.id)
    if not locus in sequences:
        sequences[locus] = str(seq_record.seq)
        coverage = (float( (locus.split("_cov_"))[1] ))
        length = (int( (((locus.split("_length_"))[1]).split("_"))[0] ))
        if length > maxLen:
            maxLen = length
            maxCov = coverage
        if not coverage in covDist:
            covDist[coverage] = length
        if coverage > maxCov and length>3000:
            maxCov = coverage


orderedCov = []
for key, value in sorted(covDist.items()): # Note the () after items!
    orderedCov.append((key, value))

print "Max coverage",maxCov
print "Max length",maxLen

#sys.stdin.read(1)
totalBases = 0
represenativeScaffolds = open("representativeScaffolds.fasta","w")
for a in range(len(orderedCov)-1,1,-1):
    if orderedCov[a][0] > maxCov/3 and orderedCov[a][1] >1000:
        totalBases += orderedCov[a][1]
        for locus in sequences:
            if str(orderedCov[a][0]) in locus:
                represenativeScaffolds.write(">"+locus+"\n"+sequences[locus]+"\n")


print "Total bases ",totalBases








