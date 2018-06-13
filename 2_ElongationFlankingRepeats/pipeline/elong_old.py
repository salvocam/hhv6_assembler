import os
import sys
from Bio import SeqIO
import biomodule

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return lcs_set



os.system("cp ./startingScaffolds/starting.fasta ./toElong.fasta")

endReached = False
cycle = 1
while endReached == False:
    os.system("./launchAlignment "+sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3])
    #collect the toElong sequence:
    for seq_record in SeqIO.parse("toElong.fasta","fasta"):
        reference = str(seq_record.seq)

    #check the presence of the original dna portion in the assembled scaffolds
    bestScaffold = ""
    lengthBestScaffold = 0
    overhang = 0
    for seq_record in SeqIO.parse("unmapped.fasta.contigs","fasta"):
        sequence = str(seq_record.seq)
        print "Lunghezza migliore Scaffold:",lengthBestScaffold
        print "Overhang:",overhang
        #Check forward contigs

        if sequence[:20] in reference:
            if len(sequence) > lengthBestScaffold:
                lengthBestScaffold = len(sequence)
                startPosition = sequence.find(reference[:20])
                elongedSequence = sequence[startPosition:]
                bestScaffold = sequence
                overhang = len(sequence[startPosition:])-len(reference)
                print "Overhang composition:",len(sequence[startPosition:]),len(reference),str(seq_record.id)




        if reference[:20] in sequence:
            if len(sequence) > lengthBestScaffold:
                lengthBestScaffold = len(sequence)
                startPosition = sequence.find(reference[:20])
                elongedSequence = sequence[startPosition:]
                bestScaffold = sequence
                overhang = len(sequence[startPosition:])-len(reference)
                print "Overhang composition:",len(sequence[startPosition:]),len(reference),str(seq_record.id)

        #check reverse contigs
        if biomodule.reverseComplement(reference[:20]) in sequence:
            if len(sequence) > lengthBestScaffold:
                lengthBestScaffold = len(sequence)
                startPosition = (biomodule.reverseComplement(sequence)).find(reference[:20])
                elongedSequence = (biomodule.reverseComplement(sequence))[startPosition:]
                bestScaffold = sequence
                overhang = len((biomodule.reverseComplement(sequence))[startPosition:])-len(reference)
                print "Overhang composition:",len((biomodule.reverseComplement(sequence))[startPosition:]),len(reference),str(seq_record.id)

    if overhang < 10:
        print "Poor elongment. Now exit...."
        print "Selected contig:",sequence
        print "Searched string:",reference[:20]
        print "Found position:",startPosition 
        print "Overhang:",overhang
        os.system("cp toElong.fasta elonged.fasta")
        exit()
    else:
        tempFile=open("tempFile.fasta","w")
        tempFile.write(">toElong\n"+elongedSequence+"\n")
        tempFile.close()
        os.system("mv tempFile.fasta toElong.fasta")

    os.system("cp tempFile.fasta ./steps/tempFile.fasta_"+str(cycle)+".txt")
    os.system("cp toElong.fasta ./steps/toElong.fasta_"+str(cycle)+".txt")
    os.system("cp unmapped.fasta.contigs ./steps/unmapped.fasta.contigs_"+str(cycle)+".txt")


    #check if the elonged sequence reached the specified termini portion
    os.system("cat toElong.fasta termini.fasta > checkTermini.fasta")
    os.system("phrap -raw checkTermini.fasta")
    
    loci = []
    for seq_record in SeqIO.parse("checkTermini.fasta.singlets","fasta"):
        loci.append(str(seq_record.id))
    if not "termini" in loci:
        endReached = True
        os.system("mv toElong.fasta elonged.fasta")
    





