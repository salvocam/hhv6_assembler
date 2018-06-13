import os
import sys
from Bio import SeqIO
import biomodule

def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = []
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = []
                    longest = c
                    lcs_set.append(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.append(S[i-c+1:i+1])

    return lcs_set

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
        if reference[-15:] in sequence and (len(fuseSequences(reference,sequence))-len(reference)) > overhang:
            elongedSequence = fuseSequences(reference,sequence)
            overhang = len(elongedSequence)-len(reference)
            print "Sequence",elongedSequence
            print "Forward"
            print "Overhang",len(elongedSequence)-len(reference)
            #print "Dove si trova la sequenza ",sequence.find(reference[-15:])
            #print "Da cercare ",reference[-15:]
            print seq_record
            #sys.stdin.read(1)
        
        #Check reverse contigs
        revSequence = biomodule.reverseComplement(sequence)
        if reference[-15:] in revSequence and (len(fuseSequences(reference,revSequence))-len(reference)) > overhang:
            elongedSequence = fuseSequences(reference,revSequence)
            overhang = len(elongedSequence)-len(reference)
            print "Sequence",elongedSequence
            print "Reverse"
            print "Overhang",len(elongedSequence)-len(reference)
            #print "Dove si trova la sequenza ",sequence.find(reference[-15:])
            #print "Da cercare ",reference[-15:]
            print seq_record
            #sys.stdin.read(1)


    if overhang < 10:
        print "Poor elongment. Now exit...."
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
    





