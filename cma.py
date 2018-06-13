import sys
import os
import datetime
from Bio import SeqIO
import biomodule as bm

conf = open(sys.argv[1])


#Read config file
sampleName = ((conf.readline().rstrip()).split("\t"))[1]
comment = conf.readline()
read1 = ((conf.readline().rstrip()).split("\t"))[1]
read2 = ((conf.readline().rstrip()).split("\t"))[1]
outputFolder = ((conf.readline().rstrip()).split("\t"))[1]
comment = conf.readline()
performConsensusMaker = ((conf.readline().rstrip()).split("\t"))[1]
consensusFile = ((conf.readline().rstrip()).split("\t"))[1]
performDeNovoAssembly = ((conf.readline().rstrip()).split("\t"))[1]
denovoAssemblyScaffolds = ((conf.readline().rstrip()).split("\t"))[1]
performRepeatFlankingElongment = ((conf.readline().rstrip()).split("\t"))[1]
elongedRepetsFile = ((conf.readline().rstrip()).split("\t"))[1]
comment = conf.readline()
bwaEditDist = ((conf.readline().rstrip()).split("\t"))[1]
bwaNumThreads = ((conf.readline().rstrip()).split("\t"))[1]

logFileName = sampleName+".log"
logFile = open(logFileName,"w")

#Start pipeline

#Create consensus if required *******************************************************************************
if performConsensusMaker == "yes" or performConsensusMaker == "Yes":
    logFile.write("1 Creation of a consensus sequence from the merlin reference\n")

    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
    os.system("./1_consensusRetriever/merlinAlignment "+sampleName+" " + read1+" "+read2+" "+outputFolder+"/ "+bwaEditDist+" "+bwaNumThreads)

    now = datetime.datetime.now()
    logFile.write("Process ended at "+now.strftime("%Y-%m-%d %H:%M"))
    logFile.write("\n\n")

    #Load consensus sequence in memory
    for seq_record in SeqIO.parse(outputFolder+"/"+sampleName+"/hcmv_genome.fasta_con.fasta","fasta"):
        consensusSequence = str(seq_record.seq)
    
    

else:
    for seq_record in SeqIO.parse(consensusFile,"fasta"):
        consensusSequence = str(seq_record.seq)
    
#**************************************************************************************************************


#Create file with repeat flanking regions fro consensus sequence (this is needed by pipeline step 2)***********
flankingSequencesFile = open("repeatsFlanking.fasta","w")
flankingSequencesFile.write(">TRLflankingStarting\n"+bm.reverseComplement(consensusSequence[1364:3000])+"\n")
flankingSequencesFile.write(">TRLflankingEnding\nCCATTCCGGGCCGTGTGCTGGGTCCCCGAGGGGCGGGGGGGTGTTTTCTGCGGGGGGGTGAAATTTGGAGTTGCGTGTGTGGACGGCGACGGCGACTAGTTGCGTGTGCTGCGGTGGGTACGGCGACGGCGAATAAAAGCGACGTGCGGCGCGCACGGCGAAAAGCAGACGCGCGTCTGTGTCTGTTTGAGTCCCCAGGGGACGGCAGCG\n")
flankingSequencesFile.write(">IRflankingStarting\n"+consensusSequence[192000:193500]+"\n")
flankingSequencesFile.write(">IRflankingEnding\n"+consensusSequence[197000:197500]+"\n")
flankingSequencesFile.write(">TRSflankingEnding\n"+consensusSequence[232000:233500]+"\n")
flankingSequencesFile.write(">TRLflankingEnding\nCCCGGCCAACACACCCCGACACACCCGGCACACGCCCGCGACACACCCGGCCAACACACCCCGACACACCCGGCACACGCCCGCGACACACCCGCGGCACACCCTGACACACCCGCCACACCCGGCACACACCCACCCCGCCGCGCCCCCGACACACCCCGACCGCCGCCGGTGCGGGACAGGGCT\n")
flankingSequencesFile.close()
os.system("mv repeatsFlanking.fasta ./2_ElongationFlankingRepeats/")
#****************************************************************************************************************

#Launch the repeats flanking sequence elongation pipeline ********************************************************
#Elong TRL flanking
if performRepeatFlankingElongment == "yes" or performRepeatFlankingElongment == "Yes":
    elongedSeqencesFile = open("elongedSequences.fasta","w")

    os.chdir("./2_ElongationFlankingRepeats/pipeline/")

    logFile.write("Trying to elong the TRL flanking region\n")
    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")

    startingSeqFile = open("starting.fasta","w")
    startingSeqFile.write(">consensuSequence [1364:3000]\n"+bm.reverseComplement(consensusSequence[1364:3000])+"\n")
    startingSeqFile.close()
    os.system("mv starting.fasta ./startingScaffolds/")

    endingSeqFile = open("termini.fasta","w")
    endingSeqFile.write(">termini\nCCATTCCGGGCCGTGTGCTGGGTCCCCGAGGGGCGGGGGGGTGTTTTCTGCGGGGGGGTGAAATTTGGAGTTGCGTGTGTGGACGGCGACGGCGACTAGTTGCGTGTGCTGCGGTGGGTACGGCGACGGCGAATAAAAGCGACGTGCGGCGCGCACGGCGAAAAGCAGACGCGCGTCTGTGTCTGTTTGAGTCCCCAGGGGACGGCAGCG\n")
    endingSeqFile.close()

    os.system("python elong.py "+sampleName+" "+read1+" "+read2)

    for seq_record in SeqIO.parse("elonged.fasta","fasta"):
        elongedSequence = str(seq_record.seq)
    elongedSeqencesFile.write(">TRL_elonged\n"+bm.reverseComplement(elongedSequence)+"\n")

    now = datetime.datetime.now()
    logFile.write("Process ended at "+now.strftime("%Y-%m-%d %H:%M")+"\n")

    #Elong IR block flanking
    logFile.write("Trying to elong the IR block flanking region\n")
    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")

    startingSeqFile = open("starting.fasta","w")
    startingSeqFile.write(">IRflankingStarting\n"+consensusSequence[192000:193500]+"\n")
    startingSeqFile.close()
    os.system("mv starting.fasta ./startingScaffolds/")

    endingSeqFile = open("termini.fasta","w")
    endingSeqFile.write(">termini\n"+consensusSequence[197000:197500]+"\n")
    endingSeqFile.close()

    os.system("python elong.py "+sampleName+" "+read1+" "+read2)

    for seq_record in SeqIO.parse("elonged.fasta","fasta"):
        elongedSequence = str(seq_record.seq)
    elongedSeqencesFile.write(">IR_elonged\n"+elongedSequence+"\n")

    now = datetime.datetime.now()
    logFile.write("Process ended at "+now.strftime("%Y-%m-%d %H:%M")+"\n")


    #Elong TRS block flanking
    logFile.write("Trying to elong the TRS block flanking region\n")
    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")

    startingSeqFile = open("starting.fasta","w")
    startingSeqFile.write(">TRSflankingStarting\n"+consensusSequence[232000:233500]+"\n")
    startingSeqFile.close()
    os.system("mv starting.fasta ./startingScaffolds/")

    endingSeqFile = open("termini.fasta","w")
    endingSeqFile.write(">termini\nCCCGGCCAACACACCCCGACACACCCGGCACACGCCCGCGACACACCCGGCCAACACACCCCGACACACCCGGCACACGCCCGCGACACACCCGCGGCACACCCTGACACACCCGCCACACCCGGCACACACCCACCCCGCCGCGCCCCCGACACACCCCGACCGCCGCCGGTGCGGGACAGGGCT\n")
    endingSeqFile.close()

    os.system("python elong.py "+sampleName+" "+read1+" "+read2)

    for seq_record in SeqIO.parse("elonged.fasta","fasta"):
        elongedSequence = str(seq_record.seq)
    elongedSeqencesFile.write(">TRS_elonged\n"+elongedSequence+"\n")

    now = datetime.datetime.now()
    logFile.write("Process ended at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
    elongedSeqencesFile.close()
    os.chdir("../../")
    os.system("cp elongedSequences.fasta ./4_ScaffoldBuilder/")
else:
    os.system("cp " +elongedRepetsFile+ " ./4_ScaffoldBuilder/")



#*******************************************************************************************************************************************


#Perform the denovo assembly step ********************************************************
if performDeNovoAssembly == "yes" or performDeNovoAssembly == "Yes":
    logFile.write("Trying to elong the Denovo assembly step \n")
    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
    os.system("~/Software/SPAdes-3.12.0-Linux/bin/spades.py -1 "+read1+" -2 "+read2+" -o ./3_DenovoAssembly/")
    now = datetime.datetime.now()
    logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
    os.system("python ./Software/representative.py ./3_DenovoAssembly/scaffolds.fasta")
    os.system("cp representativeScaffolds.fasta ./4_ScaffoldBuilder/")
else:
    os.system("python ./Software/representative.py "+denovoAssemblyScaffolds)
    os.system("mv representativeScaffolds.fasta ./4_ScaffoldBuilder/")
#***************************************************************************************************************************


#Perform Scaffold buuilder step ********************************************************
logFile.write("Performing the scaffold builder step \n")
now = datetime.datetime.now()
logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
os.system("cat ./4_ScaffoldBuilder/elongedSequences.fasta ./4_ScaffoldBuilder/representativeScaffolds.fasta >./4_ScaffoldBuilder/sequencesToOrient.fasta")
os.system("scaffold_builder_v2.py -q ./4_ScaffoldBuilder/sequencesToOrient.fasta -r ./Software/hcmv_genome.fasta -p sbSB")
os.system("mv sbSB* ./4_ScaffoldBuilder/")
now = datetime.datetime.now()
logFile.write("Process Ended at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
#***************************************************************************************************************************

#Perform Scaffold merging step ********************************************************
logFile.write("Performing the scaffold merging step \n")
now = datetime.datetime.now()
logFile.write("Process started at "+now.strftime("%Y-%m-%d %H:%M")+"\n")
os.system("cp ./4_ScaffoldBuilder/sbSB_output.txt ./5_MergeScaffolds/")
os.system("cp ./4_ScaffoldBuilder/sequencesToOrient.fasta ./5_MergeScaffolds/")
os.system("python ./5_MergeScaffolds/mergeSB.py ./5_MergeScaffolds/sequencesToOrient.fasta  ./5_MergeScaffolds/sbSB_output.txt")
os.system("mv ./5_MergeScaffolds/finalScaffold.fasta ./")
os.system("rm -f tempAssem* toAssem*")











logFile.close()






















