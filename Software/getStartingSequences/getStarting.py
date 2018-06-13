import os
import sys

os.system("cp "+sys.argv[1]+"  ./consensus.fasta")
os.system("makeblastdb -in consensus.fasta -dbtype nucl")
os.system("blastn -query IRending.fasta -db consensus.fasta -outfmt 6 -out outputBlast.txt")

