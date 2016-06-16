# This script will take a fasta file with multiple sequences and blast them one by one against the nr database
# using the Bio.Blast tools from the biopython package.  The output that will be saved is:
#

import os
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = "chadwick@caltech.edu"  #Let Entrez know who you are
file_path = raw_input("Path:")
files = raw_input("File:")
os.chdir(file_path) #Set directory to one containing fasta
file_name = files.split('.')[0]
sequence_list = list(SeqIO.parse(files,"fasta")) #Reads in fasta file
file = open(file_name+"_relatives.txt", "w") #Makes a file to store ANME relatives
file.close()
file = open(file_name+"_relatives_unique.txt", "w") #Makes a file to store ANME relatives
file.close()
file = open(file_name+"_relatives.fasta", "w") #Makes a file to store ANME relatives
file.close()
file = open(file_name+"_report.txt", "w") #Makes a file to store report
file.close()
counter = 1
for sequence in sequence_list:
  counter1 = 0
  counter2 = 0
  print counter, 'of', len(sequence_list)
  counter+=1
  print sequence.format("fasta")
  sequence_id = sequence.format("fasta").split('     ')[3].split('\n')[0]
  split_holder = sequence.format("fasta").split('     ')[1].split(' ')[0]
  bp = int(split_holder)
  if bp >= 700:
    sequence_relatives = NCBIWWW.qblast("blastn", "nr", sequence.seq, hitlist_size = 10000, perc_ident = 90, megablast='TRUE')
    save_1_file = open("test_blast.xml", "w")
    save_1_file.write(sequence_relatives.read())
    save_1_file.close()
    sequence_relatives.close()
    sequence_relatives = open("test_blast.xml")
    blast_1_record = NCBIXML.read(sequence_relatives)
    file = open(file_name+"_relatives.txt", "r")
    current_sequence = file.read()
    file.close()
    save_1_file = open(file_name+"_relatives.txt", "a")
    save_1_file2 = open(file_name+"_relatives_unique.txt", "a")
    save_1_report = open(file_name+"_report.txt", "a")
    save_1_file.write(sequence_id)
    save_1_file.write('\n')
    save_1_fasta = open(file_name+"_relatives.fasta", "a")
    for alignment in blast_1_record.alignments:
      for hsp in alignment.hsps:
        length = len(hsp.query)
        print alignment.title
        print length
        if length >= 600:
          index = current_sequence.find(alignment.title)
          print index
          counter1 += 1
          save_1_file.write(alignment.title)
          save_1_file.write('\n')
          if index == -1:
            save_1_file2.write(alignment.title)
            save_1_file2.write('\n')
            save_1_fasta.write(sequence.format("fasta"))
          else:
            counter2 += 1
    save_1_report.write(sequence_id)
    save_1_report.write(" Total Relatives: ")
    save_1_report.write(str(counter1))
    save_1_report.write(" New Relatives: ")
    save_1_report.write(str(counter1-counter2))
    save_1_report.write('\n')
    save_1_report.close()
    save_1_file.close()
    save_1_file2.close()
  else:
    save_1_report = open(file_name+"_report.txt", "a")
    save_1_report.write(sequence_id)
    save_1_report.write(" too short: ")
    save_1_report.write(split_holder[0])
    save_1_report.write(" base pairs")
    save_1_report.write('\n')
    save_1_report.close()
print "done"
