from Bio import SeqIO
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd


covid_path="original data\COVID-19\COVID.sequences.fasta"
omicron_path="original data\Omicron\Omicron.sequences.fasta"
aligned_covid_path="cov-aligned-gisaid_hcov-19_2021_12_31_12.fas"
aligned_omicron_path="omicron_aligned_gisaid_hcov-19_2021_12_31_12.fas"
consvsomicron_path="consensusomicrongisaid_hcov-19_2021_12_31_12.fas"
extension='fasta'


# Read 10 sequence records from a file using its path and extension

def read_10seq(file_path,extension):

    records=[]
    for record in SeqIO.parse(file_path, extension):
        records.append(record.seq)
    return records

# Make an alignment profile for the reference sequences to know the count of each nuclotide at every position across all sequences and choose the dominant one  

def Alignment_profile(alignment):

    n= len(alignment[0])

    profile_matrix =  { 
        'A': [0]*n,
        'C': [0]*n,
        'G': [0]*n,
        'T': [0]*n,
        'N': [0]*n,
        '-': [0]*n
    }

    # loop acroos all sequences
    for record in alignment:
        # loop across all positions in one sequence
        for i in range (len(record)):
            # fill the profile matrix with each nucleotide count at every position in the sequence
            if record[i]=="A":
                profile_matrix['A'][i] += 1
            if record[i]=="C":
                profile_matrix['C'][i] += 1
            if record[i]=="G":
                profile_matrix['G'][i] += 1
            if record[i]=="T":
                profile_matrix['T'][i] += 1
            if record[i]=="N":
                profile_matrix['N'][i] += 1
            if record[i]=="-":
                profile_matrix['-'][i] += 1

    return n,profile_matrix

# choose the most dominant nuclotide from all sequences at every position and join thm together to a new sequence
def consensus(n,profile_matrix):

    sequence = []
    nucleotides=['A','C','G','T']
    # loop across all positions.
    for i in range(n):
        count=[]
        max_count=0
        max_nucleotide=''
        # loop across the 4 types of nucleotides and choose the dominant one
        for j in range(len(nucleotides)) :
            count.append(profile_matrix[nucleotides[j]][i])
            if count[j] > max_count:
              max_count = count[j]
              max_nucleotide = nucleotides[j]
        sequence.append(max_nucleotide)
    # join all dominant nucleotides together and generating the new consensus sequence
    consensus_sequence = "".join(sequence)
    return(consensus_sequence)

# calculate the A, C, G, T and CG contents of a sequence
def calculate_content(seq):

    length_of_sequence = len(seq)
    # count  A, C, G, T and CG of a sequence
    C_count = seq.count('C')
    G_count = seq.count('G')
    A_count = seq.count('A')
    T_count = seq.count('T')
    CG_count = C_count + G_count
    # calculate the percentage of each chemical constituent of a sequence
    C_content= (C_count/length_of_sequence)*100
    G_content= (G_count/length_of_sequence)*100
    A_content= (A_count/length_of_sequence)*100
    T_content= (T_count/length_of_sequence)*100
    CG_content= (CG_count/length_of_sequence)*100
    # putting all contents in a list
    contents = ["seq", C_content, G_content,A_content, T_content, CG_content]
    return(contents)


def dissimilar_region(seq):
    index_region=[]
    for case in seq:
        for index in range(0, len(seq[0])):
            if case[index] != seq[10][index]:
                index_region.append(index)
            else:
                continue
    # Delete the repition from the index_region 
    index_region.sort()

    return index_region

def main():
  # Read 10 sequences from each type after alignment
  covid=read_10seq(aligned_covid_path,extension)
  print(len(covid[0]))
  omicron=read_10seq(aligned_omicron_path,extension)
  print(len(omicron[0]))
  consensus_omicron=read_10seq(consvsomicron_path,extension)
  #Generate profile matrix
  n, profile_matrix=Alignment_profile(covid)
  #Generate consensus_sequence of the reference data
  consensus_sequence=consensus(n,profile_matrix)
  # Generate fasta file of the new consensus_sequence
  record = SeqRecord(Seq(consensus_sequence), id ="covidCons", description ="Consensus_sequence")
  SeqIO.write(record, "Consensus.fasta", "fasta")

  # calculate the percentage of all contents of the omicron sequences (case sequences) and their average
  omicron_contents=[]
  C_content_sum=0
  G_content_sum=0
  A_content_sum=0
  T_content_sum=0
  CG_content_sum=0
  
  for i in omicron:
      contents=calculate_content(i)
      C_content_sum+=contents[1]
      G_content_sum+=contents[2]
      A_content_sum+=contents[3]
      T_content_sum+=contents[4]
      CG_content_sum+=contents[5]
      omicron_contents.append(contents)
  
  avg_omicron_contents=["average",C_content_sum/10,G_content_sum/10,A_content_sum/10,T_content_sum/10,CG_content_sum/10]
  omicron_contents.append(avg_omicron_contents)


  # Display the results of Comparison between the chemical constituents of the consensus sequence and the case sequences
  col_names = [ "sequences","C_content", "G_content","A_content", " T_content", "CG_content"]
    
  print("The chemical constituents of the case sequences:")
  print(tabulate(omicron_contents, headers=col_names, tablefmt="fancy_grid", showindex="always"))
  
  consensus_content=calculate_content(consensus_sequence)
  print("Consensus Contents:",consensus_content)
  
  comparison = {"consensus_content": consensus_content, "Average_omicron_contents": avg_omicron_contents}
  print("Comparison between the chemical constituents of the consensus sequence and the case sequences:")
  print(tabulate([[k]+comparison[k] for k in comparison], headers=col_names, tablefmt="fancy_grid"))
  
  # Find a dissimilar region between case sequances and consensus sequance
  print(np.shape(consensus_omicron))
  dissimilar_region_index=list(dict.fromkeys(dissimilar_region(consensus_omicron)))
  print('The dissimilar region index  ',dissimilar_region_index) 
  print('Number of dissimilar regions  ',len(dissimilar_region_index)) 

 


  #Print the sequances that exist in  dissimilar_region_index
  for case in consensus_omicron:
      case_sequance_list = []
      consensus_sequance_list = []
      for index in dissimilar_region_index:
          case_sequance_list.append(case[index])
          consensus_sequance_list.append((consensus_omicron[10][index]))
      print(consensus_sequance_list)

  n_dissimilar=len(dissimilar_region_index)
  n_similar=len(consensus_omicron[0])- n_dissimilar
  print('Number of dissimilar regions', n_dissimilar)
  y = np.array([ n_dissimilar,  n_similar])
  mylabels = [ r"n_dissimilar:{}".format(n_dissimilar),r"n_similar:{}".format(n_similar)]
  plt.pie(y, labels = mylabels,colors=['#C1CDC1','#8B0A50'])
  plt.savefig("DissimilarityChart")
  plt.close()

  n_dissimilarP=np.round((n_dissimilar/len(consensus_omicron[0]))*100,2)
  n_similarP=100- n_dissimilarP
  X = np.array([ n_dissimilarP,  n_similarP])
  mylabels = [ r"Dissimilarity%:{}".format(n_dissimilarP),r"Similarity%:{}".format(n_similarP)]
  
  plt.pie(X, labels = mylabels,colors=['#C1CDC1','#8B3A62'])
  plt.savefig("DissimilarityPercentage") 
  
  
 
  

if __name__=="__main__":
    main()
       
        