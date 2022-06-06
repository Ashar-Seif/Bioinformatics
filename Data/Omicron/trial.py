from Bio import SeqIO
import numpy as np
reference_sequences = [refrence.seq for refrence in SeqIO.parse("Omicron.sequences.fasta", "fasta")]
reference_sequences=np.asarray(reference_sequences)
print(len(np.unique(reference_sequences)))



