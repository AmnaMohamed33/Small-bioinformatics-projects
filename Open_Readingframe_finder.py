#----------PROJECT's STEPS----------------
# download a genome from NCBI
# find ORF on strand and rev.comp
# arrange them from longest to shortest
# calculate the mass of each protein
# ----------------------------------

from Bio import SeqIO                       # Import SeqIO for reading sequence files
from Bio.Seq import Seq                     # Import Seq object for sequence manipulation
 
# Read a single sequence from the FASTA file
record=SeqIO.read("/home/useramna/ROSALIND-Examples.py/sequence.fasta", "fasta")     
seq=record.seq                       # Store the DNA sequence

def find_orf(seq):
   orf=[]                            # Initialize empty list to store all found ORFs
   for i in range(len(seq)-2):                # Loop through every position in the sequence
      if seq[i:i+3]=='ATG':                   # Check if current position is a start codon
         protein=""                        # Initialize empty string to build the protein
         for j in range(i, len(seq)-2, +3):              # Read codons 3 by 3
            codon=seq[j:j+3]
            if codon in ['TAA', 'TGA', 'TAG']:                # Check if current codon is a stop codon
                if protein:
                   orf.append(protein)                       # Save the protein to the list
                break      
            else:
                protein += str(Seq(codon).translate())               # Translate codon to amino acid
   return orf                                                        # Return list of all found ORFs

rev_comp=seq.reverse_complement()                # Calculate the reverse complement of the sequence

orf1= find_orf(str(seq))                         # Find ORFs in the original strand
orf2= find_orf(rev_comp)                         # Find ORFs in the reverse complement strand

total_orfs= set(orf1 + orf2)                     # Combine both lists and remove duplicates using set
total_orfs_sorted= sorted(total_orfs, key=len, reverse= True)       # Sort ORFs from longest to shortest

# Amino acid mass table in Daltons
table_mass={'A':71.03711,
'C':103.00919,
'D':115.02694,
'E':129.04259,
'F':147.06841,
'G':57.02146,
'H':137.05891,
'I':113.08406,
'K':128.09496,
'L':113.08406,
'M':131.04049,
'N':114.04293,
'P':97.05276,
'Q':128.05858,
'R':156.10111,
'S':87.03203,
'T':101.04768,
'V':99.06841,
'W':186.07931,
'Y':163.06333
}


# Loop through all ORFs and calculate the mass of each protein
for orf in total_orfs_sorted:
    mass= sum(table_mass[aa] for aa in orf )             # Sum the mass of each amino acid
    print(f"length = {len(orf):3d} | Mass = {mass:.3f}| {orf}")          # Print the length, mass, and sequence of each ORF




   


