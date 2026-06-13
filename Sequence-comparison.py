#-------------PROJECT's STEPS-----------
# download 2 genomes from NCBI
# compare the 2 organisms
# calculate hamming dist (point mutation )
# global alignment for both sequences
# find the similarity percentage 
# ---------------------------------------

from Bio.Seq import Seq                   # Import Seq object for sequence manipulation
from Bio import SeqIO                     # Import SeqIO for reading sequence files
from Bio import Align                     # Import Align for sequence alignment

# Read the first organism sequence (Homo sapiens COX1)
seq1=SeqIO.read("/home/useramna/ROSALIND-Examples.py/cox1 homo.fasta", "fasta")
homo=seq1.seq                                # Store the human COX1 sequence

# Read the second organism sequence (Mus musculus COX1)
seq2=SeqIO.read("/home/useramna/ROSALIND-Examples.py/cox1 mus.fasta", "fasta")
mus=seq2.seq                                 # Store the mouse COX1 sequence

count=0                            # Initialize counter for point mutations
for a, b in zip(homo, mus):            # Loop through both sequences simultaneously
    if a!=b:                             # Check if nucleotides are different
        count+=1                           # Increment counter for each difference

# Print the total number of point mutations (Hamming Distance)
print(f"Point Mutation (Hamming Distance): {count}")


aligner=Align.PairwiseAligner()                  # Create a pairwise sequence aligner object

aligner.mode="global"                           # Set alignment mode to global (Needleman-Wunsch)
aligner.match_score= +2                         # Score for matching nucleotides
aligner.mismatch_score= -1                      # Penalty for mismatching nucleotides
aligner.open_gap_score=-1                       # Penalty for opening a new gap

alignments= aligner.align(homo,mus)                  # Perform the alignment between both sequences
best=alignments[0]                                   # Get the best alignment (highest score)
print(best)                               # Print the best alignment
print(f"Score: {best.score}")             # Print the alignment score


# Calculate similarity percentage
# max(0,...) ensures the percentage never goes below 0
similarity_score= max(0, best.score / max(len(homo), len(mus)) * 100)
print(f"Similarity: {similarity_score:.2f} %")                # Print similarity percentage


