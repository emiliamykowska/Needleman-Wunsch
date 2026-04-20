**Needleman–Wunsch Sequence Alignment Tool**

Python program for global alignment of DNA and protein sequences using the Needleman–Wunsch algorithm.

**Overview**

The program enables quantitative comparison of two biological sequences (DNA or protein). It supports multiple input methods and provides both numerical and graphical results of the alignment.
The application is interactive and guides the user step by step, making it easy to use even for someone without programming experience.

**Features**
**Sequence Input**

Sequences can be provided in three ways:

- manual input
- loading from a FASTA file
- fetching directly from the NCBI database using an accession ID

Input methods can be mixed (e.g., one sequence from FASTA, the other from NCBI).

**Validation**

The program validates all input data:

DNA sequences: only A, C, G, T allowed

Protein sequences: only A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y amino acids allowed

Invalid input results in an error message.

**Algorithm Parameters**

User can:
use default parameters:
- match = 1
- mismatch = -1
- gap = -2
or define custom values, which are then validated to ensure logical scoring.

**Alignment Output**

The program generates:

- one optimal global alignment
- alignment score
- alignment length
- percentage of:
    - identical positions
    - gaps

Alignment is displayed in a format:

ACTG

\*\*|*

ACTT

**Visualization**

Short sequences (< 100 characters)
- scoring matrix displayed as a heatmap
- optimal path marked in red

Long sequences (> 100 characters)
- simplified plot showing only the optimal path

All visualizations are saved as image files.

**Report generation**

Results are saved to a text file (needleman_report.txt) including:

- sequences used
- algorithm parameters
- alignment
- final score
- alignment length
- percentage of matches and gaps

**User Instructions**

The program is designed to be used from the command line and guides the user step by step.

- Make sure Python 3 is installed
- Run the following command in the command line:
  pip install -r requirements.txt
- Open a terminal and navigate to the folder containing the program
- Run the program with the command:
  python needleman.py
- Choose sequence type for both sequences:
  1. protein
  2. DNA
  
  For each sequence, choose input method:
  1. manual input
  2. FASTA file (provide file path)
  3. NCBI (provide accession ID)
  
  Choose scoring parameters:
  1. default values
  2. customized values
   


