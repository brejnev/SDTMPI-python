SDTMPI_Linux64 (Sequence Demarcation Tool - MPI version for Linux 64 bits)
=========================================================================

SDTMPI_Linux64 is a free Linux-based program that runs on multiple cores to allow quick calculation of DNA sequence pairwise identities. 
The program is written in python and uses the mpi4py library to parallelize the process of pairwise alignment and similarity score 
calculation. Given a FASTA file containing DNA sequences, the program aligns all possible pairs of sequences using Muscle (Edgar, 2004), 
ClustalW2 (Larkin et al., 2007)  or Mafft (Katoh et al., 2009), calculates the sequence similarity score for each pair and 
uses a rooted neighbour joining phylogenetic tree to cluster closely related sequences based on similarity scores.
It outputs a text file containing the scores and a ".sdt" file that can be open with SDTv.1.0 Windows version to visualise the
pairwise identity distribution plot 
and colour coded similarity matrix.
The indentity scores are calculated as 1-(M/N) where M is the number of mismatching nucleotides and N is the total number of positions 
along the alignment where neither sequence has a gap character. 


DOWNLOAD AND INSTALLATION
==============================
1. Software requirements:

   The program requires the following to be installed:

   - Python2.7 
   - mpi4py, see installation instruction on http://mpi4py.scipy.org/docs/usrman/install.html
   - Mafft v6.923b (optional : only if required by the user), available at http://mafft.cbrc.jp/alignment/software/source.html
   - Neigbour (from the package phylib-3.69), available at http://evolution.genetics.washington.edu/phylip.html

   NB: muscle and clustalw executable files are located in the bin directory, in the main script "SDTMPI_Linux32.py" the path to these programs can 
       be changed to where these programs are installed on your system.

2. Install mpi4py and test whether it is working on your system.

3. Dowload the SDTMPI_Linux64 from http://web.cbio.uct.ac.za/SDT

4. Extract the SDTMPI_Linux64.tar.gz file into the location you want to place the program.
   The folder contains:
      - The bin directory which contains the executable files "muscle3.8.31_i86linux64", "clustalw2" and "neighbor".
      - The Bio directory which contains the Biopython library.
      - The output directory in which the output files after each run are stored.
      - SDTMPI_Linux64.py which the parallel version of the program.
      - A sample submission shell file (used to submit a job on the cluster).
      - test.fas a sample FASTA file to test the program.

5. Please change the mode of "muscle3.8.31_i86linux64", "clustalw2" and "neighbor" in the bin directory to "executable". 


6. In the SDMPI_Linux64.py script change the paths to Neighbour and alignment programs (MUSCLE,MAFFS and CLUSTALW) if you have them already installed on your system. Other wise use those provided in the SDTMPI_Linux64/bin directory and remember to enable their executable permission. 

7. Running commands:

	The script SDTMPI_Linux64.py takes two parameters, the path to the input fasta file and the name of the alignment program to be used.
	The execution command for the parallel version is as follows: 

	mpiexec -n 8 python SDTMPI_Linux64.py test.fas muscle  

	This will result in the use of MUSCLE as the alignment program and will run on 8 cores. Replace "muscle" by "clustal" or "mafft" (if mafft is installed) 
	to change the alignment program that is used. When the pairwise alignments and similarity score calculations are completed the scores will be written 
	into a text file named after the input FASTA file and save into the output folder and a ".sdt" file will be produced which can be open by the 
	SDTv.1.0 windows version to easily produce the plot, matrix and data. 

WARNING
-------
This code has been  used succesfully on Scientific Linux 5.4, 5.5 and 5.8, OpenMPI 1.4-4 and MPI4PY 1.3.
On other platforms for which the OpenMPI implementation does not support a call to the fork() (use of popen python command) the code will generate a error.

--------------------------------------------------------------------------------------------------------------------------------------------------------------

References 
----------
1. Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high
throughput, Nucleic Acids Research 32, 1792-1797. 
 
2. Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F,
Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ, Higgins DG. 
(2007). Clustal W and Clustal X version 2.0., Bioinformatics 23, 2947-2948 
 
3. Katoh K, Asimenos G, Toh H. (2009) Multiple Alignment of DNA Sequences with MAFFT,
Methods in Molecular Biology 537:39-64  
 
4. Felsenstein, J. (1995)PHYLIP (Phylogeny Inference Package) Version 3.57c, available at
http://www.med.nyu.edu/rcr/rcr/phylip/main.html#refs 

Authors  
-------

Brejnev Muhire [1]
Darren Martin [1]
Arvind Varsani [2]

[1] Institute of Infectious Diseases and Molecular Medicine (IIDMM)
    Computational Biology Group, 
    University of Cape Town
    South Africa 
[2] School of Biological Sciences 
    University of Canterbury
    Private Bag 4800 
    Christchurch, 8140
    New Zealand
BM is funded by the University of Cape Town
website: http://web.cbio.uct.ac.za/SDT               
email: mhrbre001@myuct.ac.za
email: mubrejnev@gmail.com	



