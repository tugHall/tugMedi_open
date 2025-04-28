

Command Line Arguments
python connect_ccds.py <ccds_input_file> <ccds_output_file> <gene_file>
- ccds_input_file:  Name of an input file, CCDS data file downloaded from the National Center for Biotechnology Information.
- ccds_output_file: Name of the output file, cconnected CCDS file for the simulator.
- gene_file:        Name of an input file, IDs of tumor-related genes of your interest. An example is in Data/. 

What the Script Does
- Load the CCDS data from the ccds input file and keep only entries with the "Public" status.
- Connect all genes in each chromosome. All genes in a chromosome is written like a (artificial) single gene in the output. 
- Extract genes specified in the gene file from the CCDS data, and write them in the output. 


