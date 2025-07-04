# Manual on cloneout2vcf.py

## Purpose

`cloneout2vcf.py` is a script that converts a cloneout.txt of `tugMedi` into the VCF format, using mutation transition probabilities given, or calculated from a MAF file.

## Overview

This script interacts with the following files and directories:

-   README.md - This file.
-   cloneout2vcf.py - The main script.
-   testdata/ - Contains an example of tugMedi output files.
-   TEST_Both/ - Contains output files generated when cloneout2vcf.py is run (in the Both mode).

## Usage

This section demonstrates how to use the script with data in the ./testdata/ directory.

### 1. Extract tugMedi output files

``` bash
cd testdata/
unzip 0003.zip
```

### 2. Save Homo_sapiens.GRCh38.dna.primary_assembly.fa to the testdata directory

``` bash
$ wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ./testdata/
```

### 3. Run cloneout2vcf.py

``` bash
./cloneout2vcf.py \
    --time max \
    --cloneout ./testdata/0003/Output/cloneout.txt \
    --pointMutations ./testdata/0003/Output/Mutations/pointMutations.txt \
    --VAF ./testdata/0003/Output/VAF/VAF.txt \
    --fasta ./testdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --maf ./testdata/TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --seed \
    --output TEST_both
```

### 4. Output Files

-   TEST_both/output_max_VAF_primary.vcf
    -   The converted VCF file
-   TEST_both/testdata/TCGA-AZ-6608-01A-11D-1835-10_may_alt_weight.txt
    -   Calculated transition probabilities

### Required Libraries

-   `pyfaidx`
-   `pandas`

## Options

### Modes

This script can operate in three modes: the Both, Calcweight, and Convert modes, which are explained below.

#### a. Both Mode

The Both mode 1) calculates mutation transition probabilities (weights) based on a MAF file and then 2) uses them to convert a coneout.txt into the VCF format.

##### Input Data Specifications and Examples

-   **TCGA MAF File**: Input TCGA MAF file used for calculating mutation weights.
    -   Example: `TCGA-AZ-.....maf`
-   **tugMedi Output Files**:
    -   `cloneout.txt`: Contains clone information for mutations.
    -   `pointMutations.txt`: Contains positional information for mutations.
    -   `VAF.txt`: Contains allele frequency information for mutations.
-   **FASTA File**: Reference genome sequence data.
    -   Example: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`

##### Output Data Examples

-   **VCF File**:
    -   Example: `output_all_VAF_primary.vcf`
-   **Weight Data File**: File containing the calculated mutation weights.
    -   Example: `TEST/TCGA-AZ-....._alt_weight.txt`

##### Execution Command and Option Descriptions

``` bash
./cloneout2vcf.py \
--mode both \
--cloneout ~/cloneout.txt \
--pointMutations ~/pointMutations.txt \
--VAF ~/VAF.txt \
--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--maf TCGA-AZ-.....maf \
--output TEST
--seed
```

| Option | Description |
|:---|:---|
| `--mode both` | Specifies the 'both' mode. (Default) |
| `--cloneout` | Path to the `cloneout.txt` file. |
| `--pointMutations` | Path to the `pointMutations.txt` file. |
| `--VAF` | Path to the `VAF.txt` file. |
| `--fasta` | Path to the reference genome FASTA file. |
| `--maf` | Path to the input TCGA MAF file. |
| `--output` | Directory where output files will be saved. (Default: `TEST`) |
| `--time` | Filters data based on the 'Time'. Multiple values/conditions allowed. Example: `--time 100 200 300`. |
| `--time_filter_logic` | Logic ('and' or 'or') used to combine multiple `--time` filters. (Default: 'or') |
| `--VAF_column` | Specifies which VAF column to use ('VAF_primary' or 'VAF_metastatic'). (Default: 'VAF_primary') |
| `--seed` | To fix the seed value, use the --seed option. |

#### b. Calcweight Mode

The Calcweight mode calculates mutation transition probabilities based on a MAF file.

##### Input Data Specifications and Examples

-   **TCGA MAF File**: Input TCGA MAF file used to calculate mutation transition probabilities.
    -   Example: `TCGA-AZ-.....maf`

##### Output Data Examples

-   **Weight Data File**: File containing the calculated mutation weights.
    -   Example: `TEST/TCGA-AZ-....._alt_weight.txt`

##### Execution Command and Option Descriptions

``` bash
./cloneout2vcf.py \
--mode calcweight \
--maf TCGA-AZ-.....maf \
--output TEST
```

#### c. Convert Mode

The Convert mode converts a cloneout.txt into the VCF format, using a pre-calculated weight file.

##### Input Data Specifications and Examples

-   **tugMedi Output Files**:
    -   `cloneout.txt`: Contains clone information for mutations.
    -   `pointMutations.txt`: Contains positional information for mutations.
    -   `VAF.txt`: Contains allele frequency information for mutations.
-   **FASTA File**: Reference genome sequence data.
    -   Example: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
-   **Weight Data File**: Path to the pre-calculated mutation weight file (generated by `calcweight` mode).
    -   Example: `TCGA-AZ-....._alt_weight.txt`

##### Output Data Examples

-   **VCF File**:
    -   Example: `output_all_VAF_primary.vcf`

##### Execution Command and Option Descriptions

``` bash
./cloneout2vcf.py \
--mode convert \
--cloneout ~/cloneout.txt \
--pointMutations ~/pointMutations.txt \
--VAF ~/VAF.txt \
--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--alt_weight TCGA-AZ-....._alt_weight.txt \
--output TEST
--seed
```

## Execution Examples

### Both Mode

``` bash
# Run without specifying time (uses all time points)
echo "[INFO] start test both mode (all time points)"
./cloneout2vcf.py \
    --cloneout ./generator/0003/Output/cloneout.txt \
    --pointMutations ./generator/0003/Output/Mutations/pointMutations.txt \
    --VAF ./generator/0003/Output/VAF/VAF.txt \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --maf TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --output TEST_both

# Filter for the maximum time point only
echo "[INFO] start test both mode (max time)"
./cloneout2vcf.py \
    --time max \
    --cloneout ./generator/0003/Output/cloneout.txt \
    --pointMutations ./generator/0003/Output/Mutations/pointMutations.txt \
    --VAF ./generator/0003/Output/VAF/VAF.txt \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --maf TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --output TEST_both_max_time

# Specify a specific time point
    --time 100 \

# Specify multiple time points (default logic is "or")
    --time 100 200 \

# Specify multiple conditions using comparison operators (requires quotes)
    --time ">100" "<=105" \
# Use 'and' logic to combine time filters
    --time_filter_logic "and" \
```

### Calcweight Mode

``` bash
echo "[INFO] start test calcweight mode"
./cloneout2vcf.py \
    --mode calcweight \
    --maf TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --output TEST_calcweight
```

### Convert Mode

``` bash
# Example using the weight file generated by calcweight mode
echo "[INFO] start test convert mode with calcweight intermediate file"
./cloneout2vcf.py \
    --mode convert \
    --time max \
    --cloneout ./generator/0003/Output/cloneout.txt \
    --pointMutations ./generator/0003/Output/Mutations/pointMutations.txt \
    --VAF ./generator/0003/Output/VAF/VAF.txt \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --alt_weight ./TEST_calcweight/TCGA-AZ-6608-01A-11D-1835-10_may_alt_weight.txt \
    --output TEST_convert
```
