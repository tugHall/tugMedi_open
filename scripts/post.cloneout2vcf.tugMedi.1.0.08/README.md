# cloneout2vcf.py Manual

# ---

`cloneout2vcf.py` is a script for converting the output of tugHall into VCF format. This script can operate in the following three modes.

## Modes

1.  **Both Mode**: Performs both weight calculation and VCF conversion.
2.  **Calcweight Mode**: Calculates mutation transition probabilities (weights) based on TCGA data.
3.  **Convert Mode**: Converts tugHall output to VCF format using pre-calculated weight data.

## Required Libraries

-   `pyfaidx`
-   `pandas`

## 1. Both Mode

### Overview

-   The `Both` mode calculates mutation transition probabilities based on TCGA data and uses them to convert tugHall output to VCF format.

### Input Data Specifications and Examples

-   **TCGA MAF File**: Input TCGA MAF file used for calculating mutation weights.
    -   Example: `TCGA-AZ-.....maf`
-   **tugHall Output Files**:
    -   `cloneout.txt`: Contains clone information for mutations.
    -   `pointMutations.txt`: Contains positional information for mutations.
    -   `VAF.txt`: Contains allele frequency information for mutations.
-   **FASTA File**: Reference genome sequence data.
    -   Example: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`

### Output Data Examples

-   **VCF File**:
    -   Example: `output_all_VAF_primary.vcf`
-   **Weight Data File**: File containing the calculated mutation weights.
    -   Example: `TEST/TCGA-AZ-....._alt_weight.txt`

### Execution Command and Option Descriptions

```bash
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

| Option                | Description                                                                                                |
| :-------------------- | :--------------------------------------------------------------------------------------------------------- |
| `--mode both`         | Specifies the 'both' mode (calculate weights and convert to VCF). (Default)                                |
| `--cloneout`          | Path to the `cloneout.txt` file (tugHall output).                                                          |
| `--pointMutations`    | Path to the `pointMutations.txt` file (tugHall output).                                                    |
| `--VAF`               | Path to the `VAF.txt` file (tugHall output).                                                               |
| `--fasta`             | Path to the reference genome FASTA file.                                                                   |
| `--maf`               | Path to the input TCGA MAF file.                                                                           |
| `--output`            | Directory where output files will be saved. (Default: `TEST`)                                              |
| `--time`              | Filters data based on the 'Time'. Multiple values/conditions allowed. Example: `--time 100 200 300`.       |
| `--time_filter_logic` | Logic ('and' or 'or') used to combine multiple `--time` filters. (Default: 'or')                           |
| `--VAF_column`        | Specifies which VAF column to use ('VAF_primary' or 'VAF_metastatic'). (Default: 'VAF_primary')            |
| `--seed`              | Use fixed random seed for reproducibility. (Default: False)                                                |

## 2. Calcweight Mode

### Overview

-   The `Calcweight` mode calculates mutation transition probabilities based on TCGA data and saves them to a weight file.

### Input Data Specifications and Examples

-   **TCGA MAF File**: Input TCGA MAF file used to calculate mutation transition probabilities.
    -   Example: `TCGA-AZ-.....maf`

### Output Data Examples

-   **Weight Data File**: File containing the calculated mutation weights.
    -   Example: `TEST/TCGA-AZ-....._alt_weight.txt`

### Execution Command and Option Descriptions

```bash
./cloneout2vcf.py \
--mode calcweight \
--maf TCGA-AZ-.....maf \
--output TEST
```

| Option             | Description                                                       |
| :----------------- | :---------------------------------------------------------------- |
| `--mode calcweight` | Specifies the 'calcweight' mode (only calculate weights).         |
| `--maf`            | Path to the input TCGA MAF file.                                  |
| `--output`         | Directory where the output weight file will be saved. (Default: `TEST`) |

## 3. Convert Mode

### Overview

-   The `Convert` mode converts tugHall output to VCF format using a pre-calculated weight file.

### Input Data Specifications and Examples

-   **tugHall Output Files**:
    -   `cloneout.txt`: Contains clone information for mutations.
    -   `pointMutations.txt`: Contains positional information for mutations.
    -   `VAF.txt`: Contains allele frequency information for mutations.
-   **FASTA File**: Reference genome sequence data.
    -   Example: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
-   **Weight Data File**: Path to the pre-calculated mutation weight file (generated by `calcweight` mode).
    -   Example: `TCGA-AZ-....._alt_weight.txt`

### Output Data Examples

-   **VCF File**:
    -   Example: `output_all_VAF_primary.vcf`

### Execution Command and Option Descriptions

```bash
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

| Option                | Description                                                                                                |
| :-------------------- | :--------------------------------------------------------------------------------------------------------- |
| `--mode convert`      | Specifies the 'convert' mode (convert to VCF using existing weights).                                        |
| `--cloneout`          | Path to the `cloneout.txt` file (tugHall output).                                                          |
| `--pointMutations`    | Path to the `pointMutations.txt` file (tugHall output).                                                    |
| `--VAF`               | Path to the `VAF.txt` file (tugHall output).                                                               |
| `--fasta`             | Path to the reference genome FASTA file.                                                                   |
| `--alt_weight`        | Path to the pre-calculated mutation weight file.                                                           |
| `--output`            | Directory where the output VCF file will be saved. (Default: `TEST`)                                       |
| `--time`              | Filters data based on the 'Time'. Multiple values/conditions allowed. Example: `--time 100 200 300`.       |
| `--time_filter_logic` | Logic ('and' or 'or') used to combine multiple `--time` filters. (Default: 'or')                           |
| `--VAF_column`        | Specifies which VAF column to use ('VAF_primary' or 'VAF_metastatic'). (Default: 'VAF_primary')            |
| `--seed`              | Use fixed random seed for reproducibility. (Default: False)                                                |

## 4. Execution Examples

### Both Mode

```bash
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

```bash
echo "[INFO] start test calcweight mode"
./cloneout2vcf.py \
    --mode calcweight \
    --maf TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --output TEST_calcweight
```

### Convert Mode

```bash
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

## 5. Displaying Help

```bash
./cloneout2vcf.py -h
```

The script will display the following help message:

```
usage: cloneout2vcf.py [-h] [--time [TIME ...]] [--time_filter_logic TIME_FILTER_LOGIC] [--cloneout CLONEOUT] [--pointMutations POINTMUTATIONS] [--VAF VAF]
                       [--VAF_column VAF_COLUMN] [--maf MAF] [--alt_weight ALT_WEIGHT] [--fasta FASTA] [--mode MODE] [--output OUTPUT]

Convert three tugHall outputs to vcf

options:
  -h, --help            show this help message and exit
  --time [TIME ...]     Time corresponds to Time column in VAF_data. multiple values are allowed. if not specified, all times are used. also specify max/min as string.
                        Example: --time 1 2 3
                                  --time max
                                  --time min
                                  --time >0
                                  --time <100
  --time_filter_logic TIME_FILTER_LOGIC
                        Logic to combine time filters. 'and' or 'or'. [default=or]
  --cloneout CLONEOUT   First tugHall output data.
  --pointMutations POINTMUTATIONS
                        Second tugHall output data.
  --VAF VAF             Third tugHall output data.
  --VAF_column VAF_COLUMN
                        VAF_primary or VAF_metastatic. [default=VAF_primary]
  --maf MAF             TCGA maf filepath.
  --alt_weight ALT_WEIGHT
                        Pre-calculated alteration weights data.
  --fasta FASTA         genome fasta filepath.
  --mode MODE           Mode to run the program. [convert/calcweight/both, default=both]
  --output OUTPUT       Output directory to save output vcf file. [default=TEST]
  --seed                Use fixed random seed for reproducibility.

Usage Examples:
1. Run the script in both mode:
#both mode is default.
#please specify maf file and tugHall output files
cloneout2vcf.py
    --cloneout ~/cloneout.txt
    --pointMutations ~/pointMutations.txt
    --VAF ~/VAF.txt
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
    --maf TCGA-AZ-.....maf
    --output TEST

2. Run the script in calcweight mode:
#please specify maf file
cloneout2vcf.py
    --mode calcweight
    --maf TCGA-AZ-.....maf

3. Run the script in convert mode:
#Please specify the tugHall output files and the alt_weight intermediate file, which is produced in calcweight mode:
cloneout2vcf.py
    --mode convert
    --cloneout ~/cloneout.txt
    --pointMutations ~/pointMutations.txt
    --VAF ~/VAF.txt
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
    --alt_weight TCGA-AZ........txt
    --output TEST
```
