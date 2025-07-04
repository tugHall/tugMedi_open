# Tool for preparing observed data

**Tool Name:** prepare_observed_data.py

## Purpose

This tool formats a MAF or its derived file into files corresponding to samples.{Rint, Rother}.txt in Input/Samples/.

## Overview

This script interacts with the following files and directories:

-   README.md - This file.
-   src/prepare_observed_data.py - The main script.
-   src/{other_scripts} - Scripts imported from prepare_observed_data.py.
-   run_sample_Test1.sh - Script to run prepare_observed_data.py for sample_Test1.
-   run_sample_Test2.sh - Script to run prepare_observed_data.py for sample_Test2.
-   samples/
    -   sample_Test1: Folder containing config files and output results
    -   sample_Test2: Folder containing config files and output results

## Usage examples

### For Test1

#### 1. Command Execution

``` bash
python ./src/prepare_observed_data.py \
-i ../../inst/extdata/tugMedi.1.0/Test1/Input/Samples/Processes/samples.02.txt \
-c ./samples/sample_Test1/config/config.conf \
-o ./samples/sample_Test1/output_all
```

#### 2. Output File Description

-   samples/sample_Test1/output_all/samples.02_tumor_specific_vaf.Rint.txt
    -   Corresponds to samples.Rint.txt
-   samples/sample_Test1/output_all/samples.02_tumor_specific_vaf.Rother.txt
    -   Corresponds to samples.Rother.txt
-   samples/sample_Test1/output_all/samples.02_tumor_specific.txt
    -   Intermediate file: list of extracted tumor-specific mutations for samples.{Rint, Rother}.txt
    -   In this sample, since there is no comparison with normal samples, all mutations are extracted as tumor-specific mutations.
-   samples/sample_Test1/output_all/samples.02_tumor_specific_vaf.txt
    -   Intermediate file: list with calculated VAF values and vaf column added

#### 3. Execution by batch script

run_sample_Test1.sh contains commands over all steps. To run:

``` bash
bash run_sample_Test1.sh
```

### For Test2

#### 1. Command Execution

``` bash
python ./src/prepare_observed_data.py \
-i ../../inst/extdata/tugMedi.1.0/Test2/Input.TCGA_7903/Samples/Processes/TCGA.LUAD.55-7903.01.txt \
-c ./samples/sample_Test2/config/config.conf \
-o ./samples/sample_Test2/output_all
```

#### 2. Output File Description

-   samples/sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific_vaf.Rint.txt
    -   Corresponds to samples.Rint.txt
-   samples/sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific_vaf.Rother.txt
    -   Corresponds to samples.Rother.txt
-   samples/sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific.txt
    -   Intermediate file: list of extracted tumor-specific mutations for samples.{Rint, Rother}.txt.
    -   In this sample, since there is no comparison with normal samples, all mutations are extracted as tumor-specific mutations.
-   samples/sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific_vaf.txt
    -   Intermediate file: list with calculated VAF values and vaf column added

## Input data

### Mutation data file

-   A TSV file in the MAF format or its derived format
-   The following columns are required.

```         
  * Sample ID   
  * Gene name (Hugo_Symbol)  
  * Chromosome: Chromosome number (e.g., 1, X)  
  * Position: Position on the chromosome (e.g., 7674220)  
  * Ref: Reference allele (e.g., A)  
  * Alt: Alteration allele (e.g. G)
  * Ref_count: Count of reference allele
  * Alt_count: Count of alteration allele
  * VAF: Variant Allele Frequency
```

-   Note that columns actually required differ depending on calculation modes described below.

-   Typically, columns are selected from MAF.

Data example

```         
Hugo_Symbol    Chromosome  Start_Position  End_Position  Reference_Allele  Tumor_Seq_Allele1
TP53           17          7674220         7674220       G                 G
EGFR           7           55249071        55249071      A                 A                  
BRCA1          17          41276045        41276045      C                 C                  


Tumor_Seq_Allele2  Variant_Classification  Variant_Type  Tumor_Sample_Barcode
A                  Missense_Mutation       SNP           TCGA-AB-1234          
T                  Missense_Mutation       SNP           TCGA-CD-5678          
-                  Frame_Shift_Del         DEL           TCGA-EF-9012          


t_ref_count  t_alt_count  n_ref_count  n_alt_count
50           45           100          0
60           30           120          0
70           20           140          0
```

## sample list file

To select samples of interest in the mutation data file. In the get-tumor-specific mode, two columns in TSV are required.

Data example

```         
tumor_sample_id1  normal_sample_id1
tumor_sample_id2  normal_sample_id1
```

or

```         
tumor_sample_id1
tumor_sample_id2
```

## Configuration file

A configuration file is required to run this program. Write according to the following format. Example of config file.

```         
[SampleList]
sample_file = "xxx/sample_list.tsv"

[Columns]
sample_id_column_index = 10
hugo_symbol_column_index = 1
chromsome_column_index = 2
position_column_index = 3
t_ref_column_index = 5
t_alt_column_index = 7
t_ref_count_column_index = 11
t_alt_count_column_index = 12
vaf_column_index = 15
t_total_count_column_index = 14
other_column_index = [22, 26]

n_alt1_column_index = 20
n_alt2_column_index = 21
n_alt_count_column_index = 23

[VafMode]
rewrite_vaf_mode = "default" or "force"

[OutputFormat]
output_format = "all" or "part"

[GeneSelection]
driver_gene_include_filter = "$HugoSymbol ~ /TP53/ || $HugoSymbol ~ /NF1/"
driver_gene_exclude_filter = "$VAF <= 0.01"
duplicate_gene_mode = "vaf_max"

passenger_gene_include_filter = "$Chr >= 4"
passenger_gene_exclude_filter = "($VAF < 0.23 && $SampleID ~ /T1_WES/) || ($VAF < 0.50  && $SampleID ~ /T2_WES/)"

[PassengerGeneNameFormat]
passenger_gene_name_format = "_allGenes_$Chr_ref_$Ref_alt$Alt"
```

+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Section                 | Item                          | Description                                                                                                                                            |
+:========================+:==============================+:=======================================================================================================================================================+
| SampleList              | sample_file                   | The path to the sample list file                                                                                                                       |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| Columns                 | sample_id_column_index        | Column number for the sample ID (e.g., Tumor_Sample_Barcode)                                                                                           |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | hugo_symbol_column_index      | Column number for the gene name (Hugo_Symbol [e.g., TP53])                                                                                             |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | chromsome_column_index        | Column number for the chromosome number (e.g., 1, X, Y)                                                                                                |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | position_column_index         | Column number for the mutation start position                                                                                                          |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | t_ref_column_index            | Column number for the reference allele (e.g., A)                                                                                                       |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | t_alt_column_index            | Column number for the alteartion allele                                                                                                                |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | t_ref_count_column_index      | Column number for the count of the reference allele.                                                                                                   |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | t_alt_count_column_index      | Column number for the count of the alteration allele                                                                                                   |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | t_total_count_column_index    | Column number for the total count of the alleles                                                                                                       |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | vaf_column_index              | Column number for the VAF.                                                                                                                             |
|                         |                               |                                                                                                                                                        |
|                         |                               | If there is no column for VAF in the input sample, it can remain blank (or not listed).                                                                |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | filter_column_index           | Column number for the FILTER.                                                                                                                          |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | n_alt1_column_index           | Column number for the reference allele of normal sample, only when this allele is defined in different column from that of tumor sample                |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | n_alt2_column_index           | Column number for the alteration allele of normal sample, only when this allele is defined in different column from that of tumor sample               |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | n_alt_count_column_index      | Column number for the count of the alteration allele of normal sample, only when this count is defined in different column from that of tumor sample   |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | other_column_index            | Column numbers to include in the final output file other than the above. Multiple columns possible.                                                    |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| VafMode                 | rewrite_vaf_mode              | Only valid in all mode and get-vaf mode.                                                                                                               |
|                         |                               |                                                                                                                                                        |
|                         |                               | -   default:                                                                                                                                           |
|                         |                               |     -   If there is no VAF column in the input data file, the calculated VAF column is added                                                           |
|                         |                               |     -   If there is a VAF column, the original VAF value is retained.                                                                                  |
|                         |                               | -   force: overwrite the VAF column.                                                                                                                   |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| OutputFormat            | output_format                 | Only valid in get-driver-genes and get-passenger-genes modes.                                                                                          |
|                         |                               |                                                                                                                                                        |
|                         |                               | -   all: all columns of the input file are output in the output file.                                                                                  |
|                         |                               | -   part: only the columns specified in the [Columns] section are output.                                                                              |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| GeneSelection           | driver_gene_include_filter    | Conditions to include mutations into samples.Rint.txt in awk pattern format.                                                                           |
|                         |                               |                                                                                                                                                        |
|                         |                               | Variables shown in the subsequent table can be used to represent columns specified in [Columns].                                                       |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | driver_gene_exclude_filter    | Conditions to exclude mutations from samples.Rint.txt in awk pattern format.                                                                           |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | duplicate_gene_mode           | vaf_max (or vaf_min): If there are multiple SNVs/indels for a possible driver gene, one with the highest (or lowest) VAF is selected. Default: vaf_max |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | passenger_gene_include_filter | Conditions to include mutations into samples.Rother.txt in awk pattern format.                                                                         |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
|                         | passenger_gene_exclude_filter | Conditions to exclude mutations from samples.Rother.txt in awk pattern format.                                                                         |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+
| PassengerGeneNameFormat | PassengerGeneNameFormat       | Convert the gene name (Hugo_Symbol) of passenger genes into the name specified in this item.                                                           |
|                         |                               |                                                                                                                                                        |
|                         |                               | Variables shown in the subsequent table can be used to represent columns specified in [Columns].                                                       |
+-------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------+

### Variables usable in GeneSelection and PassengerGeneNameFormat

| Variables    | Columns                    | Description                  |
|--------------|----------------------------|------------------------------|
| \$SampleId   | sample_id_column_index     | Sample ID                    |
| \$HugoSymbol | hugo_symbol_column_index   | Gene name                    |
| \$Chr        | chromsome_column_index     | Chromosome number            |
| \$Pos        | position_column_index      | Position                     |
| \$Ref        | t_ref_column_index         | Reference allele             |
| \$Alt        | t_alt_column_index         | Alteration allele            |
| \$Rcount     | t_ref_count_column_index   | Number of reference alleles  |
| \$Acount     | t_alt_count_column_index   | Number of alteration alleles |
| \$TotalCount | t_total_count_column_index | Total read count for alleles |
| \$VAF        | vaf_column_index           | VAF value                    |
| \$Filter     | filter_column_index        | FILTER                       |

## Program execution

In the following example, the mutation data file is specified with "data_mutation.txt" and configuration file is specified with "config.conf".

[All Mode]

```         
$ python prepare_observed_data.py -i data_mutation.txt -c config.conf -o ./output/
```

-   This mode executes all steps over to obtain samples.Rint.txt and samples.Rother.txt.

[get-tumor-specific Mode]

```         
$ python prepare_observed_data.py -m get-tumor-specific
  -i data_mutation.txt -c config.conf -o ./output/
```

-   Subtracts variants in a user-specified normal from variants in a user-specified tumor.
-   Data on normal and tumor samples are required in the input files.

[get-vaf Mode]

```         
$ python prepare_observed_data.py -m get-vaf
  -i data_mutation_tumor-specific.txt -c config.conf -o ./output/
```

-   Calculates the VAF (variant allele frequency) for each mutation.
-   For each mutation in the input file, calculate VAF = t_alt_count / ( t_alt_count + t_ref_count ).

[get-driver-genes Mode]

```         
$ python prepare_observed_data.py -m get-driver-genes
  -i data_mutation_tumor-specific_vaf.txt -c config.conf -o ./output/
```

-   Outputs a file corresponding to samples.Rint.txt.
-   Select or exclude mutations that meet the conditions specified in the configuration file for columns such as gene name and VAF.

[get-passenger-genes Mode]

```         
$ python prepare_observed_data.py -m get-passenger-genes
  -i data_mutation_tumor-specific_vaf.txt -c config.conf -o ./output/
```

-   Outputs a file corresponding to samples.Rother.txt.
-   Select or exclude mutations that meet the conditions specified in the configuration file for columns such as gene name and VAF.
-   Gene names are converted into the names specified in the configuration file.

+--------------+-------------------+-----------------------------------------------------------------------+
| Option       | Required/Optional | Description                                                           |
+:=============+:=================:+:======================================================================+
| -m, --model  | Optional          | -   Default: all steps are executed in order.                         |
|              |                   | -   get-tumor-specific                                                |
|              |                   | -   get-vaf                                                           |
|              |                   | -   get-driver-genes                                                  |
|              |                   | -   get-passenger-genes                                               |
+--------------+-------------------+-----------------------------------------------------------------------+
| -i, --input  | Required          | Mutation data file                                                    |
+--------------+-------------------+-----------------------------------------------------------------------+
| -c, --config | Required          | Configuration file.                                                   |
+--------------+-------------------+-----------------------------------------------------------------------+
| -o, --output | Optional          | Output directory name. If not present, the current directory is used. |
+--------------+-------------------+-----------------------------------------------------------------------+

## Output files

The following example shows output files for all modes, when the output directory is specified with "./output/" and mutation data file is specified with "data_mutation.txt".

```         
 -./output/
    |-data_mutation_tumor-specific.txt
    |-data_mutation_tumor-specific_vaf.txt
    |-data_mutation_tumor-specific_vaf.Rint.txt
    |-data_mutation_tumor-specific_vaf.Rother.txt
    |-prepare_observed_data_20241203_121505.log
```

| Output File | Description |
|:---|:---|
| [Base name of input file]\_tumor-specific.txt | Output file of get-tumor-specific mode |
| [Base name of input file]\_vaf.txt | Output file of get-vaf mode |
| [Base name of input file].Rint.txt | Output file of get-driver-genes mode |
| [Base name of input file].Rother.txt | Output file of get-passenger-genes mode |
| xxxxxxxx.log | Log file |

-   For output files of get-driver-genes step and get-passenger-genes modes

    -   If output_format = "all",\
        all columns of the input file are output (a VAF column is added if it does not exist).
    -   If output_format = "part",\
        only the columns specified in the [Columns] section of the configuration file are output.
