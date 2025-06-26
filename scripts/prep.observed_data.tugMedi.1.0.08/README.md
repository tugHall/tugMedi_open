# Tool for preparing observed data

**Tool Name:** prepare_observed_data.py

Example Command:

`$ python prepare_observed_data.py -i data_mutation.txt -c config.conf -o ./output/`

## Input data format

### Mutation data file (TEXT)

* Input: Required  
* Number of inputs: 1  
* Extension: None specified  
* File containing mutation information of nucleotide sequences (e.g., assumed to be similar to MAF)

The following items are required. Note that the columns required will change depending on the calculation mode.
```
  * Sample ID   
  * Gene name (Hugo_Symbol)  
  * Chromosome: Chromosome number (e.g., 1, X)  
  * Position: Position on the chromosome (e.g., 7674220)  
  * Ref: Reference sequence (e.g., A)  
  * Alt: Mutation sequence of the sample (e.g. G)
  * Ref_count: Count of reference alleles in thesample
  * Alt_count: Count of mutated alleles in the sample
  * VAF: Variant Allele Frequency
  
  The following is required for formats that have normal sample data in a separate column,
  * n_alt1: Mutation sequence 1 of the normal sample
  * n_alt2: Mutation sequence 2 of the normal sample
  * n_alt_count: Count of mutant alleles in the normal sample  
```

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

## sample list file format
Used to specify samples to be used in the calculation of this script in the mutation data file.
A sample list file is in tsv format and is required in get-tumor-specific mode.
Write according to the following format.
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


## Configuration file format

A configuration file is required to run this software. Write according to the following format.
Example of config file.
```
[SampleList]
sample_file = “xxx/sample_list.tsv”

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

| Section | Item | Description |
| :---- | :---- | :---- |
| SampleList | sample_file | The path to the sample list file |
| Columns | sample_id_column_index  | Column number where the sample ID (e.g., Tumor_Sample_Barcode) is listed. |
|  | hugo_symbol_column_index | Column number where the gene name (Hugo_Symbol: e.g., TP53) is listed. |
|  | chromsome_column_index | Column number where the chromosome number (e.g., 1, X, Y) is listed.  |
|  | position_column_index | Column number where the mutation start position is listed. |
|  | t_ref_column_index | Column number where the reference allele (e.g., A) is listed. |
|  | t_alt_column_index | Column number where the mutation allele is listed. |
|  | t_ref_count_column_index | Column number where the count of the reference allele is listed. |
|  | t_alt_count_column_index | Column number where the count of the mutation allele is listed. |
|  | t_total_count_column_index | Column number where the count of the total allele is listed. |
|  | vaf_column_index | Column number where the VAF is listed. If there is no column for VAF in the input sample, it can remain blank (or not listed). |
|  | filter_column_index | Column number where the FILTER is listed. |
|  | n_alt1_column_index | Column number where the mutation allele of normal sample is listed. |
|  | n_alt2_column_index | Column number where the another mutation allele of normal sample is listed. |
|  | n_alt_count_column_index | Column number where the count of the mutation allele of normal sample is listed. |
|  | other_column_index  | Column numbers to include in the final output file other than the above. Multiple entries possible. |
| VafMode | rewrite_vaf_mode | Valid in all mode and get-vaf mode. default: If there is no VAF column in the input data file, the calculated VAF column is added; if there is a VAF column, the original VAF value is used. force: overwrite the VAF column.  |
| OutputFormat | output_format | Valid in get-driver-genes and get-passenger-genes modes. all: all columns of the input file are included in the output file. part: only the columns specified in the [Columns] section are output. |
| GeneSelection | driver_gene_include_filter | Describe the conditions that are considered as mutations of the driver gene in awk pattern format. The variables that can be used are the columns specified in [Columns].* |
|  | driver_gene_exclude_filter | Describe the conditions to be excluded from the mutation of the driver gene in awk pattern format. |
|  | duplicate_gene_mode | vaf_max(vaf_min): If there are multiple SNVs/indels for a possible driver gene, the highest (lowest) VAF is selected. |
|  | passenger_gene_include_filter | Describe the conditions that are considered as mutations of the passenger gene in awk pattern format. |
|  | passenger_gene_exclude_filter | Describe the conditions to be excluded from the mutation of the passenger gene in awk pattern format. |
| PassengerGeneNameFormat | PassengerGeneNameFormat | Convert the gene name (Hugo_Symbol) of passenger genes to the name specified in this item. The variables that can be used are the columns specified in [Columns].*|
* $SampleId: Sample ID, $HugoSymbol: Gene name, $Chr: Chromosome number, $Pos: Position, $Ref: Reference allele, $Alt: Mutated allele, $Rcount: Number of reference alleles, $Acount: Number of mutant alleles, $TotalCount: Total read count for alleles, $VAF: VAF value, $Filter: FILTER.

## Program execution image

In the following example, a mutation data file (data_mutation.txt), configuration file (config.conf) are prepared in the current directory.

[All Mode]
```
$ python prepare_observed_data.py -i data_mutation.txt -c config.conf -o ./output/
```
* This mode executes up to the creation of samples.Rint.txt and samples.Rother.txt in one go.
  1. Prepare mutation data file, sample-list file and configuration file.
  2. Run the script using the command line with the appropriate options.
  3. Execute the script to process the data and generate output files (see Output files section).
  

[get-tumor-specific Mode]
```
$ python prepare_observed_data.py -m get-tumor-specific
  -i data_mutation.txt -c config.conf -o ./output/
```
  * To identify tumor-specific mutations.
  * Normal and tumor samples are required in the input file.


[get-vaf Mode]
``` 
$ python prepare_observed_data.py -m get-vaf
  -i data_mutation_tumor-specific.txt -c config.conf -o ./output/
```
  * Execute when calculating the VAF (variant allele frequency) for each mutation.
  * For each mutation in the input file, calculate VAF = t_alt_count / ( t_alt_count + t_ref_count ) using the values of t_ref_count and t_alt_count.


[get-driver-genes Mode]
```
$ python prepare_observed_data.py -m get-driver-genes
  -i data_mutation_tumor-specific_vaf.txt -c config.conf -o ./output/
```
  * Execute when extracting mutation data corresponding to driver genes.
  * Select or exclude mutations that match the conditions specified in the configuration file for columns such as gene name and VAF.


[get-passenger-genes Mode]
```
$ python prepare_observed_data.py -m get-passenger-genes
  -i data_mutation_tumor-specific_vaf.txt -c config.conf -o ./output/
```
  * Execute when extracting mutation data corresponding to passenger genes.
  * Select or exclude mutations that match the conditions specified in the configuration file for columns such as gene name and VAF.
  * Gene names are converted to the names specified in the configuration file.


| Option | Required/Optional | Description |
| :---- | :---: | :---- |
| -m, --model | Optional | Default: all steps are executed in order. get-tumor-specific: extract tumor-specific mutations. get-vaf: calculate VAF columns. get-driver-genes: extract driver genes. get-passenger-genes: extract passenger genes. |
| -i, --input | Required | Gene mutation data file to input. (similar to MAF format) |
| -c, --config | Required | Specify the configuration file. |
| -o, --output | Optional | Specify the output directory name. If not present, it will be output to the current directory. |

## Output files

The following example shows the output files when "./output/" is specified as the output directory.

Assuming the input sample file name is data_mutation.txt, the output files for all modes are as follows.

```
 -./output/
    |-data_mutation_tumor-specific.txt
    |-data_mutation_tumor-specific_vaf.txt
    |-data_mutation_tumor-specific_vaf.Rint.txt
    |-data_mutation_tumor-specific_vaf.Rother.txt
    |-prepare_observed_data_20241203_121505.log

```

| Output File | Description |
| :---- | :---- |
| [Base name of input file]_tumor-specific.txt | Output file of get-tumor-specific step |
| [Base name of input file]_vaf.txt | Output file of get-vaf step |
| [Base name of input file].Rint.txt | Output file of get-driver-genes step |
| [Base name of input file].Rother.txt | Output file of get-passenger-genes step |
| xxxxxxxx.log | Log file |

* Data format (output files of get-driver-genes step and get-passenger-genes step)

  * Tab-separated (TSV) format.  
  * If output_format = "all" is selected, the format of the input file is output as is (a VAF column is added if it does not exist).  
  * By output_format = "part" is selected, only the columns specified in the [Columns] section of the configuration file are output.  
