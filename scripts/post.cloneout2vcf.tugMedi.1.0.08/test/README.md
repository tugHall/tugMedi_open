# Test Procedures

1. Create test data
2. Run pytest
3. Check test results

## 1. Create Test Data

As an example, use the files extracted from `testdata/0003.zip`.

**Extract `testdata/0003.zip`**

```bash
$ cd testdata
$ unzip 0003.zip
```

**Save Homo_sapiens.GRCh38.dna.primary_assembly.fa to the testdata directory**

```bash
$ wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ./testdata/
```

**Run cloneout2vcf.py**

Procedure: execute test/test_convert_to_vcf.sh
```bash
$ cd ../
$ bash test/test_convert_to_vcf.sh
```
This outputs a vcf file. Additionally, a vcf.ori file will be created as a backup of the output vcf.

Note that if the seed is not fixed in cloneout2vcf.py, the output VCF file from test_convert_to_vcf.sh will differ from the example output (TEST_both/output_max_VAF_primary.vcf).

## 2. Run pytest

**Execute `pytest`**

All tests under test/ will be executed with the following command.
(Currently, only test_convert_to_vcf.py is implemented)
```bash
$ pytest 
```
This creates a vcf file for validation as named vcf.pytest.


If necessary, change the parameters in the following part of `test/test_convert_to_vcf.py`

Note that,
`m = 10  # Data amplification ratio`
is the ratio to pseudo-increase the number of bases in pointMutations.txt used as input.

```python
m = 10  # Data amplification ratio
cloneout_filepath = "./testdata/0003/Output/cloneout.txt"
pointMutations_filepath = "./testdata/0003/Output/Mutations/pointMutations.txt"
VAF_filepath = "./testdata/0003/Output/VAF/VAF.txt"
fasta_filepath = "./testdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
alt_weight_filepath = "./TEST_both/testdata/TCGA-AZ-6608-01A-11D-1835-10_may_alt_weight.txt"
output_dir = './TEST_both'
seed = True
```
To deactivate seed fixation, configure the 'seed' flag as False.


## 3. Check Test Results

This script tests the `convert_to_vcf` function of `cloneout2vcf.py`.  
It compares whether the transition probabilities of each base match the input Weight and the transition probabilities calculated by the `convert_to_vcf` function.
For validation, the Wilson's score method is used here.
The script calculates the 95% confidence interval using the Wilson's score method and checks whether the probabilities calculated by the `convert_to_vcf` function falls within that confidence interval.

### a. Console Output

If the probabilities fall within the confidence interval, "1 passed" is displayed, confirming that there is no problem with the output of the convert_to_vcf function.
If the probabilities do not fall within the confidence interval, "1 failed" is displayed, and the calculated values of each transition, the confidence interval, and a flag indicating whether the probability falls within the confidence interval are displayed.

Example: Display example when there is a transition that does not fall within the confidence interval
```
   REF ALT  count     ratio  ...     lower    center     upper  in_interval
0    A   A      0  0.000000  ...  0.000000  0.007014  0.014028         True
1    A   T    196  0.725926  ...  0.640435  0.694901  0.749367         True
2    A   G     33  0.122222  ...  0.103241  0.144592  0.185942         True
3    A   C     41  0.151852  ...  0.123541  0.167521  0.211501         True
4    T   A    176  0.606897  ...  0.517489  0.574020  0.630550         True
5    T   T      0  0.000000  ...  0.000000  0.006537  0.013073         True
6    T   G     36  0.124138  ...  0.135597  0.179249  0.222901        False
7    T   C     78  0.268966  ...  0.203651  0.253268  0.302886         True
8    G   A    126  0.630000  ...  0.598736  0.663526  0.728316         True
9    G   T     35  0.175000  ...  0.134528  0.187814  0.241100         True
10   G   G      0  0.000000  ...  0.000000  0.009423  0.018845         True
11   G   C     39  0.195000  ...  0.108425  0.158082  0.207740         True
12   C   A     15  0.046875  ...  0.037901  0.064057  0.090212         True
13   C   T    253  0.790625  ...  0.746452  0.790629  0.834805         True
14   C   G     52  0.162500  ...  0.112446  0.151245  0.190045         True
15   C   C      0  0.000000  ...  0.000000  0.005931  0.011862         True
[16 rows x 11 columns]
================================================== short test summary info ==================================================
FAILED test/test_convert_to_vcf.py::test_convert_to_vcf - AssertionError: There are transitions not included in the confidence interval.
===================================================== 1 failed in 3.09s ==================================================
```

Here, `ratio` is the transition probability calculated by the `convert_to_vcf` function, and `center` is the center value of the confidence interval calculated using the Weight with the Wilson's score method. `in_interval` is a flag indicating whether it falls within the confidence interval.

### b. Output Result File: `result_confidence_interval.tsv`

In the output folder specified in test/test_convert_to_vcf.py
```
output_dir = './TEST_both'
```
The table column names in `result_confidence_interval.tsv` are as follows.

- REF
- ALT
- count: Number of transitions 
- ratio: Transition probability
- Weight: Weight used for input
- ref_num: Total number of bases by reference base
- expected_num: Expected number of transitions (ref_num * Weight)
- lower: Lower limit of the confidence interval calculated by the Wilson score method
- center: Center value of the confidence interval calculated by the Wilson score method
- upper: Upper limit of the confidence interval calculated by the Wilson score method
- in_interval: Flag indicating whether it falls within the confidence interval
