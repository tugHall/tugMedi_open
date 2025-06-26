

showPosts.py, then bayesOptim.R


1. showPosts.py
take in the results of the ABC scripts to make a table of ABC distance and parameter values
* Now, only for the parameters of m0, dN (from Input/parameters.txt), 
  and EF timings of M2, M3, ...          (from Input/EF.Rint.txt)

Usage: 
showPosts.py -h

Example: 
showPosts.py --nn 2000 \
  --file_ABCdist "./ABC2/TCGA-55-7903-01A-11D-2167-08/ABC_dist_tumor_content=0.80.txt" \
  --dir_baseReps "./work/TCGA-55-7903-01A-11D-2167-08/" \
  > work/prmDist.txt

Input: 
  -nn: top nn based on ABC distance (ascending order, ie, shorter, better)
  --file_ABCdist: an output of the ABC scripts, the ABC distance file
  --dir_baseReps: an output of the ABC scripts, the directory with simulation replicates such as 000001, 000002, ...
Output: 
  1st col: simulation replicate ID
  2nd col: ABC distance
  3rd col: m0
  4th col: dN
  5th col: M2
  6th col: M3
  ...    : M...


2. bayesOptim.R
take in the output of showPosts.py to perform Bayesian optimization

Usage: 
bayesOptim.R

Example: 
bayesOptim.R work/prmDist.txt 1 'c(3,4)' T > work/prmOpt1.T.txt 2>&1 
bayesOptim.R work/prmDist.txt 2 'c(3,4)' T > work/prmOpt2.T.txt 2>&1 
bayesOptim.R work/prmDist.txt 1 'c(3,4)' F > work/prmOpt1.F.txt 2>&1 
bayesOptim.R work/prmDist.txt 2 'c(3,4)' F > work/prmOpt2.F.txt 2>&1 

Input: 
  1st arg: output of showPosts.py
  2nd arg: 1, standard mode; 2, trial mode; 3, random search
  3rd arg: col number to transform into log10 in 1st arg, eg, 3rd (m0) and 4th (dN) cols
  4th arg: whether to perform nugget estimation, see the usage
Output:
  Output from R, 
    where the last 2 lines show optimal parameter values and distance 
    (tail -n2 work/prmOpt* |less)



