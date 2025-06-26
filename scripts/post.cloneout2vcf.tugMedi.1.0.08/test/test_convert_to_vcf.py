import math
import os
import shutil
import sys
from pathlib import Path

# Add the project's root directory to the path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from cloneout2vcf import convert_to_vcf, create_time_prefix


def calc_weight_new_df(vcf_df):
    """
    Function to calculate the weight of base substitutions from a VCF dataframe
    
    Args:
        vcf_df (pd.DataFrame): Dataframe in VCF format
    
    Returns:
        pd.DataFrame: Dataframe containing the occurrence count and ratio of each base substitution pair
    """
    # Extract base substitution pairs from REF and ALT columns
    ref_col = vcf_df['REF']
    alt_col = vcf_df['ALT']

    # Count the occurrence of base substitution pairs
    pair_counts = vcf_df.groupby([ref_col, alt_col]).size().reset_index(name='count')

    # Define all possible base substitution pairs (all combinations of A, T, G, C)
    bases = ['A', 'T', 'G', 'C']
    all_pairs = [(ref, alt) for ref in bases for alt in bases]

    # Set initial count of all base substitution pairs to 0
    pair_counts_dict = {pair: 0 for pair in all_pairs}

    # Update with actual occurrence counts
    for _, row in pair_counts.iterrows():
        pair = (row['REF'], row['ALT'])
        if pair in pair_counts_dict:
            pair_counts_dict[pair] = row['count']

    # Convert results to DataFrame
    ordered_pair_counts = pd.DataFrame(
        [(ref, alt, count) for (ref, alt), count in pair_counts_dict.items()],
        columns=['REF', 'ALT', 'count']
    )

    # Calculate the occurrence ratio of ALT bases for each REF base
    ref_totals = ordered_pair_counts.groupby('REF')['count'].sum()
    ordered_pair_counts['ratio'] = ordered_pair_counts.apply(
        lambda row: row['count'] / ref_totals[row['REF']] if ref_totals[row['REF']] != 0 else 0,
        axis=1
    )
    return ordered_pair_counts


def calc_confidence_interval(res_df_weight, alt_weight_df):
    """
    Function to calculate base substitution weights and the confidence interval of theoretical values
    
    Args:
        res_df_weight (pd.DataFrame): Dataframe containing the weight of base substitutions
        alt_weight_df (pd.DataFrame): Dataframe containing theoretical base substitution probabilities
    
    Returns:
        pd.DataFrame: Dataframe containing the results of confidence interval calculations
    """
    # Merge actual and theoretical data
    res_df_ratio_weight = res_df_weight.merge(
        alt_weight_df, 
        left_on=['REF', 'ALT'], 
        right_on=['Ref', 'Alt'], 
        how='left'
    )
    res_df_ratio_weight = res_df_ratio_weight.drop(columns=['Ref', 'Alt'])

    # Calculate the total occurrence count for each REF base
    ref_counts = res_df_ratio_weight.groupby('REF')['count'].sum()
    res_df_ratio_weight['ref_num'] = res_df_ratio_weight['REF'].map(ref_counts)

    # Calculate expected values
    res_df_ratio_weight['expected_num'] = res_df_ratio_weight['Weight'] * res_df_ratio_weight['ref_num']

    # Calculate confidence intervals using the Wilson's score method
    res_df_ratio_weight[['lower', 'center', 'upper']] = res_df_ratio_weight.apply(
        lambda row: pd.Series(wilson_score_interval(row['Weight'], row['ref_num'], confidence=0.95)),
        axis=1
    )

    # Adjust the boundaries of the confidence interval
    res_df_ratio_weight['lower'] = res_df_ratio_weight['lower'].apply(lambda x: 0.0 if x < 0.0000000001 else x)
    res_df_ratio_weight['upper'] = res_df_ratio_weight['upper'].apply(lambda x: 1.0 if x > 0.9999999999 else x)

    # Determine if the actual values fall within the confidence interval
    res_df_ratio_weight['in_interval'] = res_df_ratio_weight.apply(
        lambda row: row['lower'] <= row['ratio'] <= row['upper'],
        axis=1
    )
    return res_df_ratio_weight


def wilson_score_interval(successes_p, trials, confidence=0.95):
    """
    Calculate the confidence interval of a binomial distribution using the Wilson's score method.
    
    Parameters:
        successes_p (float): Success ratio (p)
        trials (int): Number of trials (n)
        confidence (float): Confidence level (default is 95%)
    
    Returns:
        tuple: Lower bound, center, and upper bound of the confidence interval
    """
    if trials == 0:
        raise ValueError("The number of trials must be greater than 0.")
    
    # Calculate the Z value corresponding to the confidence level
    z = norm.ppf(1 - (1 - confidence) / 2)
    
    # Calculate the Wilson's score method
    denominator = 1 + (z**2 / trials)
    center = (successes_p + (z**2 / (2 * trials))) / denominator
    margin = z * math.sqrt((successes_p * (1 - successes_p) + (z**2 / (4 * trials))) / trials) / denominator
    
    lower_bound = center - margin
    upper_bound = center + margin
    
    return lower_bound, center, upper_bound


def test_convert_to_vcf():
    """
    Test function for VCF conversion
    Amplify input data, perform conversion, and verify the validity of the results

    Execution command:
    python -m pytest test/test_convert_to_vcf.py
    """
    m = 1  # Data amplification factor
    cloneout_filepath = "./testdata/0003/Output/cloneout.txt"
    pointMutations_filepath = "./testdata/0003/Output/Mutations/pointMutations.txt"
    VAF_filepath = "./testdata/0003/Output/VAF/VAF.txt"
    fasta_filepath = "./testdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    alt_weight_filepath = "./TEST_both/testdata/TCGA-AZ-6608-01A-11D-1835-10_may_alt_weight.txt"
    output_dir = './TEST_both'

    # Load theoretical data and expand input data
    alt_weight_df = pd.read_csv(alt_weight_filepath, sep="\t")
    pointMutations_data = pd.read_csv(pointMutations_filepath, sep='\t')
    pointMutations_data_expanded = pd.concat([pointMutations_data] * m, ignore_index=True)
    
    # Temporarily save expanded data
    pointMutations_filepath_tmp = output_dir + "/pointMutations_tmp.txt"
    pointMutations_data_expanded.to_csv(pointMutations_filepath_tmp, sep='\t', index=False)

    # Execute VCF conversion
    result_df = convert_to_vcf(
        time=['max'],
        time_filter_logic='or',
        cloneout_filepath=cloneout_filepath,
        pointMutations_filepath=pointMutations_filepath_tmp,
        VAF_filepath=VAF_filepath,
        VAF_column='VAF_primary',
        fasta_filepath=fasta_filepath,
        output_path=output_dir,
        alt_weight_df=alt_weight_df,
        seed=True
    )

    # Get time prefix for file operations
    time_prefix = create_time_prefix(['max'], 'or')
    VAF_column = 'VAF_primary'
    
    # 1. Rename the original VCF file to include amplification factor
    vcf_file = Path(f'{output_dir}/output_{time_prefix}_{VAF_column}.vcf')
    new_vcf_file = Path(f'{output_dir}/output_{time_prefix}_{VAF_column}.vcf.pytest')
    if vcf_file.exists():
        shutil.move(str(vcf_file), str(new_vcf_file))
        print(f"Renamed {vcf_file} to {new_vcf_file}")
    
    # 2. Rename .ori file to original VCF filename
    ori_file = Path(f'{output_dir}/output_{time_prefix}_{VAF_column}.vcf.ori')
    if ori_file.exists():
        shutil.copy(str(ori_file), str(vcf_file))
        print(f"Renamed {ori_file} to {vcf_file}")


    ref_counts = result_df['REF'].value_counts()
    res_df_weight = calc_weight_new_df(result_df)

    # Merge Weight and transition probability table, and calculate confidence intervals
    res_df_confidence_interval = calc_confidence_interval(res_df_weight, alt_weight_df)

    # Output res_df_confidence_interval to output_dir
    res_df_confidence_interval.to_csv(output_dir + "/result_confidence_interval.tsv", sep='\t', index=False)
    print(res_df_confidence_interval)

    # Verify results
    assert result_df is not None, "The result DataFrame is None"
    assert not result_df.empty, "The result DataFrame is empty"

    # Check if within confidence interval
    assert res_df_confidence_interval['in_interval'].all(), "There are transitions not included in the confidence interval."

