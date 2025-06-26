#!/usr/bin/env python
import pandas as pd
import numpy as np
import random
from pyfaidx import Fasta
from pathlib import Path
import datetime
import argparse
import os
import sys
#from calc_weight_from_maf import alt_weight_df

date_now = datetime.datetime.now().strftime("%Y%m%d")
VCF_HEADER = f"""##fileformat=VCFv4.2
##fileDate={date_now}
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=Time,Number=1,Type=Integer,Description="Time">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
""".strip()


pd.set_option('display.max_rows', 200)
################### from calc_weight_from_maf.py ###################
init_alt_weight_df = pd.DataFrame({'Ref': ['A', 'A', 'A', 'A', 'T', 'T', 'T', 'T', 'G', 'G', 'G', 'G', 'C', 'C', 'C', 'C'],
                              'Alt': ['A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C'],
                              'Weight': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]})



def assign_weight(ref, alt, df, init_df):
    """Preprocess function for calc_weight_from_TCGA.
    """
    total_mutations = df.loc[df['Reference_Allele']==ref, 'count'].sum()
    ref_to_alt = df.loc[(df['Reference_Allele']==ref) & (df['Tumor_Seq_Allele2']==alt), 'count'].iloc[0]
    init_df.loc[(init_df['Ref']==ref) & (init_df['Alt']==alt), 'Weight'] = float(ref_to_alt/total_mutations)


def calc_weight_from_TCGA(maf_filepath, init_df):
    """Calculate alteration weights from TCGA data.

    Args:
        maf_filepath(maf formatted filepath): TCGA data to calculate alteration weights.
    
    Return:
        alt_weight_df(pandas.DataFrame): 16row*3col('Ref', 'Alt', 'Weight') base alteration weights data including 'A'to'A', 'T'to'T', 'G'to'G', 'C'to'C' in weight 0.0.

    """
    #maf_df = pd.read_csv(maf_filepath, delimiter="\t", comment='#', usecols=['Reference_Allele', 'Tumor_Seq_Allele2'], na_values=["-"], encoding="shift-jis")
    maf_df = pd.read_csv(maf_filepath, delimiter="\t", comment='#', na_values=["-"], encoding="shift-jis")
    maf_df.iloc[0]
    selected_df = maf_df.groupby('Reference_Allele')['Tumor_Seq_Allele2'].value_counts()
    selected_df = selected_df.rename("count")
    selected_df = selected_df.reset_index()
    # 特定の塩基へのmutationのデータがなかった場合にそのmutationのcount 0 を追加する
    if len(selected_df) < 12:
        mapping_dict = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
        reverse_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}
        selected_df['Allele_id'] = selected_df["Tumor_Seq_Allele2"].map(mapping_dict)
        A_alt_counter = [0 for _ in range(4)]
        T_alt_counter = [0 for _ in range(4)]
        G_alt_counter = [0 for _ in range(4)]
        C_alt_counter = [0 for _ in range(4)]
        for i in range(len(selected_df)):
            if selected_df.loc[i, 'Reference_Allele'] == 'A':
                A_alt_counter[selected_df.loc[i, 'Allele_id']] += 1
            elif selected_df.loc[i, 'Reference_Allele'] == 'T':
                T_alt_counter[selected_df.loc[i, 'Allele_id']] += 1
            elif selected_df.loc[i, 'Reference_Allele'] == 'G':
                G_alt_counter[selected_df.loc[i, 'Allele_id']] += 1
            elif selected_df.loc[i, 'Reference_Allele'] == 'C':
                C_alt_counter[selected_df.loc[i, 'Allele_id']] += 1
        counters = [A_alt_counter, T_alt_counter, G_alt_counter, C_alt_counter]
        for i in range(len(counters)):
            counter = counters[i]
            if min(counter) == 0:
                add_data = {'Reference_Allele':[],
                            'Tumor_Seq_Allele2':[],
                            'count':[]}
                base = reverse_dict[i]
                for idx in range(len(counter)):
                    if counter[idx] == 0:
                        add_data['Reference_Allele'].append(base)
                        add_data['Tumor_Seq_Allele2'].append(reverse_dict[idx])
                        add_data['count'].append(0)
                selected_df = pd.concat([selected_df, pd.DataFrame(add_data)])
    bases = ['A', 'T', 'G', 'C']
    for ref in bases:
        for alt in bases:
            if ref != alt:
                assign_weight(ref, alt, selected_df, init_df)
    return init_df

#####################################################################

# preprocess functions
def check_path(dir):
    if not Path(dir).is_dir():
        Path(dir).parents[0].mkdir(parents=True, exist_ok=True)
    return dir


def check_columns_for_nan(df, column_name, correspond_col, joined_data_name):
    """Function to raise exception.
    Check new columns when join data. If new column contains only nan, raise ValueError.
    """
    for col in column_name:
        if df[col].isnull().all():
            raise ValueError("ValueError: {} column contains only NaN values.Check {} column corresponds to {}:{} in {}.".format(column_name, column_name, correspond_col, df[correspond_col].tolist(), joined_data_name)) 


def flatten(lst):
    new_lst = []
    for elm in lst:
        if type(elm) == list:
            new_lst = new_lst + flatten(elm)
        else:
            new_lst.append(int(elm))
    return new_lst


def unique(lst):
    unique_lst = []
    for elm in lst:
        if elm not in unique_lst:
            unique_lst.append(elm) 
    return unique_lst


def check_df_format(df):
    """Function to raise exception.
    """
    if df.shape != (16, 3):
        raise ValueError('Base alteration weights data size should be (16, 3).')
    

def calc_weight(alt_weights_data, base):
    """Preprocess function for assign_alt_bases()

    Args:
        alt_weights_data(pandas.DataFrame): 16row*3col('Ref', 'Alt', 'Weight') base alteration weights data.
                                            Should include 'A'to'A', 'T'to'T', 'G'to'G', 'C'to'C' in weight 0.0.
        base(str): 'A' or 'T' or 'G' or 'C'

    Return:
        base_alt_bases(str list): Ordered four bases
        base_alt_weights(float list): Ordered four weights

    Exsample:
        >>> calc_alt_bases(base_alt_df, 'A')
        (['A', 'T', 'G', 'C'], [0.0, 0.2, 0.4, 0.4]) 
    """
    try:
        check_df_format(alt_weights_data)
    except ValueError as e:
        print(e)
        return
    base_alt_bases = alt_weights_data.loc[alt_weights_data['Ref'] == base, ['Alt']]['Alt'].tolist()
    base_alt_weights = alt_weights_data.loc[alt_weights_data['Ref'] == base, ['Weight']]['Weight'].tolist()
    return  base_alt_bases, base_alt_weights


def assign_alt(df, alt_weight_df): 
    A_alt_base, A_alt_weights = calc_weight(alt_weight_df, 'A') 
    T_alt_base, T_alt_weights = calc_weight(alt_weight_df, 'T')
    G_alt_base, G_alt_weights = calc_weight(alt_weight_df, 'G')
    C_alt_base, C_alt_weights = calc_weight(alt_weight_df, 'C')
    if df['REF'] == 'A':
        alt = random.choices(A_alt_base, A_alt_weights)[0]
    elif df['REF'] == 'T':
        alt = random.choices(T_alt_base, T_alt_weights)[0]
    elif df['REF'] == 'G':
        alt = random.choices(G_alt_base, G_alt_weights)[0]
    elif df['REF'] == 'C':
        alt = random.choices(C_alt_base, C_alt_weights)[0]
    return alt

#
# add time filter functions
#
def parse_time_filter(time_args, target_df=None, target_column='Time'):
    """
    when max/min is specified, target_df is needed to get max/min value.
    """
    filters = []
    for arg in map(str, time_args): #ensure all arguments are strings
        if arg.lower() == 'max':
            #get max value from target_df
            max_time = max(target_df[target_column])
            filters.append(lambda x, max_time=max_time: x == max_time)
        elif arg.lower() == 'min':
            #get min value from target_df
            min_time = min(target_df[target_column])
            filters.append(lambda x, min_time=min_time: x == min_time)
        elif arg.startswith('>='):
            threshold = int(arg[2:].strip())
            filters.append(lambda x, threshold=threshold: x >= threshold)
        elif arg.startswith('=>'):
            threshold = int(arg[2:].strip())
            filters.append(lambda x, threshold=threshold: x >= threshold)
        elif arg.startswith('<='):
            threshold = int(arg[2:].strip())
            filters.append(lambda x, threshold=threshold: x <= threshold)
        elif arg.startswith('=<'):
            threshold = int(arg[2:].strip())
            filters.append(lambda x, threshold=threshold: x <= threshold)
        elif arg.startswith('>'):
            threshold = int(arg[1:].strip())
            filters.append(lambda x, threshold=threshold: x > threshold)
        elif arg.startswith('<'):
            threshold = int(arg[1:].strip())
            filters.append(lambda x, threshold=threshold: x < threshold)
        elif arg.startswith('=='):
            threshold = int(arg[2:].strip())
            filters.append(lambda x, threshold=threshold: x == threshold)
        else:
            try:
                value = int(arg)
                filters.append(lambda x, value=value: x == value)
            except ValueError:
                raise ValueError(f"Invalid time filter argument: {arg}")
    return filters


def get_time_filters(data, filters, logic='or', target_column='Time'):
    if logic not in ['and', 'or']:
        raise ValueError("Logic must be either 'and' or 'or'.")
    
    # No filters, return all True
    if len(filters) == 0:
        return [True] * len(data)
    
    # Start with all True or all False depending on the logic
    if logic == 'and':
        mask = pd.Series([True] * len(data), index=data.index)  # All True for AND logic
        for filter_func in filters:
            mask &= data[target_column].apply(filter_func)  # AND logic
    else:  # logic == 'or'
        mask = pd.Series([False] * len(data), index=data.index)  # All False for OR logic
        for filter_func in filters:
            mask |= data[target_column].apply(filter_func)  # OR logic
    
    #just return boolean list
    return mask.to_list()
    #return mask

def create_time_prefix(time, logic='or'):
    if len(time) == 0:
        return "all"

    comparisons = {
        ">": "gt",
        "<": "lt",
        ">=": "ge",
        "<=": "le",
        "==": ""
    }

    tlist=list()
    for one in time:
        one = str(one)
        if len(one) >= 2 and one[:2] in comparisons:
            tlist.append(comparisons[one[:2]] + one[2:].strip())
        elif len(one) >= 1 and one[0] in comparisons:
            tlist.append(comparisons[one[0]] + one[1:].strip())
        else:
            tlist.append(one)
    
    if len(time) > 1:
        tlist.append(logic)

    return "_".join(tlist)

# converter function
def convert_to_vcf(
        time,
        time_filter_logic,
        cloneout_filepath, 
        pointMutations_filepath, 
        VAF_filepath, 
        VAF_column, 
        fasta_filepath, 
        output_path, 
        alt_weight_df,
        seed=False
    ):
    """Convert three tugHall outputs to vcf formatted DataFrame and save as tsv a file.

    convert_to_vcf takes five arguments and save tcv formatted table as tsv file to output_path and return pandas.DataFrame.

    Args:
        time(List(int)): Time corresponds to Time column in cloneout_data.
        time_filter_logic(str): 'and' or 'or'. [default=or]
        cloneout_filepath(tsv formatted filepath): First tugHall output data.
        pointMutations_filepath(tsv formatted filepath): Second tugHall output data.
        VAF_fileppath(tsv formatted filepath): Third tugHall output data.
        VAF_column(str): 'VAF_primary' or 'VAF_metastatic'.
        fasta_filepath(fasta filepath): Genomic data of Homo sapiense to fill 'REF' column.
        output_path(directory): Output directory to save output tsv file.
        alt_weight_df(pd.DataFrame): Output data of calc_weight_from_maf.py. Assign bases to ALT column based on this data.
        seed(bool): Use fixed random seed for reproducibility.

    Save:
        tsv file: tcv formatted tsv file

    Return:
        pandas.DataFrame: tcv formatted dataframe

    Exsample:
        >>> convert_to_vcf(3, cloneout_df, pointMutations_df, VAF_df, 'VAF_metastatic')    
        
       CHROM  ID   POSITION REF ALT  QUAL FILTER        INFO
    0      1   9   74572100   G   T    99   PASS  AF=0.12297
    1      1  11  146139901   T   A    99   PASS  AF=0.12297
    """
    # Set random seed if flag is enabled
    if seed:
        seed_value = 312
        random.seed(seed_value)
        np.random.seed(seed_value)
        print(f"[INFO] Random seed set to: {seed_value}")

    # load_data
    cloneout_data = pd.read_csv(cloneout_filepath, sep='\t')
    pointMutations_data = pd.read_csv(pointMutations_filepath, sep='\t')
    VAF_data = pd.read_csv(VAF_filepath, sep='\t')
    fasta = Fasta(fasta_filepath)
    if 'tumor_content' in VAF_data.columns:
        VAF_data = VAF_data[VAF_data['tumor_content'] == 1]
    
    # Flatten and Explode IDs based on PointMut_ID in VAF_data
    vaf_PointMut_IDs = VAF_data['PointMut_ID']
    VAF_data['flatten_PointMut_ID'] = vaf_PointMut_IDs.apply(lambda x: str(x).split(','))
    VAF_data = VAF_data.explode('flatten_PointMut_ID')
    VAF_data = VAF_data.drop('PointMut_ID', axis='columns')
    VAF_data = VAF_data.rename({'flatten_PointMut_ID' : 'PointMut_ID'}, axis='columns')
    VAF_data["PointMut_ID"] = VAF_data["PointMut_ID"].astype(str)
    # filter VAF_data by time

    VAF_data["Time"] = VAF_data["Time"].astype(int)
    time_filters = parse_time_filter(time, VAF_data)
    time_mask_list = get_time_filters(VAF_data, time_filters, time_filter_logic)
    VAF_data = VAF_data[time_mask_list]

    if len(VAF_data) == 0:
        print(f"After time filtering, VAF_data has no data, exit.")
        return None
    
    # filter cloneout_data by time
    time_filters = parse_time_filter(time, VAF_data)
    time_mask_list = get_time_filters(cloneout_data, time_filters, time_filter_logic)
    cloneout_data = cloneout_data[time_mask_list]

    # extract just PointMut_IDs from cloneout_data
    pointmut_ids = cloneout_data[['PointMut_ID']]
    
    if len(pointmut_ids) == 0:
        print(f"Time {time} has no data, exit.")
        return None
    
    lst_ids = [tid_list for tid_list in pointmut_ids['PointMut_ID'].str.split(',') if tid_list != ['-']]
    new_lst = [str(i) for i in unique(flatten(lst_ids))]
    new_df = pd.DataFrame(new_lst, columns =['PointMut_ID'])
    new_df['PointMut_ID'] = new_df['PointMut_ID'].astype(str)

    #
    # exclude PointMut_IDs that are not in VAF_data
    #
    new_df = new_df[new_df['PointMut_ID'].isin(VAF_data['PointMut_ID'])]
    if len(new_df) == 0:
        print(f"Time {time} has no data, exit.")
        return None

    pointMutations_data['PointMut_ID'] = pointMutations_data['PointMut_ID'].astype(str)
    new_df = pd.merge(new_df, pointMutations_data[['PointMut_ID', 'Chr', 'Ref_pos']], on='PointMut_ID', how='inner')
    #for develop save
    #new_df_merge = pd.merge(new_df, pointMutations_data[['PointMut_ID', 'Chr', 'Ref_pos']], on='PointMut_ID', how='inner')
    try:
        check_columns_for_nan(new_df, ['Chr', 'Ref_pos'], 'PointMut_ID', 'pointMutationsdata')
    except ValueError as e:
        print(e)
        return
    new_df['QUAL'] = 99
    new_df['FILTER'] = 'PASS'

    new_df = pd.merge(new_df, VAF_data[['PointMut_ID', VAF_column, "Time", "gene"]], on='PointMut_ID', how='left')
    try:
        check_columns_for_nan(new_df, [VAF_column], 'PointMut_ID', 'VAF_data')
    except ValueError as e:
        print(e)
        return 
    chromosome_num = new_df['Chr'].astype('str').tolist()
    index = new_df['Ref_pos'].tolist()
    new_df['REF'] = [fasta[chr][idx].seq for chr, idx in zip(chromosome_num, index)]
    new_df['ALT'] = new_df.apply(lambda row : assign_alt(row, alt_weight_df), axis=1)
    #
    # create INFO column. when add new INFO column, add the info in VCF_HEADER global variable.
    #
    new_df['INFO'] = 'AF=' + new_df[VAF_column].round(5).astype('str') + \
                    ';Time=' + new_df['Time'].astype('str') + \
                    ';GENE=' + new_df['gene']

    new_df.rename(columns={'PointMut_ID': 'ID', 'Chr': '#CHROM', 'Ref_pos': 'POS'}, inplace=True)

    new_df = new_df.reindex(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    #sort by CHROM and POS
    # CHROM, POS should be int
    new_df['#CHROM'] = new_df['#CHROM'].astype(int)
    new_df['POS'] = new_df['POS'].astype(int)
    new_df = new_df.sort_values(by=['#CHROM', 'POS', "REF", "ALT"])
    new_df = new_df.reset_index(drop=True)
    #
    # save as vcf
    #
    time_prefix = create_time_prefix(time, time_filter_logic)
    with open(check_path(str(Path(f'{output_path}/output_{time_prefix}_{VAF_column}.vcf'))), 'w') as f:
        f.write(VCF_HEADER + '\n')
        new_df.to_csv(f, sep="\t", index=False)
    return new_df

def get_args():
    parser = argparse.ArgumentParser(description='Convert three tugHall outputs to vcf',
    epilog="""\
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
""",
    formatter_class=argparse.RawTextHelpFormatter  # Allows multiline help text
    )
                                     
    parser.add_argument('--time', nargs='*', default=[], help="""Time corresponds to Time column in VAF_data. multiple values are allowed. if not specified, all times are used. \
also specify max/min as string.
Example: --time 1 2 3
          --time max
          --time min
          --time >0
          --time <100
""")
    parser.add_argument('--time_filter_logic', type=str, default='or', help="Logic to combine time filters. 'and' or 'or'. [default=or]")
    parser.add_argument('--cloneout', type=str, help='First tugHall output data.')
    parser.add_argument('--pointMutations', type=str, help='Second tugHall output data.')
    parser.add_argument('--VAF', type=str, help='Third tugHall output data.')
    parser.add_argument('--VAF_column', type=str, default="VAF_primary" ,help='VAF_primary or VAF_metastatic. [default=VAF_primary]')
    parser.add_argument('--maf', type=str, help='TCGA maf filepath.')
    parser.add_argument('--alt_weight', type=str, help='Pre-calculated alteration weights data.')
    parser.add_argument('--fasta', type=str, help='genome fasta filepath.')
    parser.add_argument('--mode', type=str, default='both', help='Mode to run the program. [convert/calcweight/both, default=both]')
    parser.add_argument('--output', type=str, default="TEST", help='Output directory to save output vcf file. [default=TEST]')
    parser.add_argument('--seed', action='store_true', help='Use fixed random seed for reproducibility.')


    #for validation. all 
    modes = ["convert", "calcweight", "both"]
    valid_list_both = ["cloneout", "pointMutations", "VAF", "fasta", "maf"]
    valid_list_convert = ["cloneout", "pointMutations", "VAF", "fasta", "alt_weight"]
    valid_list_calcweight = ["maf"]

    args = parser.parse_args()
    if args.mode not in modes:
        raise ValueError(f"Invalid mode: {args.mode}")
    
    error_list=list()
    if args.mode == "both":
        for key in valid_list_both:
            if getattr(args, key) is None:
                error_list.append(f"Please specify {key} filepath.")
                continue
            if not Path(getattr(args, key)).exists():
                error_list.append(f'No such file: {getattr(args, key)}')
    elif args.mode == "convert":
        for key in valid_list_convert:
            if getattr(args, key) is None:
                error_list.append(f"Please specify {key} filepath.")
                continue
            if not Path(getattr(args, key)).exists():
                error_list.append(f'No such file: {getattr(args, key)}')
    else:
        for key in valid_list_calcweight:
            if getattr(args, key) is None:
                error_list.append(f"Please specify {key} filepath.")
                continue
            if not Path(getattr(args, key)).exists():
                error_list.append(f'No such file: {getattr(args, key)}')
    
    if len(error_list) > 0:
        for error in error_list:
            print(error)
            sys.exit(1)

    return args

if __name__ == '__main__':
    args = get_args()
    os.makedirs(args.output, exist_ok=True)
    #export data
    if args.mode in ["both", "calcweight"]:
        #print(convert_to_vcf(1, "../Output_2404/cloneout.txt", "../Output_2404/Mutations/pointMutations.txt", "../Output_2404/VAF/VAF.txt", 'VAF_primary', "./data/Homo_sapiens.GRCh38.dna.primary_assembly.fa", "./Demo01"))
        alt_weight_df = calc_weight_from_TCGA(
            args.maf,
            init_alt_weight_df
        )
        file_name, _ = os.path.splitext(args.maf)
        alt_weight_file = f"{args.output}/{file_name}_alt_weight.txt"
        alt_weight_df.to_csv(alt_weight_file, index=False, sep="\t")
    
    if args.mode in ["both", "convert"]:
        #when just convert mode, read intermediate alt_weight file
        if args.mode == "convert":
            alt_weight_df = pd.read_csv(args.alt_weight, sep="\t")

        new_df = convert_to_vcf(
            time=args.time,
            time_filter_logic=args.time_filter_logic,
            cloneout_filepath=args.cloneout,
            pointMutations_filepath=args.pointMutations,
            VAF_filepath=args.VAF,
            VAF_column=args.VAF_column,
            fasta_filepath=args.fasta,
            output_path=args.output,
            alt_weight_df=alt_weight_df,
            seed=args.seed
        )
