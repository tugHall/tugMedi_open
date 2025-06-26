import os
import pandas as pd
from file_helper import output_to_file


def get_tumor_specific_data(dict_t_samples, dict_n_samples, config):
    """
    Get the tumor specific data.

    Args:
        dict_t_samples: DataFrame of tumor samples
        dict_n_samples: DataFrame of normal samples
        config: Configuration object containing settings

    Returns:
        pd.DataFrame: The tumor specific data
    """
    t_sample_df_result_list = []

    for t_sample_id in dict_t_samples.keys():
        t_sample_df = dict_t_samples[t_sample_id]
        n_sample_df = dict_n_samples.get(t_sample_id, None)

        if n_sample_df is not None:
            chromosome_index = config.get_index("chromosome_column_index")
            position_index = config.get_index("position_column_index")
            ref_index = config.get_index("t_ref_column_index")
            alt_index = config.get_index("t_alt_column_index")

            t_sample_df["_tkey_"] = t_sample_df.apply(
                lambda row: f"{row.iloc[chromosome_index]}_{row.iloc[position_index]}_{row.iloc[ref_index]}_{row.iloc[alt_index]}",
                axis=1,
            )
            n_sample_df["_nkey1_"] = n_sample_df.apply(
                lambda row: f"{row['chromosome']}_{row['position']}_{row['ref']}_{row['n_alt']}",
                axis=1,
            )
            n_sample_df["_nkey2_"] = n_sample_df.apply(
                lambda row: f"{row['chromosome']}_{row['position']}_{row['ref']}_{row['n_alt1']}",
                axis=1,
            )
            n_sample_df["_nkey3_"] = n_sample_df.apply(
                lambda row: f"{row['chromosome']}_{row['position']}_{row['ref']}_{row['n_alt2']}",
                axis=1,
            )

            n_keys = set(
                n_sample_df["_nkey1_"].tolist()
                + n_sample_df["_nkey2_"].tolist()
                + n_sample_df["_nkey3_"].tolist()
            )

            t_sample_df_result = t_sample_df[~t_sample_df["_tkey_"].isin(n_keys)]
            t_sample_df_result = t_sample_df_result.drop(columns=["_tkey_"])
            t_sample_df_result_list.append(t_sample_df_result)
        else:
            t_sample_df_result_list.append(t_sample_df)

    return pd.concat(t_sample_df_result_list).sort_index()


def get_tumor_specific_data_main(
    args, config, log_file, dict_t_samples, dict_n_samples, basename
):
    """
    Main function for getting tumor specific data.

    Args:
        args: Command line arguments
        config: Configuration object containing settings
        log_file: File object for logging
        dict_t_samples: DataFrame of tumor samples
        dict_n_samples: DataFrame of normal samples
        basename: Base name for output files

    Output:
        tumor-specific-data file
    """
    print("get-tumor-specific mode start")
    log_file.write("get-tumor-specific mode start\n")
    output_path = os.path.join(args.output, f"{basename}_tumor_specific.txt")
    target_data = get_tumor_specific_data(dict_t_samples, dict_n_samples, config)
    output_to_file(output_path, target_data)
    print("get-tumor-specific mode passed")
    log_file.write("get-tumor-specific mode passed\n")
    return basename + "_tumor_specific", target_data
