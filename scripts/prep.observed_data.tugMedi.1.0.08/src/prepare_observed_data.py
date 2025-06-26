import os
import argparse
from datetime import datetime
import pandas as pd
from config import Config
from tumor_specific import get_tumor_specific_data_main
from vaf import vaf_mode_main
from gene_filter import (
    process_driver_genes_mode,
    process_passenger_genes_mode,
    driver_genes_mode_main,
    passenger_genes_mode_main,
)

# Functions related to reading sample files
def load_sample_data(dfs, config, all_samples):
    if all_samples == []:
        return pd.concat(dfs)
    else:
        sample_id_col_index = config.get_index("sample_id_column_index")

        # Read only rows containing elements in all_samples in the specified column
        return pd.concat(
            df[df.iloc[:, sample_id_col_index].isin(all_samples)]
            for df in dfs
        )


def chk_sample_names(dict_samples, tumor_samples):
    # Get sample names from t_samples
    t_sample_names = list(dict_samples.keys())

    # Check if the sample names in tumor_samples are included in t_samples
    for tumor_sample in tumor_samples:
        if tumor_sample not in t_sample_names:
            raise ValueError(
                f"Error: Sample {tumor_sample} is not included in the sample file"
            )

    return True


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process gene mutation data.")
    parser.add_argument(
        "-i", "--input", required=True, help="Input gene mutation data file"
    )
    parser.add_argument("-c", "--config", required=True, help="Configuration file")
    parser.add_argument(
        "-m",
        "--mode",
        choices=[
            "all",
            "get-tumor-specific",
            "get-vaf",
            "get-driver-genes",
            "get-passenger-genes",
        ],
        default="all",
        help="Mode of operation",
    )
    parser.add_argument("-o", "--output", default="./output/", help="Output directory")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug mode"
    )
    return parser.parse_args()


def setup_log_file(args, basename):
    log_file_name = (
        f"prepare_observed_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    )
    log_file_path = f"{args.output}/{log_file_name}"
    # Create folder
    log_file_dir = os.path.dirname(log_file_path)
    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)
    log_file = open(log_file_path, "w")
    log_file.write(f"Start date and time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log_file.write(f"Sample file path: {args.input}\n")
    log_file.write(f"Config file path: {args.config}\n")
    log_file.write(
        f"Command content: python prepare_observed_data.py -i {args.input} -c {args.config} -m {args.mode} -o {args.output}\n"
    )
    return log_file


def load_config(args):
    print("Read config file")
    try:
        config = Config(args.config)
    except Exception as e:
        raise ValueError(f"Error: Failed to load configuration file: {e}")

    config.check(args.mode)
    return config


def set_tumor_samples(data, label_samples, config):
    dict_samples = {}

    if label_samples == []:
        dict_samples["_sample_id_"] = data
    else:
        sample_id_col_index = config.get_index("sample_id_column_index")
        for sample_id in label_samples:
            dict_samples[sample_id] = data[
                data.iloc[:, sample_id_col_index] == sample_id
            ]

    return dict_samples


def set_normal_sample(data, label_samples, config):
    dict_samples = {}

    # Get chromosome_column_index, position_column_index, t_ref_column_index, t_alt_column_index from the Columns section of the config file
    # If not set, set to None
    sample_id_index = config.get_index("sample_id_column_index")
    chromosome_index = config.get_index("chromosome_column_index")
    position_index = config.get_index("position_column_index")
    ref_index = config.get_index("t_ref_column_index")
    alt_index = config.get_index("t_alt_column_index")
    n_alt1_index = None
    n_alt2_index = None
    n_alt_count_index = None

    selected_columns = [
        sample_id_index,
        chromosome_index,
        position_index,
        ref_index,
        alt_index,
        n_alt1_index,
        n_alt2_index,
        n_alt_count_index,
    ]

    n_data = pd.DataFrame()
    for column_index in selected_columns:
        if column_index is not None:
            # Add columns to n_data from the first column
            n_data = pd.concat([n_data, data.iloc[:, column_index]], axis=1)
        else:
            # Add an empty DataFrame
            empty_data = pd.DataFrame(
                [""] * len(data), index=data.index, columns=["_n_col_"]
            )
            n_data = pd.concat([n_data, empty_data], axis=1)

    n_data.columns = [
        "sample_id",
        "chromosome",
        "position",
        "ref",
        "n_alt",
        "n_alt1",
        "n_alt2",
        "n_alt_count",
    ]

    for sample_id in label_samples:
        dict_samples[sample_id] = n_data[n_data.iloc[:, 0] == sample_id]

    return dict_samples


def set_normal_samples_parallel(data, config):
    dict_samples = {}

    # Get chromosome_column_index, position_column_index, t_ref_column_index, t_alt_column_index from the Columns section of the config file
    sample_id_index = config.get_index("sample_id_column_index")
    chromosome_index = config.get_index("chromosome_column_index")
    position_index = config.get_index("position_column_index")
    ref_index = config.get_index("t_ref_column_index")
    alt_index = None  # Set alt_index to None as we only want the Normal Sample's alt column
    n_alt1_index = config.get_index("n_alt1_column_index")
    n_alt2_index = config.get_index("n_alt2_column_index")
    n_alt_count_index = config.get_index("n_alt_count_column_index")

    selected_columns = [
        sample_id_index,
        chromosome_index,
        position_index,
        ref_index,
        alt_index,
        n_alt1_index,
        n_alt2_index,
        n_alt_count_index,
    ]

    n_data = pd.DataFrame()
    for column_index in selected_columns:
        if column_index is not None:
            # Add columns to n_data from the first column
            n_data = pd.concat([n_data, data.iloc[:, column_index]], axis=1)
        else:
            # Add an empty DataFrame
            empty_data = pd.DataFrame(
                [""] * len(data), index=data.index, columns=["_n_col_"]
            )
            n_data = pd.concat([n_data, empty_data], axis=1)

    n_data.columns = [
        "sample_id",
        "chromosome",
        "position",
        "ref",
        "n_alt",
        "n_alt1",
        "n_alt2",
        "n_alt_count",
    ]

    return n_data


def have_normal_sample_column(config):
    # Check if n_alt1_column_index, n_alt2_column_index, n_alt_count_column_index is set
    n_alt1_index = config.get_index("n_alt1_column_index")
    n_alt2_index = config.get_index("n_alt2_column_index")
    n_alt_count_index = config.get_index("n_alt_count_column_index")

    return (
        n_alt1_index is not None
        or n_alt2_index is not None
        or n_alt_count_index is not None
    )


def read_sample_files(args, config):
    """
    Read the sample file.

    Args:
        args: Command line arguments.
        config: Configuration object containing settings.

    Returns:
        dict_t_samples: DataFrame of tumor samples
        dict_n_samples: DataFrame of normal samples
    """
    print("Read sample file")
    tumor_sample_labels = config.get_values_dict().get("tumor_sample_labels", [])
    normal_sample_labels = config.get_values_dict().get("normal_sample_labels", [])

    # Read only rows containing the specified sample ID
    df = pd.read_csv(args.input, sep="\t", chunksize=100, comment="#")
    data = load_sample_data(df, config, tumor_sample_labels + normal_sample_labels)

    # For tumor samples, split data by sample ID
    dict_t_samples = set_tumor_samples(data, tumor_sample_labels, config)
    if not chk_sample_names(dict_t_samples, tumor_sample_labels):
        raise ValueError(
            f"Error: There are samples not included in the sample file"
        )

    # For normal samples, split data by sample ID
    dict_n_samples = {}
    if len(normal_sample_labels) > 0:
        # If normal_sample_labels is set in the config file
        # Read data with normal_sample_label as key
        dict_n_samples_tmp = set_normal_sample(data, normal_sample_labels, config)
        if not chk_sample_names(dict_n_samples_tmp, normal_sample_labels):
            raise ValueError(
                f"Error: There are samples not included in the sample file"
            )

        # Convert to a dictionary with tumor_sample_label as key
        for sample_id in dict_t_samples.keys():
            # Get the normal_sample_id corresponding to tumor_sample_id from config
            normal_sample_id = (
                config.get_values_dict().get("tumor_normal_pair", {}).get(sample_id)
            )
            dict_n_samples[sample_id] = dict_n_samples_tmp[normal_sample_id]

    elif have_normal_sample_column(config):
        # If there is data for Normal Sample in another column
        for sample_id in dict_t_samples.keys():
            data = dict_t_samples[sample_id]
            dict_n_samples[sample_id] = set_normal_samples_parallel(data, config)

    return dict_t_samples, dict_n_samples


def main():
    args = parse_arguments()
    basename = os.path.splitext(os.path.basename(args.input))[0]

    # Process log file
    log_file = setup_log_file(args, basename)

    # Load configuration file
    config = load_config(args)
    config.set_output_dir(args.output)

    # Read sample file
    dict_t_samples, dict_n_samples = read_sample_files(args, config)
    if dict_t_samples is None:
        return
    target_data = pd.concat(
        [t_sample for t_sample in dict_t_samples.values()]
    ).sort_index()

    # get-tumor-specific mode
    if args.mode in ["all", "get-tumor-specific"]:
        basename, target_data = get_tumor_specific_data_main(
            args, config, log_file, dict_t_samples, dict_n_samples, basename
        )

    # get-vaf mode
    if args.mode in ["all", "get-vaf"]:
        basename, target_data = vaf_mode_main(args, config, log_file, target_data, basename)

    # get-driver-genes mode
    if args.mode in ["all", "get-driver-genes"]:
        driver_genes_mode_main(args, config, log_file, target_data, basename)

    # get-passenger-genes mode
    if args.mode in ["all", "get-passenger-genes"]:
        passenger_genes_mode_main(args, config, log_file, target_data, basename)


if __name__ == "__main__":
    main()