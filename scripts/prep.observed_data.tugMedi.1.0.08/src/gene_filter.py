import os
import pandas as pd
from file_helper import output_to_file


def driver_genes_mode_main(args, config, log_file, target_data, basename):
    """
    Main function for processing driver genes mode.

    Args:
        args: Command line arguments.
        config: Configuration object containing settings.
        log_file: File object for logging.
        target_data (pd.DataFrame): The input data.
        basename (str): Base name for output files.

    Output:
        driver-gene-data file
    """
    print("get-driver-genes mode start")
    log_file.write("get-driver-genes mode start\n")
    output_path = os.path.join(args.output, f"{basename}.Rint.txt")
    driver_data = process_driver_genes_mode(target_data, config)
    if config.get_values_dict().get("output_format", None) == "part":
        driver_data = get_part_df(driver_data, config, rename_f=False)
    output_to_file(output_path, driver_data)
    print("get-driver-genes mode passed")
    log_file.write("get-driver-genes mode passed\n")


def process_driver_genes_mode(target_data, config):
    """
    Process the driver genes based on the given configuration.

    Args:
        target_data (pd.DataFrame): The input data.
        config: Configuration object containing settings.

    Returns:
        pd.DataFrame: The filtered driver genes.
    """
    renamed_part_df = get_part_df(target_data, config, rename_f=True)
    din_filter = config.get_values_dict().get("driver_gene_include_filter", None)
    dex_filter = config.get_values_dict().get("driver_gene_exclude_filter", None)

    if din_filter is None and dex_filter is None:
        return []

    din_condition = gene_filter(renamed_part_df, din_filter, dex_filter, config)
    driver_df = get_driver(target_data, din_condition, config)

    return driver_df


def get_driver(target_data, din_condition, config):
    """ 
    Extract driver genes based on the given condition and configuration.

    Args:
        target_data (pd.DataFrame): The input data.
        din_condition (pd.Series): The condition for filtering driver genes.
        config: Configuration object containing settings.

    Returns:
        pd.DataFrame: The filtered driver genes.
    """
    df = target_data.reset_index(drop=True)  # Reset index to ensure alignment
    driver_df = df[din_condition]

    vaf_index = config.get_index("vaf_column_index")
    hugo_symbol_index = config.get_index("hugo_symbol_column_index")
    sample_id_index = config.get_index("sample_id_column_index")
    vaf_column_name = driver_df.columns[vaf_index]
    hugo_symbol_column_name = driver_df.columns[hugo_symbol_index]

    if config.get_values_dict().get("duplicate_gene_mode", None) == "vaf_max":
        # When duplicate_gene_mode is vaf_max, extract the row with the maximum VAF value
        # Extract the row with the maximum value in the column indicated by self.vaf_index

        # If sample_id_column_index is set, also consider sample_id
        if sample_id_index is not None:
            sample_id_column_name = driver_df.columns[sample_id_index]
            driver_df = driver_df.sort_values(
                by=vaf_column_name, ascending=False
            ).drop_duplicates(
                subset=[hugo_symbol_column_name, sample_id_column_name], keep="first"
            )
        else:
            driver_df = driver_df.sort_values(
                by=vaf_column_name, ascending=False
            ).drop_duplicates(subset=hugo_symbol_column_name, keep="first")

    elif config.get_values_dict().get("duplicate_gene_mode", None) == "vaf_min":
        # When duplicate_gene_mode is vaf_min, extract the row with the minimum VAF value
        # Extract the row with the minimum value in the column indicated by self.vaf_index

        # If sample_id_column_index is set, also consider sample_id
        if sample_id_index is not None:
            sample_id_column_name = driver_df.columns[sample_id_index]
            driver_df = driver_df.sort_values(
                by=vaf_column_name, ascending=True
            ).drop_duplicates(
                subset=[hugo_symbol_column_name, sample_id_column_name], keep="first"
            )
        else:
            driver_df = driver_df.sort_values(
                by=vaf_column_name, ascending=True
            ).drop_duplicates(subset=hugo_symbol_column_name, keep="first")

    return driver_df


def passenger_genes_mode_main(args, config, log_file, target_data, basename):
    """
    Main function for processing passenger genes mode.

    Args:
        args: Command line arguments.
        config: Configuration object containing settings.
        log_file: File object for logging.
        target_data (pd.DataFrame): The input data.
        basename (str): Base name for output files.

    Output:
        passenger-gene-data file
    """
    print("get-passenger-genes mode start")
    log_file.write("get-passenger-genes mode start\n")
    passenger_data = process_passenger_genes_mode(target_data, config)
    driver_data = process_driver_genes_mode(target_data, config)

    # Limit to records in passenger_data that are not included in driver_data
    if len(driver_data) > 0:
        passenger_data = passenger_data[
            ~passenger_data.isin(driver_data.to_dict(orient="list")).all(axis=1)
        ]

    # If a gene name format is set, change the gene name
    p_gene_name_format = config.get_values_dict().get(
        "passenger_gene_name_format", None
    )
    if p_gene_name_format is not None:
        passenger_data = rename_hugo_symbol(passenger_data, p_gene_name_format, config)
    else:
        if config.get_values_dict().get("output_format", None) == "part":
            passenger_data = get_part_df(passenger_data, config, rename_f=False)

    output_path = os.path.join(args.output, f"{basename}.Rother.txt")
    output_to_file(output_path, passenger_data)
    print("get-passenger-genes mode passed")
    log_file.write("get-passenger-genes mode passed\n")


def process_passenger_genes_mode(target_data, config):
    """
    Process the passenger genes based on the given configuration.

    Args:
        target_data (pd.DataFrame): The input data.
        config: Configuration object containing settings.

    Returns:
        pd.DataFrame: The filtered passenger genes.
    """
    renamed_part_df = get_part_df(target_data, config, rename_f=True)

    # Load passenger gene filters
    pin_filter = config.get_values_dict().get("passenger_gene_include_filter", None)
    pex_filter = config.get_values_dict().get("passenger_gene_exclude_filter", None)
    pin_condition = gene_filter(renamed_part_df, pin_filter, pex_filter, config)

    # Reset index to ensure alignment
    target_data = target_data.reset_index(drop=True)
    pin_condition = pin_condition.reset_index(drop=True)
    passenger_df = target_data[pin_condition]

    return passenger_df


def get_part_df(sample_df, config, rename_f=False):
    """
    Extract and optionally rename parts of the DataFrame based on configuration.

    Args:
        sample_df (pd.DataFrame): The input data.
        config: Configuration object containing settings.
        rename_f (bool): Flag to determine if renaming is required.

    Returns:
        pd.DataFrame: The part of the DataFrame, renamed if specified.
    """
    part_df = pd.DataFrame()
    renamed_part_df = pd.DataFrame()

    indices = {
        "HugoSymbol": config.get_index("hugo_symbol_column_index"),
        "SampleID": config.get_index("sample_id_column_index"),
        "Chr": config.get_index("chromosome_column_index"),
        "Pos": config.get_index("position_column_index"),
        "Ref": config.get_index("ref_column_index"),
        "Alt": config.get_index("alt_column_index"),
        "TotalCount": config.get_index("total_count_column_index"),
        "RCount": config.get_index("ref_count_column_index"),
        "ACount": config.get_index("alt_count_column_index"),
        "VAF": config.get_index("vaf_column_index"),
        "Filter": config.get_index("filter_column_index"),
    }
    # Convert chromosome data. To handle it as a sequence, treat X, Y, etc. as 23, 24
    transform_chr = lambda x: int(
        str(x)
        .replace("chr", "")
        .replace("x", "23")
        .replace("X", "23")
        .replace("y", "24")
        .replace("Y", "24")
        .replace("mt", "25")
        .replace("MT", "25")
        .replace("m", "25")
        .replace("M", "25")
    )

    for name, index in indices.items():
        transform_func = transform_chr if name == "Chr" else None
        part_df, renamed_part_df = process_column(
            sample_df, part_df, renamed_part_df, index, name, transform_func
        )

    other_indices = config.get_values_dict().get("other_column_index", [])
    for i, index in enumerate(other_indices):
        part_df, renamed_part_df = process_column(
            sample_df, part_df, renamed_part_df, int(index) - 1, f"Other{i}"
        )

    return renamed_part_df if rename_f else part_df


def process_column(
    sample_df, part_df, renamed_part_df, index, name, transform_func=None
):
    """
    Process a single column of the DataFrame, optionally transforming it.
    """
    if index is not None:
        pcol = sample_df.iloc[:, index]
        part_df = pd.concat([part_df, pcol], axis=1)
        renamed_pcol = pcol.rename(name)
        if transform_func:
            renamed_pcol = renamed_pcol.apply(transform_func)
        renamed_part_df = pd.concat([renamed_part_df, renamed_pcol], axis=1)
    return part_df, renamed_part_df


def gene_filter(renamed_part_df, in_filter, ex_filter, config):
    """
    Apply gene filters to the DataFrame and return the condition for filtering.

    Args:
        renamed_part_df (pd.DataFrame): The DataFrame with renamed columns.
        in_filter (str): Include filter condition.
        ex_filter (str): Exclude filter condition.
        config: Configuration object containing settings.

    Returns:
        pd.Series: Boolean series indicating which rows match the filter conditions.
    """
    # Save renamed_part_df to "tmp_df.csv" in output_dir
    output_dir = config.get_values_dict().get("output_dir", None)
    if output_dir is not None:
        renamed_part_df.to_csv(
            os.path.join(output_dir, "tmp_df.csv"), index=False, header=False
        )

    # Store the combination of column names and column numbers of rename_part_df in a dictionary
    column_name_to_index = {col: i + 1 for i, col in enumerate(renamed_part_df.columns)}

    in_gene_df = pd.DataFrame()
    if in_filter is not None and len(in_filter) > 0:
        in_filter_num = convert_awk_column_number(in_filter, column_name_to_index)

        # Execute the awk command on tmp_df.csv and save the result to tmp_driver.csv
        os.system(
            f"awk -F',' '{in_filter_num}' {os.path.join(output_dir, 'tmp_df.csv')} > {os.path.join(output_dir, 'tmp_in_gene.csv')}"
        )
        # Load tmp_driver.csv into driver_df. Treat all as text
        if os.path.getsize(os.path.join(output_dir, "tmp_in_gene.csv")) > 0:
            in_gene_df = pd.read_csv(
                os.path.join(output_dir, "tmp_in_gene.csv"), header=None, dtype=str
            )

    ex_gene_df = pd.DataFrame()
    if ex_filter is not None and len(ex_filter) > 0:
        ex_filter_num = convert_awk_column_number(ex_filter, column_name_to_index)

        # Execute the awk command on tmp_df.csv and save the result to tmp_gene.csv
        os.system(
            f"awk -F',' '{ex_filter_num}' {os.path.join(output_dir, 'tmp_df.csv')} > {os.path.join(output_dir, 'tmp_ex_gene.csv')}"
        )
        # Load tmp_driver.csv into driver_df. Treat all as text
        ex_gene_df = pd.DataFrame()
        if os.path.getsize(os.path.join(output_dir, "tmp_ex_gene.csv")) > 0:
            ex_gene_df = pd.read_csv(
                os.path.join(output_dir, "tmp_ex_gene.csv"), header=None, dtype=str
            )
            # Exclude ex_gene_df from in_gene_df
            os.system(
                f"awk 'NR==FNR{{b[$0];next}} !($0 in b)' {os.path.join(output_dir, 'tmp_ex_gene.csv')} {os.path.join(output_dir, 'tmp_in_gene.csv')} > {os.path.join(output_dir, 'tmp_in_ex_gene.csv')}"
            )
            in_gene_df = pd.read_csv(
                os.path.join(output_dir, "tmp_in_ex_gene.csv"), header=None, dtype=str
            )

    gene_df = pd.read_csv(
        os.path.join(output_dir, "tmp_df.csv"), header=None, dtype=str
    )

    # Reset index to ensure alignment
    gene_df = gene_df.reset_index(drop=True)
    in_gene_df = in_gene_df.reset_index(drop=True)

    in_condition = pd.Series(False, index=gene_df.index)

    in_gene_df_set = set(["_".join(str(item) for item in row) for _, row in in_gene_df.iterrows()])
    # Set records in gene_df that match records in in_gene_df to True
    for index, row in gene_df.iterrows():
        # Create a string that concatenates all of row
        row_str = "_".join(str(item) for item in row)
        if row_str in in_gene_df_set:
            in_condition[index] = True

    # Delete tmp files
    if os.path.exists(os.path.join(output_dir, "tmp_df.csv")):
        os.remove(os.path.join(output_dir, "tmp_df.csv"))
    if os.path.exists(os.path.join(output_dir, "tmp_in_gene.csv")):
        os.remove(os.path.join(output_dir, "tmp_in_gene.csv"))
    if os.path.exists(os.path.join(output_dir, "tmp_ex_gene.csv")):
        os.remove(os.path.join(output_dir, "tmp_ex_gene.csv"))
    if os.path.exists(os.path.join(output_dir, "tmp_in_ex_gene.csv")):
        os.remove(os.path.join(output_dir, "tmp_in_ex_gene.csv"))

    return in_condition


def rename_hugo_symbol(df, p_gene_name_format, config):
    """
    Rename the HugoSymbol column based on the given format.
    """
    renamed_part_df = get_part_df(df, config, rename_f=True)
    part_df = get_part_df(df, config, rename_f=False)

    new_gene_name_series = renamed_part_df.apply(
        lambda row: generate_column_template(row, p_gene_name_format), axis=1
    )

    # If the HugoSymbol column exists, replace that column with new_gene_name_series.
    hugo_symbol_column_index = config.get_index("hugo_symbol_column_index")
    if hugo_symbol_column_index is not None:
        hugo_symbol_column_name = df.columns[hugo_symbol_column_index]
        df.loc[:, hugo_symbol_column_name] = new_gene_name_series
        part_df.loc[:, hugo_symbol_column_name] = new_gene_name_series
    else:
        # If the HugoSymbol column does not exist in target_data, add the "HugoSymbol" column from part_df to target_data.
        df["HugoSymbol"] = new_gene_name_series
        part_df["HugoSymbol"] = new_gene_name_series

    if config.get_values_dict().get("output_format", None) == "part":
        return part_df
    return df


def generate_column_template(row, template):
    """
    Replace variables in the template string with corresponding column values.
    """
    result = template
    result = result.strip().strip('"')
    for column in row.index:
        result = result.replace(f"${column}", str(row[column]))
    return result


def convert_awk_column_number(awk_condition, column_name_to_index):
    """
    Replace $column_name with $column_number in the awk condition.
    """
    for column_name, index in column_name_to_index.items():
        awk_condition = awk_condition.replace(f"${column_name}", f"${index}")
    return awk_condition
