import os
from file_helper import output_to_file


def calculate_vaf(row, config):
    t_ref_count_column_index = config.get_index("t_ref_count_column_index")
    t_alt_count_column_index = config.get_index("t_alt_count_column_index")
    t_total_count_column_index = config.get_index("t_total_count_column_index")

    # データを数値型に変換
    t_alt_count = float(row.iloc[t_alt_count_column_index])
    if t_ref_count_column_index is not None:
        t_ref_count = float(row.iloc[t_ref_count_column_index])
        vaf = t_alt_count / (t_ref_count + t_alt_count)
    else:
        t_total_count = float(row.iloc[t_total_count_column_index])
        vaf = t_alt_count / t_total_count

    return vaf


def process_vaf_mode(target_data, config):
    # If the vaf column index is specified,
    # If rewrite_vaf_mode is force, overwrite the vaf column
    # If rewrite_vaf_mode is default, leave the vaf column as is
    opt_vaf = config.get_values_dict().get("rewrite_vaf_mode", "default")
    vaf_column_index = config.get_index("vaf_column_index")

    if vaf_column_index is not None:
        # If the vaf column exists
        if opt_vaf == "force":
            # Overwrite the vaf column
            target_data.iloc[:, vaf_column_index] = target_data.apply(
                lambda row: calculate_vaf(row, config), axis=1
            )
        elif opt_vaf == "default":
            # Leave the vaf column as is
            pass
    else:
        # If the vaf column does not exist, add a vaf column at the end
        target_data["vaf"] = target_data.apply(
            lambda row: calculate_vaf(row, config), axis=1
        )

        # Update the vaf_column_index in the config
        vaf_column_index = target_data.columns.get_loc("vaf")
        config.get_values_dict()["vaf_column_index"] = vaf_column_index + 1

    return target_data


def vaf_mode_main(args, config, log_file, target_data, basename):
    """
    Main function for the VAF mode.

    Args:
        args: Command line arguments
        config: Configuration object containing settings
        log_file: File object for logging
        target_data: DataFrame containing the target data
        basename: Base name for output files

    Output:
        vaf-data file
    """
    print("get-vaf mode start")
    log_file.write("get-vaf mode start\n")
    output_path = os.path.join(args.output, f"{basename}_vaf.txt")
    target_data = process_vaf_mode(target_data, config)
    output_to_file(output_path, target_data)
    print("get-vaf mode passed")
    log_file.write("get-vaf mode passed\n")
    return basename + "_vaf", target_data
