import os
import csv

def output_to_file(output_path, output_data):
    """
    Output the data to a file.

    Args:
        output_path: The path to the output file
        output_data: The data to output
    """
    if "//" in output_path:
        output_path = output_path.replace("//", "/")

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        with open(output_path, "w") as f:
            # If a cell contains a comma, replace only that cell with a value enclosed in quotes
            output_data = output_data.apply(
                lambda col: col.map(lambda x: f'"{str(x)}"' if "," in str(x) else x)
            )

            output_data.to_csv(f, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    except IOError as e:
        raise ValueError(f"Error: Failed to write to file: {e}")
