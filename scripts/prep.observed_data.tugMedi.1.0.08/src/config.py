import toml
import re
from collections import defaultdict


class Config:
    def __init__(self, config_file):
        try:
            self.config = toml.load(config_file)
        except toml.TomlDecodeError as e:
            raise ValueError(
                f"The configuration file {config_file} contains parts that are not in TOML format: {e}"
            )

        self.values_dict = self.init_values_dict()

    def init_values_dict(self):
        values_dict = defaultdict(lambda: None)
        for section, items in self.config.items():
            for key, value in items.items():
                # If the parameter value is a list, handle it accordingly
                values_list = [item.strip() for item in str(value).strip("[] ").split(",")]
                if len(values_list) > 1:
                    values_dict[key] = values_list
                else:
                    values_dict[key] = value

        sample_list_file = values_dict.get("sample_file", None)
        if sample_list_file is not None:
            values_dict = self.load_sample_file(sample_list_file, values_dict)

        return values_dict

    def load_sample_file(self, sample_list_file, values_dict):
        tumor_sample_list = []
        normal_sample_list = []

        with open(sample_list_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tumor_sample_list.append(line.split("\t")[0].strip(' "\n'))
                if len(line.strip().split("\t")) > 1:
                    normal_sample_list.append(line.split("\t")[1].strip(' "\n'))

        dict_t_n = {tumor: normal for tumor, normal in zip(tumor_sample_list, normal_sample_list)}
        values_dict["tumor_normal_pair"] = dict_t_n

        if len(tumor_sample_list) > 0:
            values_dict["tumor_sample_labels"] = tumor_sample_list
            if len(normal_sample_list) > 0:
                values_dict["normal_sample_labels"] = normal_sample_list
            else:
                if "normal_sample_labels" in values_dict:
                    del values_dict["normal_sample_labels"]

        return values_dict

    def get_values_dict(self):
        return self.values_dict

    def get_index(self, key):
        x = self.values_dict.get(key, None)
        if x is None:
            return None
        elif isinstance(x, int):
            return x - 1
        else:
            return None

    def get_keys_from_section(self, section):
        try:
            return list(self.config[section].keys())
        except KeyError:
            return []

    def check(self, mode):
        if "tumor_sample_labels" in self.values_dict:
            if "sample_id_column_index" not in self.values_dict:
                raise ValueError("sample_id_column_index is not set")

        if mode == "all":
            self.check_tumor_specific_conditions(mode)
            self.check_vaf_conditions()
            self.check_driver_genes_conditions(mode)
            self.check_passenger_genes_conditions(mode)

        if mode == "get-tumor-specific":
            self.check_tumor_specific_conditions(mode)
        elif mode == "get-vaf":
            self.check_vaf_conditions()
        elif mode == "get-driver-genes":
            self.check_driver_genes_conditions(mode)
        elif mode == "get-passenger-genes":
            self.check_passenger_genes_conditions(mode)

    def check_tumor_specific_conditions(self, mode):
        # [1] In get-tumor-specific mode, all of the following A, B, C must be satisfied
        # A: When tumor_sample_labels is set
        # A1: sample_id_column_index must be set
        # B: When both tumor_sample_labels and normal_sample_labels are set
        # B1: The same number of samples must be set
        # B2: chromosome_column_index, position_column_index, t_ref_column_index, t_alt_column_index must be set
        # C: When normal_sample_labels is not set, (n_alt1_column_index or n_alt2_column_index) and n_alt_count_column_index must be set

        if (
            "normal_sample_labels" in self.values_dict
            and "tumor_sample_labels" in self.values_dict
        ):
            if len(self.values_dict["tumor_sample_labels"]) != len(
                self.values_dict["normal_sample_labels"]
            ):
                raise ValueError(
                    "The number of samples in tumor_sample_labels and normal_sample_labels do not match"
                )

            required_columns = [
                "chromosome_column_index",
                "position_column_index",
                "t_ref_column_index",
                "t_alt_column_index",
            ]
            for column in required_columns:
                if column not in self.values_dict:
                    raise ValueError(f"{column} is not set")

        if (
            mode == "get-tumor-specific"
            and "normal_sample_labels" not in self.values_dict
        ):
            if (
                "n_alt1_column_index" not in self.values_dict
                and "n_alt2_column_index" not in self.values_dict
                and "n_alt_count_column_index" not in self.values_dict
            ):
                raise ValueError(
                    "n_alt1_column_index, n_alt2_column_index, n_alt_count_column_index are not set"
                )

        if (
            mode == "get-tumor-specific"
            and "normal_sample_labels" not in self.values_dict
            and "n_alt_count_column_index" not in self.values_dict
        ):
            required_columns = [
                "chromosome_column_index",
                "position_column_index",
                "t_ref_column_index",
                "t_alt_column_index",
            ]
            for column in required_columns:
                if column not in self.values_dict:
                    raise ValueError(f"{column} is not set")

    def check_vaf_conditions(self):
        # [2] In get-vaf mode, all of the following A, B must be satisfied
        # A: When rewrite_vaf_mode = force, (t_ref_count_column_index or t_total_count_column_index) and t_alt_count_column_index must be set
        # B: When rewrite_vaf_mode = default, vaf_column_index or (ref_count_column_index or total_count_column_index) and alt_count_column_index must be set

        rewrite_vaf_mode = self.values_dict.get("rewrite_vaf_mode", "default")

        if rewrite_vaf_mode == "force":
            if ("t_ref_count_column_index" not in self.values_dict) and (
                "t_total_count_column_index" not in self.values_dict
            ):
                raise ValueError(
                    "t_ref_count_column_index or t_total_count_column_index is not set"
                )
            if "t_alt_count_column_index" not in self.values_dict:
                raise ValueError("t_alt_count_column_index is not set")

        elif rewrite_vaf_mode == "default":
            if "vaf_column_index" not in self.values_dict:
                if ("t_ref_count_column_index" not in self.values_dict) and (
                    "t_total_count_column_index" not in self.values_dict
                ):
                    raise ValueError(
                        "t_ref_count_column_index or t_total_count_column_index is not set"
                    )
                if "t_alt_count_column_index" not in self.values_dict:
                    raise ValueError("t_alt_count_column_index is not set")

    def check_columns_in_filter(self, filter, columns_section_keys, mode):
        variables = re.findall(r"\$(\w+)", filter)
        column_mapping = {
            "SampleID": "sample_id_column_index",
            "HugoSymbol": "hugo_symbol_column_index",
            "Chr": "chromosome_column_index",
            "Pos": "position_column_index",
            "Ref": "t_ref_column_index",
            "Alt": "t_alt_column_index",
            "VAF": "vaf_column_index",
            "TotalCount": "t_total_count_column_index",
            "ACount": "t_alt_count_column_index",
            "RCount": "t_ref_count_column_index",
            "Filter": "filter_column_index",
        }

        for var in variables:
            column_index_name = column_mapping.get(var)
            if not column_index_name:
                raise ValueError(f"${var} is an invalid variable name")

            if column_index_name not in columns_section_keys:
                if not mode == "all" or column_index_name != "vaf_column_index":
                    raise ValueError(f"{column_index_name} is not set")

    def check_driver_genes_conditions(self, mode):
        # [3] In get-driver-genes mode, all of the following A, B must be satisfied
        # A: driver_gene_include_filter or driver_gene_exclude_filter must be set
        # B: The column numbers corresponding to the variables appearing in driver_gene_include_filter and driver_gene_exclude_filter must be set in the Columns section of the config
        # C: When duplicate_gene_mode = vaf_max or vaf_min, vaf_column_index, hugo_symbol_column_index must be set

        if not (
            self.values_dict.get("driver_gene_include_filter")
            or self.values_dict.get("driver_gene_exclude_filter")
        ):
            raise ValueError(
                "driver_gene_include_filter or driver_gene_exclude_filter is not set"
            )

        if self.values_dict.get("duplicate_gene_mode", "all") in ["vaf_max", "vaf_min"]:
            if not self.values_dict.get("hugo_symbol_column_index"):
                raise ValueError("hugo_symbol_column_index is not set")
            if mode == "get-driver-genes" and not (
                self.values_dict.get("vaf_column_index")
            ):
                raise ValueError("vaf_column_index is not set")

        filter1 = self.values_dict.get("driver_gene_include_filter", "")
        filter2 = self.values_dict.get("driver_gene_exclude_filter", "")
        columns_section_keys = self.get_keys_from_section("Columns")
        self.check_columns_in_filter(
            filter1 + filter2, columns_section_keys, mode
        )

    def check_passenger_genes_conditions(self, mode):
        # [4] In get-passenger-genes mode, all of the following A, B must be satisfied
        # A: passenger_gene_include_filter or passenger_gene_exclude_filter must be set
        # B: The column numbers corresponding to the variables appearing in passenger_gene_include_filter and passenger_gene_exclude_filter must be set in the Columns section of the config

        if not (
            self.values_dict.get("passenger_gene_include_filter")
            or self.values_dict.get("passenger_gene_exclude_filter")
        ):
            raise ValueError(
                "passenger_gene_include_filter or passenger_gene_exclude_filter is not set"
            )

        if self.values_dict.get("passenger_gene_name_format"):
            if not self.values_dict.get("hugo_symbol_column_index"):
                raise ValueError("hugo_symbol_column_index is not set")

        filter1 = self.values_dict.get("passenger_gene_include_filter", "")
        filter2 = self.values_dict.get("passenger_gene_exclude_filter", "")
        filter3 = self.values_dict.get("passenger_gene_name_format", "")
        columns_section_keys = self.get_keys_from_section("Columns")
        self.check_columns_in_filter(
            filter1 + filter2 + filter3, columns_section_keys, mode
        )

    def set_output_dir(self, output_dir):
        self.values_dict["output_dir"] = output_dir
