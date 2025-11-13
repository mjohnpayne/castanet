import os
import random
import pandas as pd
from app.utils.utility_fns import read_fa
from app.utils.shell_cmds import loginfo, logerr, stoperr


class MappingRefConverter:
    def __init__(self, payload, sneaky_mode=True):
        self.payload = payload
        self.default_aggregation_val = "unaggregated"
        if not sneaky_mode:
            self.in_file = payload["InFile"]
            self.out_file = payload["OutFile"]
        else:
            self.in_file = self.payload["RefStem"]
            self.out_file = f"{payload['SaveDir']}/{payload['ExpName']}/ref.fa"
            self.payload["RefStem"] = self.out_file
        self.sneaky_mode = sneaky_mode

    def make_csv(self) -> pd.DataFrame:
        try:
            fastas = read_fa(self.in_file)
        except Exception as e:
            raise ValueError(
                f"Error reading fasta file {self.in_file}. Please ensure it is in valid fasta format.") from e
        agg_headers, descriptions, seqs, organisms = [], [], [], []
        for fasta in fastas:
            if len(fasta[0].split("_")) < 2:
                logerr(f"Mapping reference aggregator {fasta[0]} has no underscores, so will not aggregate with any other references! Please refer to documentation. "
                       f"I'm setting this to '{self.default_aggregation_val}'.")
                agg_headers.append(self.default_aggregation_val)
                descriptions.append(fasta[0].replace(">", ""))
            else:
                agg_headers.append(fasta[0].split("_")[0].replace(">", ""))
                descriptions.append(
                    "_".join(fasta[0].split("_")[1:]).replace(",", ""))
            organisms.append(fasta[0].split(
                "_")[0].replace(">", "").split("-")[0])
            seqs.append(fasta[1])
        df = pd.DataFrame({"organism": organisms, "probetype": agg_headers,
                           "description": descriptions, "sequence": seqs})
        return df

    def input_checks(self, df) -> pd.DataFrame:
        aggregation_headers = df["probetype"].unique().tolist()
        disallowed_chars = [" ", "/", "\\", ":",
                            "*", "?", "\"", "<", ">", "|", ",", "_"]

        for header in aggregation_headers:
            for char in disallowed_chars:
                if char in header:
                    raise ValueError(f"Disallowed character '{char}' found in probetype value '{header}'"
                                     f"Please remove or replace it, then try again.")
        return df

    def generate_hash(self, df) -> pd.DataFrame:
        df["key"] = df.apply(
            lambda x: f"{random.getrandbits(128)}", axis=1)
        return df

    def generate_fasta(self, df) -> list:
        fasta = []
        for _, row in df.iterrows():
            fasta.append(
                [f">{row['probetype']}_{row['key']}", row['sequence']])
        return fasta

    def save_output(self, df, fasta) -> None:
        if not self.sneaky_mode:
            '''Save CSV'''
            df[["organism", "probetype", "description", "key"]].to_csv(
                f"{self.out_file}", index=False)
        else:
            if not os.path.exists(f"{self.payload['SaveDir']}/{self.payload['ExpName']}/"):
                os.makedirs(
                    f"{self.payload['SaveDir']}/{self.payload['ExpName']}/")
            self.payload["MappingRefTable"] = f"{self.payload['SaveDir']}/{self.payload['ExpName']}/MappingRefTable.csv"
            df = df.applymap(lambda s: s.lower() if type(s) == str else s)
            df[["organism", "probetype", "description", "key"]].to_csv(
                self.payload["MappingRefTable"], index=False)
            with open(self.out_file, "w") as f:
                for header, seq in fasta:
                    f.write(f"{header}\n{seq}\n")

    def validate_user_csv(self, df):
        if df.isnull().values.any():
            stoperr(f"Your input MappingRefTable has empty values in the probetype and/or description columns. "
                    f"Castanet can't proceed as it needs names for each target we map to. "
                    f"Please manually edit these, or re-generate the mapping reference with the /convert_mapping_ref/ function.")

    def join_seqs_to_df(self, df, users_df):
        out_df = users_df.copy()
        out_df["sequence"] = df["sequence"]
        return out_df

    def main(self):
        '''Convert an input CSV or FASTA mapping reference description file to a Castanet-compatible RefStem'''
        df = self.make_csv()
        if not self.sneaky_mode:
            loginfo(f"Converting mapping reference file: {self.in_file}")
            df = self.input_checks(df)
            df = self.generate_hash(df)
            fasta = ""

        if self.sneaky_mode:
            if os.path.isfile(self.payload["MappingRefTable"]):
                loginfo(
                    f'Parsing user supplied MappingRefTable: {self.payload["MappingRefTable"]}.')
                users_df = pd.read_csv(
                    self.payload["MappingRefTable"], index_col=None)
                self.validate_user_csv(users_df)
                df = self.join_seqs_to_df(df, users_df)
            else:
                loginfo(
                    f"No MappingRefTable supplied. Generating one automatically.")
                df = self.input_checks(df)
                df = self.generate_hash(df)
            fasta = self.generate_fasta(df)

        self.save_output(df, fasta)
        complete_msg = f"Conversion complete! Output saved to: {self.out_file}"
        if not self.sneaky_mode:
            loginfo(complete_msg)
            return complete_msg
        else:
            return self.payload
