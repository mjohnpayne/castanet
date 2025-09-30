import os
import random
import pandas as pd
from app.utils.utility_fns import read_fa
from app.utils.shell_cmds import loginfo


def parse_csv(fpath) -> pd.DataFrame:
    return pd.read_csv(fpath, index_col=False)


def make_csv(fpath) -> pd.DataFrame:
    fastas = read_fa(fpath)
    agg_headers, descriptions, seqs = [], [], []
    for fasta in fastas:
        agg_headers.append(fasta[0].split("_")[0].replace(">", ""))
        descriptions.append("_".join(fasta[0].split("_")[1:]))
        seqs.append(fasta[1])
    df = pd.DataFrame({"aggregate_to": agg_headers,
                      "description": descriptions, "sequence": seqs})
    return df


def input_checks(df) -> pd.DataFrame:
    aggregation_headers = df["aggregate_to"].unique().tolist()
    disallowed_chars = [" ", "/", "\\", ":",
                        "*", "?", "\"", "<", ">", "|", ",", "_"]
    for header in aggregation_headers:
        for char in disallowed_chars:
            if char in header:
                raise ValueError(f"Disallowed character '{char}' found in aggregate_to value '{header}'"
                                 f"Please remove or replace it, then try again.")
    return df


def generate_hash(df) -> pd.DataFrame:
    df["unique_label"] = df.apply(
        lambda x: f"{random.getrandbits(128)}", axis=1)
    return df


def generate_fasta(df) -> list:
    fasta = []
    for _, row in df.iterrows():
        fasta.append(
            [f">{row['aggregate_to']}_{row['unique_label']}", row['sequence']])
    return fasta


def save_output(df, fasta, out_fpath) -> None:
    with open(out_fpath, "w") as f:
        for header, seq in fasta:
            f.write(f"{header}\n{seq}\n")
    df.to_csv(f"{out_fpath.replace('.fasta', '.csv')}", index=False)


def main(in_file, out_file):
    '''Convert an input CSV or FASTA mapping reference description file to a Castanet-compatible RefStem'''
    loginfo(f"Converting mapping reference file: {in_file}")
    if in_file.endswith(".csv"):
        df = parse_csv(in_file)
    elif in_file.endswith(".fasta") or in_file.endswith(".fa"):
        df = make_csv(in_file)
    else:
        raise ValueError("Input file must be a .csv or .fasta/.fa file")

    df = input_checks(df)
    df = generate_hash(df)
    fasta = generate_fasta(df)
    save_output(df, fasta, out_file)
    complete_msg = f"Conversion complete! Output saved to: {out_file} and {out_file.replace('.fasta', '.csv')}"
    loginfo(complete_msg)
    return complete_msg


if __name__ == "__main__":
    in_file = "test_mapping_ref.csv"
    # in_file = "test_mapping_ref_converted.fasta"
    out_file = "test_mapping_ref_converted.fasta"

    _ = main(in_file, out_file)
