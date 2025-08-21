import pandas as pd
import os

from app.utils.shell_cmds import loginfo


def depth_csvs(fnames, out_fname):
    beeg_df = pd.DataFrame()
    for fname in fnames:
        try:
            df = pd.read_csv(fname, index_col=0)
            df["sample"] = fname.split("/")[-1].split("_")[0]
            beeg_df = pd.concat([beeg_df, df])
        except:
            print(
                f"Warning: Failed appending summary dataframe to batch csv: {fname}")
            continue
    beeg_df.to_csv(out_fname.replace(".csv", "_depth.csv"))


def cov_csvs(fnames, out_fname):
    beeg_df = pd.DataFrame()  # TODO This really should be one function.
    for fname in fnames:
        try:
            df = pd.read_csv(fname.replace("_depth", "_coverage"), index_col=0)
            df["sample"] = fname.split("/")[-1].split("_")[0]
            beeg_df = pd.concat([beeg_df, df])
        except:
            print(
                f"Warning: Failed appending summary dataframe to batch csv: {fname}")
            continue
    beeg_df = beeg_df.reindex(sorted(beeg_df.columns), axis=1)
    beeg_df.to_csv(out_fname.replace(".csv", "_coverage.csv"))


def combine_output_csvs(fnames, out_fname):
    depth_csvs(fnames, out_fname)
    cov_csvs(fnames, out_fname)
    return f"Combined analytical output from all samples in batch saved to {out_fname}"


def combine_output_from_endpoint(payload):
    loginfo(
        f"Looking for depth and coverage csvs in {payload['DataFolder']}. This may take a moment while I scan the directory recursively.")
    fnames = [os.path.join(dp, f) for dp, _, fn in os.walk(
        os.path.expanduser(payload['DataFolder'])) for f in fn if "depth.csv" in f]
    out_fname = f"{payload['DataFolder']}/{payload['DataFolder'].split('/')[-1]}.csv"
    depth_csvs(fnames, out_fname)
    cov_csvs(fnames, out_fname)
    loginfo(
        f"Combined analytical output from all samples in batch saved to {out_fname}")
    return f"Combined analytical output from all samples in batch saved to {out_fname}"
