import pandas as pd


def combine_output_csvs(fnames, out_fname):
    beeg_df = pd.DataFrame()
    for fname in fnames:
        try:
            df = pd.read_csv(fname)
            df["sample"] = fname.split("/")[-1].split("_")[0]
            beeg_df = pd.concat([beeg_df, df])
        except:
            print(
                f"Warning: Failed appending summary dataframe to batch csv: {fname}")
            continue
    beeg_df.to_csv(out_fname)
    return f"Combined analytical output from all samples in batch saved to {out_fname}"
