import os
import pandas as pd
from app.utils.shell_cmds import stoperr, logerr


def make_exp_dir(ExpNameAndRoot):
    '''Checks if dir exists; creates if not'''
    if not os.path.exists(ExpNameAndRoot):
        os.makedirs(ExpNameAndRoot)


def get_gene_orgid(target_id):
    '''Find gene and orgid at specific ref; return (gene, orgid).'''
    parts = target_id.split('_')
    orgid = parts[-1] if parts[0].lower().startswith('bact') else parts[0]
    return (parts[0], orgid)


def read_fa(fpath):
    '''Read fasta file to list of lists in format [[>name, seq], [...]]'''
    seqs = []
    try:
        with open(fpath, "r") as f:
            for l in f:
                if l[0] == ">":
                    seqs.append(f"?{l}")
                else:
                    seqs.append(l.replace("\n", ""))
    except:
        '''Current release of ViralConsensus will sometimes dump binary into the output file - this is a temp fix'''
        print(
            f"*** WARNING: Encoding on your input fasta file ({fpath}) is messed up. Attempting to rescue it...")
        seqs = []  # RESET SEQS as it doesn't get wiped by transition to exception
        bases = ["A", "T", "C", "G", "N", "-"]

        with open(fpath, "r", encoding="ISO-8859-1") as f:
            for l in f.readlines():
                if l[0] not in bases:
                    seqs.append(f"?{l}")
                else:
                    seqs.append(l)

        for i in range(len(seqs)):
            if seqs[i][0] == "?":
                continue
            else:
                seqs[i] = "".join([j for j in seqs[i] if j in bases])

    seqs_split = [i.split("\n") for i in "".join(seqs).split("?")]
    return [i for i in seqs_split if not i == [""]]


def save_fa(fpath, pat):
    with open(fpath, "w") as f:
        f.write(pat)


def trim_long_fpaths(key, max_len=100):
    '''Curtail very long probe names that can't be used as folder names'''
    if len(key) > max_len:
        return key[0:max_len]
    else:
        return key


def enumerate_read_files(exp_dir, single_ended_reads, batch_name=None):
    if not exp_dir[-1] == "/":
        exp_dir = f"{exp_dir}/"
    accepted_formats = [".fq", ".fastq", ".gz"]
    if batch_name:
        exp_dir = f"{batch_name}/{exp_dir}"
    try:
        f_full = [f"{exp_dir}/{i}" for i in os.listdir(
            exp_dir) if any(subst in i for subst in accepted_formats)]
        f_full = [i for i in f_full if any(
            subst in f'.{i.split(".")[-1]}' for subst in accepted_formats)]
    except FileNotFoundError:
        stoperr(
            f"The directory you specified to look for files doesn't exist, check your ExpDir and re-run.")
    if len(f_full) == 2:
        return f_full
    elif len(f_full) == 1 and single_ended_reads:
        return f_full
    else:
        raise stoperr(
            f"I didn't find exactly 2 .fq/.fastq[.gz] read files in folder (ExpDir), so I'm skipping it: {exp_dir}.\n"
            f"If you're using single ended reads (e.g. Nanopore), set the 'SingleEndedReads' parameter to 'true', and try again.")


def enumerate_bam_files(exp_dir):
    accepted_formats = [".bam"]
    f_full = [f"{exp_dir}/{i}" for i in os.listdir(
        exp_dir) if any(subst in f".{i.split('.')[-1]}" for subst in accepted_formats)]
    assert len(
        f_full) == 1, f"ERROR: Please ensure there is a single .bam file in your experiment directory (ExpDir). I detected these .bam files: {f_full if not len(f_full) == 0 else 'None'}"
    return f_full[0]
