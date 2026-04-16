import os
import shutil

from test.utils import get_random_str, make_rand_dir, get_default_args
from app.utils.utility_fns import read_fa
from app.src.consensus import Consensus
from app.src.analysis import Analysis
from app.src.generate_counts import run_counts
from app.src.map_reads_to_ref import run_map
from app.utils.mapping_ref_convert import MappingRefConverter


def init_map(p):
    run_map(p)


def init_generate_counts(p):
    run_counts(p)


def init_analysis(p, start_with_bam=False):
    clf = Analysis(p, start_with_bam)
    clf.main()


def init_consensus(infile=None):
    '''Init start files'''
    p = get_default_args()
    bamstart = False
    fstem = f"{p['SaveDir']}/{p['ExpName']}/"
    if infile:
        p["RefStem"] = infile
        p["ExpDir"] = "./data/eval/"
        p["MappingRefTable"] = infile.replace(".fasta", ".csv")
        bamstart = True
    if os.path.exists(fstem):
        shutil.rmtree(fstem)
    os.mkdir(fstem)
    p = MappingRefConverter(p, sneaky_mode=True).main()
    if not infile:
        init_map(p)
    run_counts(p)
    init_analysis(p)

    clf = Consensus(p, start_with_bam=bamstart)
    clf.main()
    if infile:
        '''Check aggregation on bacterial targets'''
        cons = read_fa(
            f"{p['folder_stem']}/consensus_sequences/streptococcus-pneumoniae_remapped_consensus_sequence.fasta")
        assert len(cons[0][1]) > 1500
    else:
        '''Else check a consensus was made of any length'''
        assert os.stat(
            f"{p['folder_stem']}/consensus_seq_stats.csv").st_size > 2

    shutil.rmtree(fstem)


if __name__ == "__main__":
    for i in ["./data/eval/agg_refs.fasta", None]:
        init_consensus(infile=i)
