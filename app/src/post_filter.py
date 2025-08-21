from app.utils.shell_cmds import shell, loginfo
from app.utils.system_messages import end_sec_print
from app.src.analysis import Analysis
from app.src.consensus import Consensus

import os
import pandas as pd


def run_analysis(p):
    cls = Analysis(p, start_with_bam=False, is_post_filt=True)
    cls.main()


def run_consensus(p):
    clf = Consensus(p, start_with_bam=False)
    clf.main()


def run_post_filter(p):
    end_sec_print(f"INFO: Starting BAM post filter")
    p["ExpDir"] = f"{p['ExpDir']}/"
    bamview_fname = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_bamview.txt"
    filter_fname = f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_reads_to_drop.csv"

    assert os.path.exists(filter_fname)
    if not pd.read_csv(filter_fname).empty:

        if not os.path.exists(bamview_fname):
            loginfo(f"Regenerating bamview for post-analysis filter")
            shell(
                f"""samtools view -F2048 -F4 {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam > {bamview_fname}""")

        loginfo(f"Running post-analysis filter. This can be a slow step.")
        shell(
            f"python3 -m app.src.parse_bam -SingleEnded False -SaveDir {p['SaveDir']} -SeqName {p['ExpName']} -MatchLength {p['MatchLength']} -Mode filter -ExpDir {p['ExpDir']} -FilterFile {filter_fname} -ExpName {p['ExpName']}")

        assert os.path.exists(bamview_fname.replace("bamview.txt", "bamview_tmp.txt")
                              ), f"Post-filter didn't produce a BAM file. Does your reads to drop file contain anything?"

        os.remove(bamview_fname)
        os.rename(
            f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_bamview_tmp.txt", bamview_fname)
        loginfo(f"Post-filter complete. Re-running analysis and consensus generation.")
        shell(
            f"python3 -m app.src.parse_bam -SingleEnded False -SaveDir {p['SaveDir']} -SeqName {p['ExpName']} -MatchLength {p['MatchLength']} -Mode parse -ExpDir {p['ExpDir']} -ExpName {p['ExpName']} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g  > {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_PosCounts.csv")

        run_analysis(p)

    else:
        os.remove(filter_fname)
        with open(f"error_{p['ExpName']}", "w") as f:
            f.write("This file had no reads after applying your filter")

    run_consensus(p)
    # TODO Parameterise matchlen
    if os.path.exists(bamview_fname):
        os.remove(bamview_fname)
    end_sec_print(f"INFO: Post filter complete")
