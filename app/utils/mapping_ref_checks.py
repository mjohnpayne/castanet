from app.utils.utility_fns import read_fa, stoperr, logerr
from app.utils.shell_cmds import loginfo


def check_mapping_ref(refstem):
    loginfo(f"Checking mapping reference file: {refstem}...")
    refs = read_fa(refstem)
    if len(refs) == 0:
        stoperr(
            f"Mapping reference file {refstem} appears to be empty or not in fasta format. Please check your mapping reference file.")

    seen = []
    seen_short = []
    for ref in refs:
        if ref[0] in seen:
            stoperr(
                f"Mapping reference {ref[0]} appears more than once in your mapping reference file. Please ensure all reference names are unique.")
        seen.append(ref[0])

        if len(ref[0]) > 100:
            if ref[0][0:100] in seen_short:
                stoperr(
                    f"Mapping reference {ref[0]} appears to have a non-unique first 100 characters in your mapping reference file. Please ensure all reference names are unique, and less than 100 characters long to avoid issues with downstream tools.")
            seen_short.append(ref[0][0:100])

        if not ref[0][0] == ">":
            stoperr(
                f"Mapping reference header {ref[0]} is malformed, please fix and re-run.")
        if len(ref[1]) < 100:
            logerr(
                f"Mapping reference {ref[0]} is very short ({len(ref[1])}bp) and may cause mapping issues.")
        if ' ' in ref[0] or '\t' in ref[0]:
            stoperr(
                f"Mapping reference {ref[0]} contains spaces or tabs. Please remove these from the header line of your fasta file. Suggest replacing with underscores, but note that aggregation will happen on first underscore.")
        if len(ref) > 100:
            stoperr(
                f"Mapping reference {ref[0]} has a very long header line (>100 characters). Please shorten this to avoid issues with downstream tools.")
    loginfo(
        f"Castanet has verified that your mapping reference file (RefStem) appears valid.")
