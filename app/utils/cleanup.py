from app.utils.shell_cmds import shell


def files_to_kill():
    return [
        "./SAVEDIR/EXPNAME/EXPNAME_PosCounts.csv",
        "./SAVEDIR/EXPNAME/EXPNAME_rawreadnum.p",
        "./SAVEDIR/EXPNAME/EXPNAME.bam",
        "./SAVEDIR/EXPNAME/EXPNAME.bai",
        "./SAVEDIR/EXPNAME/probe_aggregation.csv",
        "./SAVEDIR/EXPNAME/probe_lengths.csv",

    ]


def clean_intermediates(payload):
    """
    Clean up intermediate files created during the analysis.
    This function is called when DebugMode is enabled.
    """
    shell(f"rm -rf {payload['SaveDir']}/{payload['ExpName']}/grouped_reads")
    for file in files_to_kill():
        shell(
            f"rm {file.replace('SAVEDIR', payload['SaveDir']).replace('EXPNAME', payload['ExpName'])}")
