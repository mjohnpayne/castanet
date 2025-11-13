import os
import pickle
import shutil
from app.utils.utility_fns import enumerate_read_files
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import stoperr, logerr
from app.utils.error_handlers import error_handler_cli


def run_map(p, is_test=False):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    '''Default in_files are created by the trimming step'''
    in_files = [f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq",
                f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"]
    CLEAN_UP = True

    # TODO < So much redundancy here. Can be refactored into single loop

    '''Check input files exist and are non-empty'''
    for fn in in_files:
        '''If default infiles not present, look in ExpDir for user-specified ones'''
        if not os.path.exists(fn):
            in_files = enumerate_read_files(
                f"{p['ExpDir']}/", p["SingleEndedReads"])
            CLEAN_UP = False
    for fn in in_files:
        if os.stat(fn).st_size < 2:
            stoperr(
                f"Castanet found an input file: {fn}, but it's empty. Please check your input file have been processed appropriately for input to BWA-mem2.")

    '''Calculate total n trimmed reads, for later calculation of read proportions (in analysis.py)'''
    try:
        out = int(
            shell(f"echo $(cat {in_files[0]}|wc -l)/4|bc", ret_output=True).decode("utf-8"))
        pickle.dump(
            out, open(f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_rawreadnum.p", "wb"))
    except Exception as e:
        logerr(f"Failed to calculate n trimmed reads from fastq files. Defaulting to calculating from total reads in bam.")

    if p["Mapper"] == "bwa":
        end_sec_print(
            f"INFO: Beginning initial mapping using BWA\nThis may take a while for large files")
        out = shell(f"bwa-mem2 index {p['RefStem']}", is_test=True)
        if p["SingleEndedReads"]:
            shell(
                f"bwa-mem2 mem -t {p['NThreads']} {p['RefStem']} {in_files[0]} | samtools view -F4 -Sb - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")
            shell(  # This is done out of sync with CLEAN_UP as we can't assume user has not transferred single file to exp directory
                f"rm {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq")

        else:
            shell(
                f"bwa-mem2 mem -t {p['NThreads']} {p['RefStem']} {in_files[0]} {in_files[1]} | samtools view -F4 -Sb - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")

        error_handler_cli(
            out, f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam", "bwa-mem2", test_f_size=True)

        if CLEAN_UP:
            '''Remove mapping ref indices'''
            shell(f"rm {p['SaveDir']}/{p['ExpName']}/ref.fa.*")

    elif p["Mapper"] == "bowtie2":
        end_sec_print(
            f"INFO: Beginning initial mapping using Bowtie2\nThis may take a while for large files")
        ref_dir = f"{p['SaveDir']}/reference_indices"
        if not os.path.exists(ref_dir):
            os.mkdir(ref_dir)
        out = shell(
            f"bowtie2-build --large-index {p['RefStem']} {p['SaveDir']}/reference_indices", is_test=True)
        if p["SingleEndedReads"]:
            shell(f"bowtie2 -x {p['SaveDir']}/reference_indices -U {in_files[0]} -p {p['NThreads']} --local -I 50 --maxins 2000 --no-unal | samtools view -@ {p['NThreads']} -h -Sb -F4 -F2048 - | samtools sort -@ {p['NThreads']} - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")
        else:
            shell(f"bowtie2 -x {p['SaveDir']}/reference_indices -1 {in_files[0]} -2 {in_files[1]} -p {p['NThreads']} --local -I 50 --maxins 2000 --no-unal | samtools view -@ {p['NThreads']} -h -Sb -F4 -F2048 - | samtools sort -@ {p['NThreads']} - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")
        error_handler_cli(
            out, f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam", "bowtie2", test_f_size=True)

    elif p["Mapper"] == "minimap2":
        end_sec_print(
            f"INFO: Beginning initial mapping using Minimap2\nThis may take a while for large files")

        if p["SingleEndedReads"]:
            out = shell(
                f"minimap2 -ax map-ont {p['RefStem']} {in_files[0]} | samtools view -F4 -Sb - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam",  is_test=True)
            shell(  # This is done out of sync with CLEAN_UP as we can't assume user has not transferred single file to exp directory
                f"rm {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq")

        else:
            out = shell(
                f"minimap2 -a {p['RefStem']} {in_files[0]} {in_files[1]} | samtools view -F4 -Sb - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam",  is_test=True)

        error_handler_cli(
            out, f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam", "minimap2", test_f_size=True)

    else:
        stoperr(
            f"User option for mapping software ('Mapper' parameter) is not 'bwa' or 'bowtie2' (you specified '{p['Mapper']}'), so I can't proceed")

    if CLEAN_UP:
        shell(
            f"rm {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq")

    end_sec_print(f"INFO: Mapping complete")
