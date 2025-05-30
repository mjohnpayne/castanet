import os
import pickle
from app.utils.utility_fns import enumerate_read_files
from app.utils.shell_cmds import shell
from app.utils.system_messages import end_sec_print
from app.utils.shell_cmds import stoperr, logerr
from app.utils.error_handlers import error_handler_cli


def run_map(p):
    '''Use BWA and Samtools to map reads from each sample to targets'''
    '''Default in_files are created by the trimming step'''
    in_files = [f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_1_clean.fastq",
                f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_2_clean.fastq"]
    CLEAN_UP = True
    for fn in in_files:
        '''If default infiles not present, look in ExpDir for user-specified ones'''
        if not os.path.exists(fn):
            in_files = enumerate_read_files(f"{p['ExpDir']}/")
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
        shell(f"bwa-mem2 mem -t {p['NThreads']} {p['RefStem']} {in_files[0]} {in_files[1]} | samtools view -F4 -Sb - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")
        error_handler_cli(
            out, f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam", "bwa-mem2", test_f_size=True)

    elif p["Mapper"] == "bowtie2":
        end_sec_print(
            f"INFO: Beginning initial mapping using Bowtie2\nThis may take a while for large files")
        ref_dir = f"{p['SaveDir']}/reference_indices"
        if not os.path.exists(ref_dir):
            os.mkdir(ref_dir)
        out = shell(
            f"bowtie2-build --large-index {p['RefStem']} {p['SaveDir']}/reference_indices", is_test=True)
        shell(f"bowtie2 -x {p['SaveDir']}/reference_indices -1 {in_files[0]} -2 {in_files[1]} -p {p['NThreads']} --local -I 50 --maxins 2000 --no-unal | samtools view -h -Sb -F4 -F2048 - | samtools sort - 1> {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam")
        error_handler_cli(
            out, f"{p['SaveDir']}/{p['ExpName']}/{p['ExpName']}.bam", "bowtie2", test_f_size=True)

    else:
        stoperr(
            f"User option for mapping software ('Mapper' parameter) is not 'bwa' or 'bowtie2' (you specified '{p['Mapper']}'), so I can't proceed")

    if CLEAN_UP:
        shell(
            f"rm {p['SaveDir']}/{p['ExpName']}/{p['ExpName']}_[12]_clean.fastq")
    end_sec_print(f"INFO: Mapping complete")
