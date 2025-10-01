import pickle as p
import sys
import re
import os

from app.utils.error_handlers import error_handler_parse_bam_positions, error_handler_cli
from app.utils.argparsers import parse_args_bam_parse
from app.utils.shell_cmds import make_dir, shell, loginfo
from app.utils.utility_fns import get_gene_orgid, trim_long_fpaths
from app.utils.basic_cli_calls import samtools_index


class Parse_bam_positions:
    '''
    Parse contents of bam file, via reading shell commands passed in.
    Count reads mapped to each target, generate read groupings for consensus calling.
    Should only get called by another Python script due to requirement for shell input.
    '''

    def __init__(self, argies) -> None:
        '''N.b. argies are a namespace (not a dict) because called from cli!!'''
        self.p = argies
        self.min_match_length = int(self.p['MatchLength'])
        self.n = 3  # Min n reads to decide we want to make a consensus
        self.minimum_n_filter = 1  # Filter unique reads if PostFilt
        self.reads_by_hit = {}
        self.fnames = {
            "bam": f"{self.p['SaveDir']}/{self.p['ExpName']}/{self.p['ExpName']}.bam",
            "bamview": f"{self.p['SaveDir']}/{self.p['ExpName']}/{self.p['ExpName']}_bamview.txt",
            "delreads": f"{self.p['SaveDir']}/{self.p['ExpName']}/{self.p['ExpName']}_reads_to_del.txt",
            "bamfilt": f"{self.p['SaveDir']}/{self.p['ExpName']}/{self.p['ExpName']}_filtered.bam",
            # TODO < Harmonise with fnames.py
            "grouped_reads": f"{self.p['SaveDir']}/{self.p['ExpName']}/grouped_reads.p",
        }

    def getmatchsize(self, cigar):
        '''Find matches in cigar string with regex, return count'''
        matches = re.findall(r'([0-9]+)M', cigar)
        if not len(matches):
            return 0
        else:
            return sum(int(x) for x in matches)

    def build_target_dbs(self, ref, seq, id):
        if not ref in self.reads_by_hit.keys():
            self.reads_by_hit[ref] = [[id, seq]]
        else:
            self.reads_by_hit[ref].append([id, seq])

    def parse_bam_position(self, fields):
        '''GENERATE COUNTS STAGE: For each line passed in, parse fields of interest; identify matches and print back to stdout.'''
        id, ref, pos, ref2, tlen, cigar, seq = fields[0], fields[2], fields[3], fields[6], int(
            fields[8]), fields[5], fields[9]

        ref_name_match = True if get_gene_orgid(
            ref2)[0] == "=" else get_gene_orgid(ref) == get_gene_orgid(ref2)
        match = tlen >= self.min_match_length
        # RM < TODO CHECK EVALUATES AND NOT PASSING A STRING
        if not self.p['SingleEndedReads']:
            improper_match = (tlen == 0) and (self.getmatchsize(
                cigar) >= self.min_match_length) and ref_name_match
        else:
            '''Experimental, for use with single ended sets (e.g. when Sequencer explodes mid-run)'''
            improper_match = (tlen == 0) and (self.getmatchsize(
                cigar) >= self.min_match_length)

        if improper_match and tlen == 0:
            try:
                tlen = self.getmatchsize(cigar)
            except:
                tlen = 0

        if match or improper_match:
            '''Properly paired and match is of decent mapped length OR
            Improperly paired BUT same gene AND match is of decent mapped length (via CIGAR string lookup) AND RNAME ref organism is same to RNEXT ref org'''
            self.build_target_dbs(ref, seq, id)
            return [ref, pos, tlen, id]
        else:
            return

    def save_hit_dbs(self):
        p.dump(self.reads_by_hit, open(
            self.fnames['grouped_reads'], "wb"), protocol=p.HIGHEST_PROTOCOL)
        grp_aln_f = f"{self.p['SaveDir']}/{self.p['ExpName']}/grouped_reads/"
        # TODO < Check exists, do with OS for safety
        make_dir(f"mkdir {grp_aln_f}")

        for key in self.reads_by_hit.keys():
            '''Save list of grouped read QNAME ids for calling consensuses, if more than n reads'''
            if len(self.reads_by_hit[key]) < self.n:  # TODO < Check still necessary
                '''Don't save grouped reads if less than n'''
                continue

            short_key = trim_long_fpaths(key)
            make_dir(f'mkdir "{grp_aln_f}{short_key}"')
            with open(f"{grp_aln_f}{short_key}/{short_key}.lst", "w") as file:
                [file.write(f"{self.reads_by_hit[key][i][0]}\n")
                    for i in range(len(self.reads_by_hit[key]))]

    def get_reads(self):
        '''Read serialised BAM file into memory, create unique indexes for vectorised matching with filter list'''
        headers = []
        dat = {}
        for l in open(self.fnames['bamview']):
            if l.startswith('@'):  # ignore headers
                continue

            res = self.parse_bam_position(l.split())
            if not res:
                continue
            if not f"{res[0]}_{res[1]}_{res[2]}" in dat.keys():
                dat[f"{res[0]}_{res[1]}_{res[2]}"] = []
            dat[f"{res[0]}_{res[1]}_{res[2]}"].append(
                [f'{res[0]},{res[1]},{res[2]},{self.p["ExpName"]}', res[3]])

        return dat, headers  # currently no use for headers

    def main(self):
        '''Entrypoint. Multi functional across generate counts and post filter.'''
        loginfo(f"Parsing BAM file {self.fnames['bam']}")
        dat, _ = self.get_reads()
        if self.p['PostFilt']:
            loginfo(
                f"Post filtering reads with less than {self.minimum_n_filter} reads")
            '''Filter data if < n reads (default = 1/no unique)'''
            '''Make new data dictionary'''
            dat = {k: v for k, v in dat.items() if len(v) >
                   self.minimum_n_filter}

            '''Kill original BAM and replace with one that's been filtered'''
            with open(self.fnames["delreads"], "w") as f:
                [f.write(f"{v[0][1]}\n") for v in dat.values()
                 if len(v) <= self.minimum_n_filter]
            out = shell(
                f"samtools view -@ {self.p['NThreads']} -b -N {self.fnames['delreads']} {self.fnames['bam']} > {self.fnames['bamfilt']}", is_test=True)
            error_handler_cli(out, self.fnames['bamfilt'], "samtools")
            os.remove(self.fnames['bam'])
            os.rename(self.fnames['bamfilt'], self.fnames['bam'])
            os.remove(self.fnames['delreads'])
            samtools_index(self.fnames['bam'])

        if len(self.reads_by_hit) == 0:
            '''No hits found between input BAM file and reference sequences'''
            return
        else:
            '''Save data for consensus call fns'''
            self.save_hit_dbs()

        loginfo(
            f"Parsed {len(self.reads_by_hit)} hits from BAM file {self.fnames['bam']}. Saving results...")
        p.dump(dat, open(self.p['counts_pickle'], "wb"))


if __name__ == '__main__':
    cls = Parse_bam_positions(parse_args_bam_parse())
    cls.main()
