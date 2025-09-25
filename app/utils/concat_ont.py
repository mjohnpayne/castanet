import os
from app.utils.shell_cmds import shell, stoperr, loginfo


class ConcatOnt:
    def __init__(self, in_dir, out_dir, allowed_format):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.allowed_format = allowed_format
        if self.in_dir == self.out_dir:
            raise ValueError(
                "Input and output directories must be different. This prevents you from overwriting your input sequence data.")

    def enumerate_files(self, sub_dir) -> list:
        '''List all files in subdirectory with allowed_format (default .fastq.gz).'''
        return [os.path.join(sub_dir, file) for file in os.listdir(sub_dir) if file.endswith(self.allowed_format)]

    def make_check_dirs(self):
        '''Check if output directory exists, if not create it.'''
        if os.path.exists(self.out_dir):
            raise stoperr(
                f"Output directory {self.out_dir} already exists. Please remove it or choose another output directory. This prevents you from overwriting your sequence data.")
        os.mkdir(self.out_dir)

    def main(self):
        '''Iterate over subdirectories in input directory, concatenate all files with allowed_format in each subdirectory, write to output directory.'''
        self.make_check_dirs()
        dirs = [i for i in os.listdir(self.in_dir) if os.path.isdir(
            os.path.join(self.in_dir, i))]
        for dir in dirs:
            files = self.enumerate_files(os.path.join(self.in_dir, dir))
            if len(files) == 0:
                print(
                    f"No {self.allowed_format} files found in {dir}, skipping.")
                continue
            else:
                new_subdir = os.path.join(self.out_dir, dir)
                new_fname = f"{self.out_dir}/{dir}/{dir}_concatenated{self.allowed_format}"
                os.mkdir(new_subdir)
                shell(f"cat {' '.join(files)} > {new_fname}")
                print(f"Concatenated {len(files)} files into {new_fname}")
        loginfo(
            f"Finished concatenating files. Concatenated files are in {self.out_dir}.")
