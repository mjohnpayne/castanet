import os
import pytest
import shutil

from test.utils import make_rand_dir, get_default_args
from app.utils.mapping_ref_convert import MappingRefConverter


def init_payload(sneaky_mode):
    p = get_default_args()
    rand_dir = make_rand_dir()
    p["SaveDir"] = rand_dir
    p["ExpName"] = "testexp"
    fasta_file = f"{rand_dir}/test.fa"
    p["RefStem"] = fasta_file
    if not sneaky_mode:
        p["InFile"] = p["RefStem"]
        p["OutFile"] = f"{rand_dir}/ref.csv"

    return p, rand_dir, fasta_file


def init_mapping_ref_convert_not_a_fasta(sneaky_mode):
    '''Test MappingRefConverter with invalid fasta file'''
    p, rand_dir, fasta_file = init_payload(sneaky_mode)
    fasta_file = f"{rand_dir}/invalid.file"
    p["RefStem"] = fasta_file
    if not sneaky_mode:
        p["InFile"] = p["RefStem"]

    with open(fasta_file, "w") as f:
        f.write("""This is not a fasta file""")
    mrc = MappingRefConverter(p, sneaky_mode=sneaky_mode)
    with pytest.raises(SystemError):
        mrc.make_csv()
    shutil.rmtree(rand_dir)


def init_mapping_ref_convert_invalid_fasta(sneaky_mode):
    '''Test MappingRefConverter with invalid fasta file'''
    p, rand_dir, fasta_file = init_payload(sneaky_mode)
    with open(fasta_file, "w") as f:
        f.write("""This is not a fasta file""")
    mrc = MappingRefConverter(p, sneaky_mode=sneaky_mode)
    with pytest.raises(SystemError):
        mrc.make_csv()
    shutil.rmtree(rand_dir)


def init_mapping_ref_convert_disallowed_chars(sneaky_mode):
    '''Test MappingRefConverter input checks with disallowed characters'''
    p, rand_dir, fasta_file = init_payload(sneaky_mode)

    with open(fasta_file, "w") as f:
        f.write(""">probetype/with:disallowed*chars\nACGTACGTACGT""")
    p["RefStem"] = fasta_file
    mrc = MappingRefConverter(p, sneaky_mode=sneaky_mode)
    df = mrc.make_csv()
    with pytest.raises(SystemError):
        mrc.input_checks(df)
    shutil.rmtree(rand_dir)


def init_mapping_ref_convert_valid_fasta(sneaky_mode):
    '''Test MappingRefConverter with valid fasta file'''
    p, rand_dir, fasta_file = init_payload(sneaky_mode)

    with open(fasta_file, "w") as f:
        f.write(""">probetype_valid_description\nACGTACGTACGT""")
    mrc = MappingRefConverter(p, sneaky_mode=sneaky_mode)
    df = mrc.make_csv()
    df_checked = mrc.input_checks(df)
    assert not df_checked.empty
    shutil.rmtree(rand_dir)


def test_sneaky_mode():
    for i in [True, False]:
        init_mapping_ref_convert_not_a_fasta(i)
        init_mapping_ref_convert_invalid_fasta(i)
        init_mapping_ref_convert_disallowed_chars(i)
        init_mapping_ref_convert_valid_fasta(i)


if __name__ == "__main__":
    test_sneaky_mode()
