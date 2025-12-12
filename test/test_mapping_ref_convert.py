import pytest
import shutil

from test.utils import make_rand_dir, create_test_file, get_default_args
from app.utils.mapping_ref_convert import MappingRefConverter


def test_mapping_ref_convert_invalid_fasta():
    '''Test MappingRefConverter with invalid fasta file'''
    p = get_default_args()
    rand_dir = make_rand_dir()
    invalid_fasta = f"{rand_dir}/invalid.fa"
    create_test_file(invalid_fasta, "This is not a fasta file.")
    p["RefStem"] = invalid_fasta
    p["SaveDir"] = rand_dir
    p["ExpName"] = "testexp"
    mrc = MappingRefConverter(p, sneaky_mode=True)
    with pytest.raises(ValueError):
        mrc.make_csv()
    shutil.rmtree(rand_dir)


def test_mapping_ref_convert_disallowed_chars():
    '''Test MappingRefConverter input checks with disallowed characters'''
    p = get_default_args()
    rand_dir = make_rand_dir()
    fasta_file = f"{rand_dir}/test.fa"
    fasta_content = """>probetype/with:disallowed*chars\nACGTACGTACGT"""
    create_test_file(fasta_file, fasta_content)
    p["RefStem"] = fasta_file
    p["SaveDir"] = rand_dir
    p["ExpName"] = "testexp"
    mrc = MappingRefConverter(p, sneaky_mode=True)
    df = mrc.make_csv()
    with pytest.raises(ValueError):
        mrc.input_checks(df)
    shutil.rmtree(rand_dir)


def test_mapping_ref_convert_valid_fasta():
    '''Test MappingRefConverter with valid fasta file'''
    p = get_default_args()
    rand_dir = make_rand_dir()
    fasta_file = f"{rand_dir}/valid.fa"
    fasta_content = """>probetype_valid_description\nACGTACGTACGT"""
    create_test_file(fasta_file, fasta_content)
    p["RefStem"] = fasta_file
    p["SaveDir"] = rand_dir
    p["ExpName"] = "testexp"
    mrc = MappingRefConverter(p, sneaky_mode=True)
    df = mrc.make_csv()
    df_checked = mrc.input_checks(df)
    assert not df_checked.empty


if __name__ == "__main__":
    test_mapping_ref_convert_invalid_fasta()
    test_mapping_ref_convert_disallowed_chars()
    test_mapping_ref_convert_valid_fasta()
