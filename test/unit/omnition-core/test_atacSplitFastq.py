# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import ray
import os

# load the module as a spec
spec = importlib.util.spec_from_file_location("split", "bin/atacSplitFastq.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# create parallel pool for ray
ray.init(num_cpus=1, ignore_reinit_error=True)

@pytest.fixture
def constants():
    pytest.readcounts_file = "test/unit/omnition-core/files/test_atacSplitFastq_readcounts.tsv"
    pytest.read1_file = "test/unit/omnition-core/files/test_atacSplitFastq_R1_001.complete_debarcoded.fastq.gz"
    pytest.sample_id = "test"

def test_parse_counts_file(constants):
    """
    Read test file and ensure that it's a list with two TIs
    """
    parsed_file = mod.parse_counts_file(pytest.readcounts_file, "test")
    assert isinstance(parsed_file, list)
    assert len(parsed_file) == 1

def test_parse_fastq_no_override(constants):
    """
    Test parsing a FASTQ file without error override.
    Expect ValueError because TI for second read in FASTQ
    is not present in the counts file.
    """
    parsed_file = mod.parse_counts_file(pytest.readcounts_file, "test")
    with pytest.raises(ValueError):
        result = ray.get(mod.parse_fastq.remote(pytest.sample_id, pytest.read1_file, parsed_file, False))


def test_parse_fastq_override(constants):
    """
    Test parsing a FASTQ file with error override.
    Expect no error because TI for second read in FASTQ.
    Expect one read to be written to output file.
    Expect no file to exist for read with ignored TI.
    """
    parsed_file = mod.parse_counts_file(pytest.readcounts_file, "test")
    result = ray.get(mod.parse_fastq.remote(pytest.sample_id, pytest.read1_file, parsed_file, True))
    lines = count_lines(f"{pytest.sample_id}-{parsed_file[0]}_R1.complete_debarcoded.split.fastq")
    assert lines == 4
    assert not os.path.exists(f"{pytest.sample_id}-TTTTT_R1.complete_debarcoded.split.fastq")


def count_lines(filename):
    """
    Count the number of lines in a file
    """
    with open(filename, "r") as f:
        for i, l in enumerate(f):
            pass
    return i + 1
