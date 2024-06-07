# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest

# load the module as a spec
spec = importlib.util.spec_from_file_location("bias", "bin/atacAssembleFragments.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

# test functions


def test_get_read_count_list():
    test_bam = "test/unit/omnition-core/files/test.bam"
    result = mod.get_read_count_list(sample = 'XC', chromosome = 'XC', input_count = 30, output_count = 20)
    assert len(result) == 3
    assert type(result) is list
    assert result[1] == ['XC', "assemble_fragments_XC", "input", 30]
    assert result[2] == ['XC', "assemble_fragments_XC", "output", 20]
