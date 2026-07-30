# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import types

# load the module as a spec
spec = importlib.util.spec_from_file_location("bias", "bin/publicCoreHelperFunctions.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

# test functions
@pytest.fixture
def constants():
    pytest.test_bam1 = "test/public/unit/omnition-core/files/test_publicRnaCalcUmiBias_bam.bam"
    pytest.test_dict = {
        "a": 1, "b": 2, "c": 3, "d": 4, "e": 5, "f": 6, "g": 7, "h": 8
    }
    pytest.chunks1 = 3
    pytest.chunks2 = 10

def test_get_chromosomes_with_alignments(constants):
    """
    Test that chromosomes with alignments are correctly identified.
    """
    chroms_test = mod.get_chromosomes_with_alignments(pytest.test_bam1)
    assert chroms_test == ['10']


def test_split_dict(constants):
    """
    Test that dictionaries are split evenly
    """
    split_dict1 = mod.split_dict(pytest.test_dict, pytest.chunks1)
    split_dict2 = mod.split_dict(pytest.test_dict, pytest.chunks2)
    assert isinstance(split_dict1, types.GeneratorType)
    assert isinstance(split_dict2, types.GeneratorType)
    split_dict_len1 = 0
    split_dict_len2 = 0
    for out_dict in split_dict1:
        assert len(out_dict.items()) <= 3
        split_dict_len1 += 1
    for out_dict in split_dict2:
        assert len(out_dict.items()) == 1
        split_dict_len2 += 1
    assert split_dict_len1 == pytest.chunks1
    assert split_dict_len2 != pytest.chunks2
    assert split_dict_len2 == 8

