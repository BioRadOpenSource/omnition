# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest

# load the module as a spec
spec = importlib.util.spec_from_file_location("bias", "bin/rnaCalcUmiBias.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

# test functions


def test_bam_to_dict():
    test_bam = "test/unit/omnition-core/files/test.bam"
    result = mod.inputbam(test_bam, barcode_tag = 'XC', umi_tag = 'XM')
    assert len(result) == 3
    assert type(result) is tuple
