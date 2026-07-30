# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import os
import ray

# load the module as a spec
spec = importlib.util.spec_from_file_location("tag", "bin/publicAtacBamTagger.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# create parallel pool for ray
ray.init(num_cpus=1, ignore_reinit_error = True)

# execute test functions
def test_element_readname():

    name = "AAAAAAA_AV232901:CZ_TagCentricATAC_59CZ_07162024:2320448247:1:10102:0097:0058 1:N:0:TAGGCATG"
    
    (barcode, readname) = mod.extract_barcode(name, 0)

    assert(barcode == "AAAAAAA")
    assert(readname == "AV232901:CZ_TagCentricATAC_59CZ_07162024:2320448247:1:10102:0097:0058 1:N:0:TAGGCATG")


def test_illumina_readname():

    name = "AAAAAAA_NB551371:158:HFWKGAFX2:1:21209:7118:19361"
    
    (barcode, readname) = mod.extract_barcode(name, 0)

    assert(barcode == "AAAAAAA")
    assert(readname == "NB551371:158:HFWKGAFX2:1:21209:7118:19361")

