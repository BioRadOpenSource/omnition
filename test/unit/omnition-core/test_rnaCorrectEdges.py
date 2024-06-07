# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
from umi_tools import UMIClusterer
import math
import ray

# load the module as a spec
spec = importlib.util.spec_from_file_location("correcting", "bin/rnaCorrectEdges.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

# create parallel pool for ray
ray.init(num_cpus=1, ignore_reinit_error=True)

@pytest.fixture
def constants():
    pytest.edge_file = "test/unit/omnition-core/files/F39-Rep1_S4_edges.tsv"

def test_make_edge_dict(constants):
    bead_dict = mod.make_edge_dict(pytest.edge_file)
    assert isinstance(bead_dict, dict)
    assert len(bead_dict) == 1648

def test_correct_umis(constants):
    bead_dict = mod.make_edge_dict(pytest.edge_file)
    c_umi_edges,corrected_umi_count = ray.get(mod.correct_umis.remote(bead_dict, 1))
    assert len(c_umi_edges) == len(bead_dict)
    assert corrected_umi_count == 2
