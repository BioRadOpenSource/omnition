# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("merging", "bin/rnaMergeBeads.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

@pytest.fixture
def constants():
    pytest.edge_file = "test/unit/omnition-core/files/DemoRnaMerge_S1_edges.tsv"
    pytest.large_edge_file = "test/unit/omnition-core/files/F39-Rep1_S4_edges.tsv"
    pytest.small_edge_file = "test/unit/omnition-core/files/test_edges.tsv"
    pytest.bead_file = "test/unit/omnition-core/files/test_bead_file.tsv"

def test_read_edgelist(constants):
    bead_dict = mod.read_edgelist(pytest.large_edge_file)
    assert isinstance(bead_dict, dict)
    # number of beads that detected another bead with a non-ambiguous umi
    assert len(bead_dict) == 1748

def test_build_droplets_1umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 1) # no edges are removed
    assert(set(bead_dict.keys()) == set(['1', '2', '3', '4', '5', '6', '7', '8', '10', '9']))
    bead_dict = mod.format_edges(bead_dict)
    droplets = mod.build_droplets(bead_dict, False)
    assert(droplets[1] == mod.Droplet(beads={'1','2'}))
    assert(droplets[2] == mod.Droplet(beads={'3','4'}))
    assert(droplets[3] == mod.Droplet(beads={'5','6'}))
    assert(droplets[4] == mod.Droplet(beads={'7','8','9','10'}))

def test_build_droplets_2umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 2)
    assert(set(bead_dict.keys()) == set(['1', '2', '3', '4', '5', '6', '7', '8', '10', '9']))
    # bead 5 and 6 lose all edges, edge 7 and 10 are disconnected
    assert(bead_dict['5'] == {})
    assert(bead_dict['6'] == {})
    assert(set(bead_dict['7'].keys()) == set('8'))
    bead_dict = mod.format_edges(bead_dict)
    droplets = mod.build_droplets(bead_dict, False)
    assert(droplets[1] == mod.Droplet(beads={'1','2'}))
    assert(droplets[2] == mod.Droplet(beads={'3','4'}))
    assert(droplets[3] == mod.Droplet(beads={'5'}))
    assert(droplets[4] == mod.Droplet(beads={'6'}))
    assert(droplets[5] == mod.Droplet(beads={'7','8','9','10'}))

def test_build_droplets_3umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 3)
    assert(set(bead_dict.keys()) == set(['1', '2', '3', '4', '5', '6', '7', '8', '10', '9']))
    # all beads disconnected except 1 and 2, 7 and 8
    assert(bead_dict['1'].keys() == set('2'))
    assert(bead_dict['7'].keys() == set('8'))
    bead_dict = mod.format_edges(bead_dict)
    droplets = mod.build_droplets(bead_dict, False)
    assert(droplets[1] == mod.Droplet(beads={'1','2'}))
    assert(droplets[2] == mod.Droplet(beads={'3'}))
    assert(droplets[3] == mod.Droplet(beads={'4'}))
    assert(droplets[4] == mod.Droplet(beads={'5'}))
    assert(droplets[5] == mod.Droplet(beads={'6'}))
    assert(droplets[6] == mod.Droplet(beads={'7','8'}))
    assert(droplets[7] == mod.Droplet(beads={'10'}))
    assert(droplets[8] == mod.Droplet(beads={'9'}))

def test_build_droplets_4umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 4)
    assert(set(bead_dict.keys()) == set(['1', '2', '3', '4', '5', '6', '7', '8', '10', '9']))
    # all beads disconnected
    for key in bead_dict.keys():
        bead_dict.get(key) == {}
    bead_dict = mod.format_edges(bead_dict)
    droplets = mod.build_droplets(bead_dict, False)
    assert(droplets[1] == mod.Droplet(beads={'1'}))
    assert(droplets[2] == mod.Droplet(beads={'2'}))
    assert(droplets[3] == mod.Droplet(beads={'3'}))
    assert(droplets[4] == mod.Droplet(beads={'4'}))
    assert(droplets[5] == mod.Droplet(beads={'5'}))
    assert(droplets[6] == mod.Droplet(beads={'6'}))
    assert(droplets[7] == mod.Droplet(beads={'7'}))
    assert(droplets[8] == mod.Droplet(beads={'8'}))
    assert(droplets[9] == mod.Droplet(beads={'10'}))
    assert(droplets[10] == mod.Droplet(beads={'9'}))

def test_order_barcode_umi():
    b1 = 'A'
    u1 = 'UMI1'
    b2 = 'B'
    u2 = 'UMI2'
    bead1, bead2, umi = mod.order_barcode_umi(b1, u1, b2, u2)
    assert(bead1 == b1)
    assert(bead2 == b2)
    assert(umi == f"{u1}{u2}")

def test_build_one_to_one_droplets_1umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 1)
    droplets = mod.build_one_to_one_droplets(bead_dict)
    assert(len(droplets) == len(bead_dict))

def test_build_one_to_one_droplets_2umi(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 2)
    droplets = mod.build_one_to_one_droplets(bead_dict)
    assert(len(droplets) == len(bead_dict))


def test_compare_bead_files(constants):
    bead_dict = mod.read_edgelist(pytest.small_edge_file)
    bead_dict = mod.prune_edges(bead_dict, 1)
    droplets = mod.build_one_to_one_droplets(bead_dict)
    bead_list = mod.read_bead_file(pytest.bead_file)
    new_bead_list = mod.compare_bead_files(bead_list, droplets)
    assert(len(new_bead_list) == 2)
