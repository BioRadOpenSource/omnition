# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("correcting", "bin/rnaFilterBeads.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)



# test filtering for edge sequence
def test_filter_edgelist():
    edgelist = mod.read_edgelist(
        "test/unit/omnition-core/files/test_filter_edge_file.tsv"
    ).collect()
    allowlist = mod.read_allowlist(
        "test/unit/omnition-core/files/test_barcode_file.csv", True
    )
    filtered_edgelist = mod.filter_edgelist(
        edgelist.lazy(), allowlist.lazy()
    )
    # check we return the right type
    assert isinstance(filtered_edgelist, pl.DataFrame)
    # check that the allowlist has different beads than the edgelist
    assert set(edgelist["bead1"]) != set(allowlist["droplet"])
    assert set(edgelist["bead2"]) != set(allowlist["droplet"])
    # check that after filtering, all bead barcodes in edgelist are on allowlist
    assert all(item in set(allowlist["droplet"]) for item in filtered_edgelist["bead1"])
    assert all(item in set(allowlist["droplet"]) for item in filtered_edgelist["bead2"])
