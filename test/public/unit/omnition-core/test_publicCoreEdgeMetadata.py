# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("correcting", "bin/publicCoreEdgeMetadata.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

@pytest.fixture
def constants():
    pytest.allowlist = "test/public/unit/omnition-core/files/test_publicCoreEdgeMetadata_allowlist.csv"
    pytest.corrected_edges = "test/public/unit/omnition-core/files/test_publicCoreEdgeMetadata_corrected_edges.tsv"
    pytest.edges = "test/public/unit/omnition-core/files/test_publicCoreEdgeMetadata_edges.tsv"
    pytest.sample_id = "test_sample"
    pytest.edge_type1 = "raw"
    pytest.edge_type2 = "corrected"
    pytest.edge_type3 = "above_knee"


def test_read_allowlist(constants):
    """
    Test that allowlist is read in correctly and column is labeled correctly
    """
    allowlist = mod.read_allowlist(pytest.allowlist)
    assert isinstance(allowlist, pl.DataFrame)
    assert allowlist.shape == (10, 2)
    assert allowlist.columns == ["rank", "bead"]


def test_read_edgelist(constants):
    """
    Test that edges are read in and reordered lexographically
    """
    edges = mod.read_edgelist(pytest.edges).collect()
    assert isinstance(edges, pl.DataFrame)
    assert edges.shape == (37, 6)
    assert edges.columns == ["edge_id", "edge_umi", "bead1", "umi1", "bead2", "umi2"]
    edges = edges.with_columns(
        pl.when(pl.col("bead1") <= pl.col("bead2"))
        .then(1)
        .otherwise(0)
        .alias("bead_test")
    )
    assert edges["bead_test"].sum() == 37
    assert edges.filter(pl.col("bead1") == pl.col("bead2")).shape[0] == 19
    edges = edges.filter(pl.col("bead1") == pl.col("bead2")).with_columns(
        pl.when(pl.col("umi1") <= pl.col("umi2"))
        .then(1)
        .otherwise(0)
        .alias("umi_test")
    )
    assert edges["umi_test"].sum() == 19


def test_filter_edgelist(constants):
    """
    Test that below knee edges are filtered out
    """
    edgelist = mod.read_edgelist(pytest.corrected_edges)
    allowlist = mod.read_allowlist(pytest.allowlist)
    filtered_edgelist = mod.filter_edgelist(
        edgelist, allowlist.lazy()
    )
    # check we return the right type
    assert isinstance(filtered_edgelist, pl.DataFrame)
    # check that the allowlist has different beads than the edgelist
    assert set(edgelist.collect()["bead1"]) != set(allowlist["bead"])
    assert set(edgelist.collect()["bead2"]) != set(allowlist["bead"])
    # check that after filtering, all bead barcodes in edgelist are on allowlist
    assert all(item in set(allowlist["bead"]) for item in filtered_edgelist["bead1"])
    assert all(item in set(allowlist["bead"]) for item in filtered_edgelist["bead2"])


def test_count_edge_types(constants):
    """
    Test that the edge types are counted correctly
    """
    counts = dict()
    edgelist = mod.read_edgelist(pytest.edges).collect()
    counts = mod.count_edge_types(edgelist, pytest.edge_type1, counts)
    corrected_edgelist = mod.read_edgelist(pytest.corrected_edges).collect()
    counts = mod.count_edge_types(corrected_edgelist, pytest.edge_type2, counts)
    allowlist = mod.read_allowlist(pytest.allowlist)
    filtered_edgelist = mod.filter_edgelist(
        corrected_edgelist.lazy(), allowlist.lazy()
    )
    counts = mod.count_edge_types(filtered_edgelist, pytest.edge_type3, counts)
    assert counts["total_edges"] == 37
    assert counts[f"{pytest.edge_type1}_unique_edges"] == 32
    assert counts[f"{pytest.edge_type1}_ambiguous_umis"] == 2
    assert counts[f"{pytest.edge_type1}_redundant_edges"] == 18
    assert counts[f"{pytest.edge_type2}_unique_edges"] == 27
    assert counts[f"{pytest.edge_type2}_ambiguous_umis"] == 1
    assert counts[f"{pytest.edge_type2}_redundant_edges"] == 22
    assert counts[f"{pytest.edge_type3}_edges"] == 12
    assert counts[f"{pytest.edge_type3}_unique_edges"] == 9
    assert counts[f"{pytest.edge_type3}_ambiguous_umis"] == 1
    assert counts[f"{pytest.edge_type3}_redundant_edges"] == 7


def test_create_final_summary(constants):
    """
    Test that we create a final summary dataframe
    Test that summary percents are calculated correctly
    """
    counts = dict()
    edgelist = mod.read_edgelist(pytest.edges).collect()
    counts = mod.count_edge_types(edgelist, pytest.edge_type1, counts)
    corrected_edgelist = mod.read_edgelist(pytest.corrected_edges).collect()
    counts = mod.count_edge_types(corrected_edgelist, pytest.edge_type2, counts)
    allowlist = mod.read_allowlist(pytest.allowlist)
    filtered_edgelist = mod.filter_edgelist(
        corrected_edgelist.lazy(), allowlist.lazy()
    )
    counts = mod.count_edge_types(filtered_edgelist, pytest.edge_type3, counts)
    summary = mod.create_final_summary(counts, pytest.sample_id)
    assert isinstance(summary, pl.DataFrame)
    assert summary.shape == (14, 4)
    assert summary.columns == [
        "sample",
        "process",
        "metric",
        "value"
    ]
    assert summary.filter(pl.col("metric") == "percent_raw_unique_deconvolution_reads")["value"][0] == 86.49
    assert summary.filter(pl.col("metric") == "percent_corrected_unique_deconvolution_reads")["value"][0] == 72.97
    assert summary.filter(pl.col("metric") == "percent_above_knee_unique_deconvolution_reads")["value"][0] == 24.32
