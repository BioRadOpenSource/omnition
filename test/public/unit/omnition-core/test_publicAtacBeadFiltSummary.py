# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("split", "bin/publicAtacBeadFiltSummary.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.fixture
def constants():
    pytest.input_dir = "test/public/unit/omnition-core/files/"
    pytest.sample_id = "SampleA"
    pytest.qc_stats_suffix = ".QCstats.csv"
    pytest.fastq_ti_counts_file = "test_publicAtacBeadFiltSummary.fastqTIreadcounts.csv"
    pytest.expected_names = ['test_publicAtacBeadFiltSummary', 'test_publicAtacBeadFiltSummary2']
    pytest.ti_read_counts = {
        'DemoAtacCombinatorial_S1-CGCATA': 38343,
        'DemoAtacCombinatorial_S1-AGCACG': 37227
    }
    pytest.index_metrics_dict = {
        'sample': ['test_publicAtacBeadFiltSummary', 'test_publicAtacBeadFiltSummary2'],
        'process': ['deconvolution', 'deconvolution'],
        'metric': ['mean_valid_read_pairs_per_cell', 'mean_valid_read_pairs_per_cell'],
        'value': ['116.3', '101.0']
    }
    pytest.full_metrics_dict = {
        'sample': ['test_publicAtacBeadFiltSummary', 'test_publicAtacBeadFiltSummary2', 'SampleA', 'SampleA', 'SampleA', 'SampleA', 'SampleA', 'SampleA'],
        'process': ['deconvolution', 'deconvolution', 'summary', 'summary', 'summary', 'summary', 'summary', 'summary'],
        'metric': ['mean_valid_read_pairs_per_cell', 'mean_valid_read_pairs_per_cell', 'total_reads', 'median_reads_per_index', 'total_cells', 'mean_valid_read_pairs_per_cell', 'median_unique_nuclear_fragments_per_cell', 'median_estimated_library_size_per_cell'],
        'value': ['116.3', '101.0', 75570, 37785.0, 185, 408.5, 210.0, 899.0]
    }


def test_find_fastq_ti_names(constants):
    """
    Test the find_fastq_ti_names function to ensure it correctly identifies fastq+TI names.
    """
    names = mod.find_fastq_ti_names(pytest.input_dir, pytest.qc_stats_suffix)
    assert names == pytest.expected_names, f"Expected {pytest.expected_names}, but got {names}"


def test_read_ti_counts(constants):
    """
    Test the read_ti_counts function to ensure it correctly reads the fastqTIreadcounts file.
    """
    ti_read_counts = mod.read_ti_counts(pytest.sample_id, pytest.input_dir, pytest.fastq_ti_counts_file)
    
    # Check if the returned dictionary has the expected structure
    assert isinstance(ti_read_counts, dict), "Expected a dictionary from read_ti_counts"
    for name, count in ti_read_counts.items():
        assert isinstance(name, str), "Expected keys to be strings"
        assert isinstance(count, int), "Expected values to be integers"
    print(ti_read_counts)
    assert ti_read_counts == pytest.ti_read_counts, f"Expected counts to match predefined values, but got {ti_read_counts}"


def test_collect_fastq_ti_metrics(constants):
    """
    Test the collect_fastq_ti_metrics function to ensure it collects metrics correctly.
    """
    metrics_dict = {
        "sample": [],
        "process": [],
        "metric": [],
        "value": [],
    }
    for name in pytest.expected_names:
        metrics_dict, qc_stats = mod.collect_fastq_ti_metrics(metrics_dict, name, pytest.input_dir, 10000, pytest.qc_stats_suffix)
        assert isinstance(qc_stats, pl.DataFrame), "Expected qc_stats to be a Polars DataFrame"
    assert metrics_dict == pytest.index_metrics_dict, f"Expected metrics_dict to match predefined values, but got {metrics_dict}"


def test_collect_sample_metrics(constants):
    """
    Test the collect_sample_metrics function to ensure it collects sample-level metrics correctly.
    """
    ti_read_counts = mod.read_ti_counts(pytest.sample_id, pytest.input_dir, pytest.fastq_ti_counts_file)
    metrics_dict = {
        "sample": [],
        "process": [],
        "metric": [],
        "value": [],
    }
    qc_stats = []
    for name in pytest.expected_names:
        metrics_dict, qc_stats_df = mod.collect_fastq_ti_metrics(metrics_dict, name, pytest.input_dir, 10000, pytest.qc_stats_suffix)
        qc_stats.append(qc_stats_df)
    
    


    metrics_dict = mod.collect_sample_metrics(
        pytest.sample_id, metrics_dict, ti_read_counts, qc_stats
    )
    assert isinstance(metrics_dict, dict), "Expected metrics_dict to be a dictionary"
    assert metrics_dict == pytest.full_metrics_dict, f"Expected metrics_dict to match predefined values, but got {metrics_dict}"
