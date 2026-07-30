# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("counts", "bin/publicAtacMergeReannCounts.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

@pytest.fixture
def constants():
    pytest.suffix = "_read_counts.csv"
    pytest.input_dir = "test/public/unit/omnition-core/files/"
    pytest.sample_id = "DemoAtacNormal_S1"
    pytest.process = "assemble_fragments"
    pytest.metric = "input"


def test_collect_counts(constants):
    """
    Test the collect_counts function to ensure it reads and concatenates CSV files correctly.
    """
    result = mod.collect_counts(pytest.input_dir, pytest.suffix)
    assert isinstance(result, pl.DataFrame), "Result should be a Polars DataFrame"
    assert result.shape == (12, 4), "Result should have 12 rows and 4 columns"
    assert result.columns == ["sample", "process", "metric", "value"], "Columns should match expected names"
    assert result["sample"].unique().to_list() == [pytest.sample_id], "Sample ID should match the expected value"
    assert result["metric"].unique(maintain_order=True).to_list() == ["input", "output"], "Metrics should match expected values"


def test_sum_counts(constants):
    """
    Test the sum_counts function to ensure it sums the counts correctly for a specific process and metric.
    """
    df = mod.collect_counts(pytest.input_dir, pytest.suffix)
    result = mod.sum_counts(df, pytest.process, pytest.metric)
    assert isinstance(result, int), "Result should be an integer"
    assert result == 49935, "Sum of counts for 'assemble_fragments' process and 'input' metric should be 49935"


def test_summarize_counts(constants):
    """
    Test the summarize_counts function to ensure it summarizes read counts correctly by sample ID.
    """
    all_counts = mod.collect_counts(pytest.input_dir, pytest.suffix)
    result = mod.summarize_counts(all_counts, pytest.sample_id)
    
    assert isinstance(result, pl.DataFrame), "Result should be a Polars DataFrame"
    assert result.shape == (6, 4), "Result should have 6 rows and 4 columns"
    assert result.columns == ["sample", "process", "metric", "value"], "Columns should match expected names"
    
    # Check if the sample ID is correct
    assert result["sample"].unique().to_list() == [pytest.sample_id], "Sample ID should match the expected value"
    
    # Check if the processes and metrics are correctly summarized
    expected_processes = ["assemble_fragments", "reannotate_fragments", "reannotate_bam"]
    expected_metrics = ["input", "output"]
    
    assert set(result["process"].unique().to_list()) == set(expected_processes), "Processes should match expected values"
    assert set(result["metric"].unique(maintain_order=True).to_list()) == set(expected_metrics), "Metrics should match expected values"