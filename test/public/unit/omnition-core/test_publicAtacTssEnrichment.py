# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest

# load the module as a spec
spec = importlib.util.spec_from_file_location("enrich", "bin/publicAtacTssEnrichment.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.fixture
def constants():
    pytest.tss_count_matrix = "test/public/unit/omnition-core/files/test_publicAtacTssEnrichment_matrix.gz"
    pytest.tss_window_size = 2001


def test_calculate_total_tss_coverage(constants):
    count, rolling_total = mod.calculate_total_tss_coverage(
        pytest.tss_count_matrix, pytest.tss_window_size
    )
    assert count == 7522
    assert len(rolling_total) == 2001


def test_normalize_tss_coverage(constants):
    count, rolling_total = mod.calculate_total_tss_coverage(
        pytest.tss_count_matrix, pytest.tss_window_size
    )
    norm, tss_enrichment_score = mod.normalize_tss_coverage(
        rolling_total, count, pytest.tss_window_size
    )
    assert len(norm) == 2001
    assert tss_enrichment_score == 7.2
