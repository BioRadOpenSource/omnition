# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import polars as pl

# load the module as a spec
spec = importlib.util.spec_from_file_location("merging", "bin/atacCallDoublets.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

@pytest.fixture
def constants():
    pytest.quant_file = "test/unit/omnition-core/files/merged.barcodeQuantSimple.csv"
    pytest.allowlist_file = "test/unit/omnition-core/files/merged_barcode_allowlist.csv"
    pytest.ti_len = 6
    pytest.csv_dir = "test/unit/omnition-core/files/"
    pytest.name = "DemoAtacCombinatorial_S1"
    pytest.one_to_one = False
    pytest.barcoded_tn5 = True

def test_load_overlap_df(constants):
    overlap_df = mod.load_overlap_df(pytest.csv_dir)
    assert isinstance(overlap_df, pl.DataFrame)
    assert overlap_df.shape[0] == 11
    assert overlap_df.shape[1] == 3

def test_substr_right(constants):
    overlap_df = mod.load_overlap_df(pytest.csv_dir)
    overlap_df = mod.substr_right(overlap_df, "barc1", int(pytest.ti_len))
    assert overlap_df.shape[0] == 11
    assert overlap_df.shape[1] == 4
    overlap_df = mod.substr_right(overlap_df, "barc2", int(pytest.ti_len))
    assert overlap_df.shape[0] == 11
    assert overlap_df.shape[1] == 5

def test_read_hq_bc_file(constants):
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    assert allowlist_bc.shape[0] == 2000
    assert allowlist_bc.shape[1] == 1

def test_read_n_bc_file(constants):
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    quantification_df = mod.read_n_bc_file(pytest.quant_file, allowlist_bc['bc'].to_list())
    assert quantification_df.shape[0] == 2000
    assert quantification_df.shape[1] == 2

def test_create_count_df(constants):
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    quantification_df = mod.read_n_bc_file(pytest.quant_file, allowlist_bc['bc'].to_list())
    count_df = mod.create_count_df(quantification_df)
    assert count_df.shape[0] == 2000
    assert count_df.shape[1] == 2

def test_create_implicated_df(constants):
    overlap_df = mod.load_overlap_df(pytest.csv_dir)
    overlap_df = mod.substr_right(overlap_df, "barc1", int(pytest.ti_len))
    overlap_df = mod.substr_right(overlap_df, "barc2", int(pytest.ti_len))
    overlap_df = overlap_df.filter(pl.col("match_barc1") == pl.col("match_barc2")).select([pl.col("barc1"), pl.col("barc2"), pl.col("n_both")])
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    quantification_df = mod.read_n_bc_file(pytest.quant_file, allowlist_bc['bc'].to_list())
    count_df = mod.create_count_df(quantification_df)
    implicated_df = mod.create_implicated_df(overlap_df, count_df)
    assert implicated_df.shape[0] == 11
    assert implicated_df.shape[1] == 6


def test_get_density_threshold(constants):
    overlap_df = mod.load_overlap_df(pytest.csv_dir)
    overlap_df = mod.substr_right(overlap_df, "barc1", int(pytest.ti_len))
    overlap_df = mod.substr_right(overlap_df, "barc2", int(pytest.ti_len))
    overlap_df = overlap_df.filter(pl.col("match_barc1") == pl.col("match_barc2")).select([pl.col("barc1"), pl.col("barc2"), pl.col("n_both")])
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    quantification_df = mod.read_n_bc_file(pytest.quant_file, allowlist_bc['bc'].to_list())
    count_df = mod.create_count_df(quantification_df)
    implicated_df = mod.create_implicated_df(overlap_df, count_df)
    jaccard_results = mod.get_density_threshold(implicated_df.select([pl.col("jaccard_frag")]).head(1000000),"jaccard", logTransform = True)
    assert jaccard_results[0] == 0.019485749533907873
    assert jaccard_results[1] == 0.019485749533907873

def test_group_beads(constants):
    overlap_df = mod.load_overlap_df(pytest.csv_dir)
    overlap_df = mod.substr_right(overlap_df, "barc1", int(pytest.ti_len))
    overlap_df = mod.substr_right(overlap_df, "barc2", int(pytest.ti_len))
    overlap_df = overlap_df.filter(pl.col("match_barc1") == pl.col("match_barc2")).select([pl.col("barc1"), pl.col("barc2"), pl.col("n_both")])
    allowlist_bc = mod.read_hq_bc_file(pytest.allowlist_file)
    quantification_df = mod.read_n_bc_file(pytest.quant_file, allowlist_bc['bc'].to_list())
    count_df = mod.create_count_df(quantification_df)
    implicated_df = mod.create_implicated_df(overlap_df, count_df)
    jaccard_results = mod.get_density_threshold(implicated_df.select([pl.col("jaccard_frag")]).head(1000000),"jaccard", logTransform = True)
    jaccard_min_frag = jaccard_results[0]
    implicated_df = implicated_df.with_column(pl.col("jaccard_frag")
                        .apply(lambda x: x > jaccard_min_frag).alias("merged"))
    barcode_filtered_df = implicated_df.filter(pl.col("merged"))
    barcode_translate_df = quantification_df.with_column(pl.lit("").alias("droplet_barcode"))
    barcode_translate_df = mod.group_beads(quantification_df,
                            barcode_filtered_df, barcode_translate_df, pytest.one_to_one, pytest.name, int(pytest.ti_len), pytest.barcoded_tn5)
    assert barcode_translate_df.shape[0] == 2000
    assert barcode_translate_df.shape[1] == 3
