# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import sys
import polars as pl

# load the module as a spec
sys.path.append("bin")  # add bin directory so that we can import functions
spec = importlib.util.spec_from_file_location("merging", "bin/rnaCorrectCellBarcodes.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# module functions are now in scope

@pytest.fixture
def constants():
    pytest.filtered_bead_file = "test/unit/omnition-core/files/correct_cell_barcode_filtered_bead_file.csv"
    pytest.barcode_translate_file = "test/unit/omnition-core/files/correct_cell_barcode_translate.tsv"
    pytest.barcode_count = 7
    pytest.sampleid = "TestNBC"
    pytest.new_padding = 6

def test_build_barcode_list(constants):
    bead_list = mod.build_barcode_list(pytest.filtered_bead_file, pytest.sampleid)
    assert(len(bead_list) == 7)

def test_format_bead_list(constants):
    bead_list = mod.build_barcode_list(pytest.filtered_bead_file, pytest.sampleid)
    barcode_df = mod.read_barcode_translate(pytest.barcode_translate_file)
    old_barcode_count = barcode_df.select(pl.col("DropBarcode")).n_unique()
    additional_barcode_count = bead_list.select(pl.col("BeadBarcode")).n_unique()
    updated_barcode_count = old_barcode_count + additional_barcode_count
    padding_length = len(str(updated_barcode_count))
    new_barcode_list = mod.format_bead_list(bead_list, old_barcode_count, updated_barcode_count, padding_length)
    new_barcode_list.row(6)
    assert(new_barcode_list.row(0) == tuple(['ATTCCCTATGTGAGCCATAAT', 'TestNBCBC08N1', 1]))
    assert(new_barcode_list.row(6) == tuple(['TTCCCCAAAGTGCGGTCGGTT', 'TestNBCBC14N1', 1]))

def reformat_barcode_translate(constants):
    barcode_translate = mod.read_barcode_translate(pytest.barcode_translate_file)
    corrected_barcode_translate = mod.reformat_barcode_translate(barcode_translate,
                                                               pytest.new_padding)
    assert( corrected_barcode_translate.row(6)[1] == "TestNBCBC0000007N1")
