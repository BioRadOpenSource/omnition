# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import tempfile
import json
from pathlib import Path

# load the module as a spec
spec = importlib.util.spec_from_file_location("publicCytometryAdtDebarcoderRef", "bin/publicCytometryAdtDebarcoderRef.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.fixture
def constants():
    pytest.basic_csv_content = """CTGGGCAATTACTCG,Hu.CD19
CTCATTGTAACTCCT,Hu.CD3
AAGTTCACTCTTTGC,Hu.CD16
TGTTCCCGCTCAACT,Hu.CD4
TACGCCTATAACTTG,Hu.CD11c"""

    pytest.complex_csv_content = """CTGGGCAATTACTCG,"Hu.CD19, B-cell marker"
CTCATTGTAACTCCT,"Hu.CD3, T-cell marker"
AAGTTCACTCTTTGC,Hu.CD16
,Empty sequence row
TGTTCCCGCTCAACT,Hu.CD4
  TACGCCTATAACTTG  ,Hu.CD11c with spaces
TCCTTTCCTGATAGG"""

    pytest.single_csv_content = """SINGLE_BARCODE,Hu.CD19"""

    pytest.expected_basic_barcodes = [
        "CTGGGCAATTACTCG",
        "CTCATTGTAACTCCT", 
        "AAGTTCACTCTTTGC",
        "TGTTCCCGCTCAACT",
        "TACGCCTATAACTTG"
    ]

    pytest.expected_complex_barcodes = [
        "CTGGGCAATTACTCG",
        "CTCATTGTAACTCCT",
        "AAGTTCACTCTTTGC",
        "TGTTCCCGCTCAACT",
        "TACGCCTATAACTTG"
    ]

    pytest.expected_single_barcode = ["SINGLE_BARCODE"]

def test_read_csv_barcodes_basic(constants):
    """
    Test the read_csv_barcodes function with a basic CSV file.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(pytest.basic_csv_content)
        csv_file = f.name

    try:
        barcodes = mod.read_csv_barcodes(csv_file)
        assert barcodes == pytest.expected_basic_barcodes, f"Expected {pytest.expected_basic_barcodes}, but got {barcodes}"
        assert len(barcodes) == 5, f"Expected 5 barcodes, but got {len(barcodes)}"
        assert isinstance(barcodes, list), "Expected barcodes to be a list"
        for barcode in barcodes:
            assert isinstance(barcode, str), "Expected each barcode to be a string"
    finally:
        Path(csv_file).unlink()


def test_read_csv_barcodes_complex(constants):
    """
    Test the read_csv_barcodes function with a complex CSV file containing edge cases.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(pytest.complex_csv_content)
        csv_file = f.name

    try:
        barcodes = mod.read_csv_barcodes(csv_file)
        assert barcodes == pytest.expected_complex_barcodes, f"Expected {pytest.expected_complex_barcodes}, but got {barcodes}"
        assert len(barcodes) == 5, f"Expected 5 barcodes (after filtering empty/incomplete rows), but got {len(barcodes)}"

        # Verify no empty strings or whitespace-only sequences
        for barcode in barcodes:
            assert barcode.strip() == barcode, f"Barcode '{barcode}' should not have leading/trailing whitespace"
            assert len(barcode) > 0, "Barcode should not be empty"
    finally:
        Path(csv_file).unlink()


def test_read_csv_barcodes_single(constants):
    """
    Test the read_csv_barcodes function with a single barcode CSV file.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(pytest.single_csv_content)
        csv_file = f.name

    try:
        barcodes = mod.read_csv_barcodes(csv_file)
        assert barcodes == pytest.expected_single_barcode, f"Expected {pytest.expected_single_barcode}, but got {barcodes}"
        assert len(barcodes) == 1, f"Expected 1 barcode, but got {len(barcodes)}"
    finally:
        Path(csv_file).unlink()


def test_read_csv_barcodes_empty_file(constants):
    """
    Test the read_csv_barcodes function with an empty CSV file.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("")
        csv_file = f.name

    try:
        barcodes = mod.read_csv_barcodes(csv_file)
        assert barcodes == [], f"Expected empty list for empty file, but got {barcodes}"
        assert len(barcodes) == 0, f"Expected 0 barcodes for empty file, but got {len(barcodes)}"
    finally:
        Path(csv_file).unlink()


def test_read_csv_barcodes_nonexistent_file(constants):
    """
    Test the read_csv_barcodes function with a non-existent file.
    """
    with pytest.raises(SystemExit):
        mod.read_csv_barcodes("nonexistent_file.csv")


def test_create_json_config_basic(constants):
    """
    Test the create_json_config function with basic barcode list.
    """
    config = mod.create_json_config(pytest.expected_basic_barcodes)

    # Verify overall structure
    assert isinstance(config, dict), "Expected config to be a dictionary"
    assert "primary_fastq" in config, "Expected 'primary_fastq' key in config"
    assert "sequences" in config, "Expected 'sequences' key in config"
    assert config["primary_fastq"] == 0, "Expected primary_fastq to be 0"

    # Verify sequences structure
    sequences = config["sequences"]
    assert isinstance(sequences, list), "Expected sequences to be a list"
    assert len(sequences) == 3, f"Expected 3 sequence sections, but got {len(sequences)}"

    # Verify Barcode section
    barcode_section = sequences[0]
    assert barcode_section["type"] == "Barcode", "Expected first section to be Barcode type"
    assert barcode_section["max_dist"] == 1, "Expected max_dist to be 1 for Barcode section"
    assert barcode_section["values"] == pytest.expected_basic_barcodes, f"Expected barcode values to match input"

    # Verify Umi section
    umi_section = sequences[1]
    assert umi_section["type"] == "Umi", "Expected second section to be Umi type"
    assert umi_section["length"] == 1, "Expected length to be 1 for Umi section"

    # Verify Const section
    const_section = sequences[2]
    assert const_section["type"] == "Const", "Expected third section to be Const type"
    assert const_section["max_dist"] == 1, "Expected max_dist to be 1 for Const section"
    assert const_section["value"] == "AAAAAAAAAA", "Expected const value to be 'AAAAAAAAAA'"


def test_create_json_config_empty(constants):
    """
    Test the create_json_config function with empty barcode list.
    """
    config = mod.create_json_config([])

    assert isinstance(config, dict), "Expected config to be a dictionary"
    assert config["sequences"][0]["values"] == [], "Expected empty values list for empty input"
    assert len(config["sequences"]) == 3, "Expected 3 sections even with empty barcodes"


def test_create_json_config_single(constants):
    """
    Test the create_json_config function with single barcode.
    """
    config = mod.create_json_config(pytest.expected_single_barcode)

    assert isinstance(config, dict), "Expected config to be a dictionary"
    assert config["sequences"][0]["values"] == pytest.expected_single_barcode, "Expected single barcode in values"
    assert len(config["sequences"][0]["values"]) == 1, "Expected exactly one barcode in values"


def test_json_structure_consistency(constants):
    """
    Test that the JSON structure is consistent regardless of input size.
    """
    test_cases = [
        [],
        pytest.expected_single_barcode,
        pytest.expected_basic_barcodes,
        pytest.expected_complex_barcodes
    ]

    for test_barcodes in test_cases:
        config = mod.create_json_config(test_barcodes)

        # All configs should have same structure
        assert "primary_fastq" in config, "Missing primary_fastq key"
        assert "sequences" in config, "Missing sequences key"
        assert len(config["sequences"]) == 3, f"Expected 3 sequences, got {len(config['sequences'])}"

        # Check each section has required keys
        barcode_section = config["sequences"][0]
        assert "type" in barcode_section and barcode_section["type"] == "Barcode"
        assert "max_dist" in barcode_section and barcode_section["max_dist"] == 1
        assert "values" in barcode_section and isinstance(barcode_section["values"], list)

        umi_section = config["sequences"][1]
        assert "type" in umi_section and umi_section["type"] == "Umi"
        assert "length" in umi_section and umi_section["length"] == 1

        const_section = config["sequences"][2]
        assert "type" in const_section and const_section["type"] == "Const"
        assert "max_dist" in const_section and const_section["max_dist"] == 1
        assert "value" in const_section and const_section["value"] == "AAAAAAAAAA"


def test_json_serialization(constants):
    """
    Test that the generated config can be properly serialized to JSON.
    """
    config = mod.create_json_config(pytest.expected_basic_barcodes)

    # Test JSON serialization
    try:
        json_string = json.dumps(config)
        assert isinstance(json_string, str), "Expected JSON serialization to produce string"
    except (TypeError, ValueError) as e:
        pytest.fail(f"JSON serialization failed: {e}")

    # Test JSON deserialization
    try:
        parsed_config = json.loads(json_string)
        assert parsed_config == config, "Expected parsed JSON to match original config"
    except (TypeError, ValueError) as e:
        pytest.fail(f"JSON deserialization failed: {e}")


def test_barcode_validation(constants):
    """
    Test that barcodes are properly validated and cleaned.
    """
    # Create a CSV with problematic data
    csv_content = "  CTGGGCAATTACTCG  ,Hu.CD19\nCTCATTGTAACTCCT,Hu.CD3\n,Empty\n   ,Whitespace\nAAGTTCACTCTTTGC,Hu.CD16"

    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name

    try:
        barcodes = mod.read_csv_barcodes(csv_file)
        # Should extract 3 valid barcodes after cleaning
        expected_clean = ["CTGGGCAATTACTCG", "CTCATTGTAACTCCT", "AAGTTCACTCTTTGC"]
        assert len(barcodes) == 3, f"Expected 3 clean barcodes, got {len(barcodes)}"
        assert barcodes == expected_clean, f"Expected {expected_clean}, got {barcodes}"

        for barcode in barcodes:
            assert barcode.strip() == barcode, f"Barcode should be stripped: '{barcode}'"
            assert len(barcode) > 0, "Barcode should not be empty"
    finally:
        Path(csv_file).unlink()
