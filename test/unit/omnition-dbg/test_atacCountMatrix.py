# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import os
import ray
import pysam
import shutil
import pandas as pd
from scipy.sparse import lil_matrix, spmatrix, hstack

# load the module as a spec
spec = importlib.util.spec_from_file_location("dedup", "bin/atacCountMatrix.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# create parallel pool for ray
ray.init(num_cpus=1, ignore_reinit_error = True)

# execute test functions

def test_frip_is_not_zero():
    bam_file = "test/unit/omnition-dbg/files/DemoAtacNormal_S1.final.bam"
    peak_bed = "test/unit/omnition-dbg/files/DemoAtacNormal_S1.fixedwidthpeaks.bed"
    test_allowlist = "test/unit/omnition-dbg/files/DemoAtacNormal_S1.QCstats.csv"
    contigs_to_use = ["DemoRefAtachg38.10", "DemoRefAtacmm10.10"]
    chromvar_compat = False
    verbose = False

    wl = pd.read_csv(test_allowlist, delimiter = ",")
    goodbc = wl.DropBarcode

    # Process contigs in parallel
    results = ray.get(
        [
            mod.detect_overlaps.remote(
                bam_file, contig, chromvar_compat, goodbc, peak_bed, verbose
            )
            for contig in contigs_to_use
        ]
    )

    # Merge all the cell metadata
    meta = mod.list_extract(results, 1)
    merged_dict = meta[0]
    for i in range(1, len(meta)):
        for barcode in goodbc:
            merged_dict[barcode] += meta[i][barcode]

    # Create data frame for storing FRIP
    frips = pd.DataFrame({'DropBarcode': pd.Series(dtype='str'),
                          'FRIP': pd.Series(dtype='float')})

    for cellid, meta_data in merged_dict.items():
        # calculate per-cell frip
        frip = round(meta_data.reads_in_peaks/meta_data.total_reads,3)
        # append to FRIP store
        frips = frips.append(pd.DataFrame({"DropBarcode": [cellid], "FRIP": [frip]}))

    # Add FRIP-per-cell to metadata
    wl_out = mod.format_metadata(wl, frips)

    # Check that data frame has proper column names
    expected_columns = ['DropBarcode', 'totalNuclearFrags', 'uniqueNuclearFrags',
       'totalMitoFrags', 'uniqueMitoFrags', 'totalDemoRefAtachg38Frags',
       'uniqueDemoRefAtachg38Frags', 'totalDemoRefAtacmm10Frags',
       'uniqueDemoRefAtacmm10Frags', 'duplicateProportion', 'librarySize',
       'meanInsertSize', 'medianInsertSize', 'tssProportion', 'FRIP',
       'beadsInDrop']

    assert(list(wl_out.columns) == expected_columns)

    # Check that FRIP column isn't just zeros
    assert(all(x > 0 for x in list(wl_out['FRIP'])))

