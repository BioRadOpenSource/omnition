# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import os
import ray
import pysam
import shutil

# load the module as a spec
spec = importlib.util.spec_from_file_location("summarize", "bin/atacSummarizeAlignments.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# create parallel pool for ray
ray.init(num_cpus=1, ignore_reinit_error = True)

# execute test functions

def test_parse_flagstat():
	# import valid test flagstat file
	flagstat_file = "test/unit/omnition-dbg/files/DemoAtacNormal_S1_flagstat.txt"
	flagstat = mod.parse_flagstat(flagstat_file)

	# check if all expected flagstat metrics are present
	keys=['total', 'secondary', 'supplementary', 'duplicates', 'mapped', 'paired in sequencing',
		'read1', 'read2', 'properly paired', 'with itself and mate mapped', 'singletons',
		'with mate mapped to a different chr', 'with mate mapped to a different chr (mapQ >= 5)']
	all_keys = [key in flagstat.keys() for key in keys]
	assert(all(all_keys))

	# check if all flagstat metrics are integers
	all_types = [isinstance(item, int) for item in flagstat.values()]
	assert(all(all_types))

def test_bam_count():
	# count from test bam file
	test_bam = "test/unit/omnition-dbg/files/feature_aware_test.bam"
	contig = "10"
	result = ray.get(
	    mod.get_stats.remote(
	        test_bam, contig
	    )
	)
	# reads in bam are a duplicated second read in pair
	assert(result == (0,2))

def test_AlignmentStatistics():
	flagstat_file = "test/unit/omnition-dbg/files/feature_aware_flagstat.txt"
	flagstat = mod.parse_flagstat(flagstat_file)

	test_bam = "test/unit/omnition-dbg/files/feature_aware_test.bam"
	contig = "10"
	result = [ray.get(
	    mod.get_stats.remote(
	        test_bam, contig
	    )
	)]
	
	stats = mod.AlignmentStatistics(flagstat, result)

	alignment_dict = stats.stats_dict()

	# check if we've returned a dictionary
	assert(isinstance(alignment_dict, dict))

	# check if our dictionary has the expected keys
	keys = ['CATEGORY', 'TOTAL_READS', 
		'PF_READS_ALIGNED', 'READS_ALIGNED_IN_PAIRS', 
		'PF_READS_IMPROPER_PAIRS']

	all_keys = [key in alignment_dict.keys() for key in keys]

	assert(all(all_keys))

	# check if each value in the dictionary is a list of length 3

	all_length = list()
	for item in alignment_dict.values():
		all_length.append(len(item) == 3)

	assert(all(all_length))




