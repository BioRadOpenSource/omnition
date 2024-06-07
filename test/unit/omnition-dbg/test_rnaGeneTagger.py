# importlib is used to load the script as an arbitrary module
import importlib.util
import pytest
import os
import pysam
import shutil

# load the module as a spec
spec = importlib.util.spec_from_file_location("tag_features", "bin/rnaGeneTagger.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def test_tagging_strategy_exclude_introns():
	include_introns = False
	gene_tag = "XT"
	feature_tag = "XF"
	assignment_tag = "XS"
	bam = pysam.AlignmentFile("test/unit/omnition-dbg/files/test_rna_tag_features.bam")
	Itr = bam.fetch()
	reads = list(Itr)
	for read in reads:
		tag = read.get_tag(gene_tag)
		mod.assign_read_tags(read,
			tag,
			include_introns,
			gene_tag,
			feature_tag,
			assignment_tag
		)
	assert(reads[0].get_tag("XT") == "INTRON")
	assert(reads[0].has_tag("XF") is False)
	assert(reads[1].get_tag("XT") == "GENE1")
	assert(reads[1].get_tag("XF") == "UTR")
	assert(reads[2].get_tag("XT") == "GENE1")
	assert(reads[2].get_tag("XF") == "UTR")
	assert(reads[3].get_tag("XT") == "GENE2")
	assert(reads[3].get_tag("XF") == "EXON")
	assert(reads[4].get_tag("XT") == "GENE1")
	assert(reads[4].get_tag("XF") == "UTR")
	assert(reads[4].get_tag("XN") == 1)
	assert(reads[5].get_tag("XT") == "GENE1")
	assert(reads[5].get_tag("XF") == "EXON")
	assert(reads[6].get_tag("XT") == "GENE1")
	assert(reads[6].get_tag("XF") == "UTR")
	assert(reads[7].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[7].has_tag("XF") is False)
	assert(reads[8].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[8].has_tag("XF") is False)
	assert(reads[9].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[9].has_tag("XF") is False)
	assert(reads[10].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[10].has_tag("XF") is False)
	assert(reads[11].get_tag("XT") == "GENE1")
	assert(reads[11].get_tag("XF") == "UTR")
	assert(reads[12].get_tag("XT") == "GENE1")
	assert(reads[12].get_tag("XF") == "EXON")

def test_tagging_strategy_include_introns():
	include_introns = True
	gene_tag = "XT"
	feature_tag = "XF"
	assignment_tag = "XS"
	bam = pysam.AlignmentFile("test/unit/omnition-dbg/files/test_rna_tag_features.bam")
	Itr = bam.fetch()
	reads = list(Itr)
	for read in reads:
		tag = read.get_tag(gene_tag)
		mod.assign_read_tags(read,
			tag,
			include_introns,
			gene_tag,
			feature_tag,
			assignment_tag
		)
	assert(reads[0].get_tag("XT") == "GENE1")
	assert(reads[0].get_tag("XF") == "INTRON")
	assert(reads[1].get_tag("XT") == "GENE1")
	assert(reads[1].get_tag("XF") == "UTR")
	assert(reads[2].get_tag("XT") == "GENE1")
	assert(reads[2].get_tag("XF") == "UTR")
	assert(reads[3].get_tag("XT") == "GENE2")
	assert(reads[3].get_tag("XF") == "EXON")
	assert(reads[4].get_tag("XT") == "GENE1")
	assert(reads[4].get_tag("XF") == "UTR")
	assert(reads[4].get_tag("XN") == 1)
	assert(reads[5].get_tag("XT") == "GENE1")
	assert(reads[5].get_tag("XF") == "EXON")
	assert(reads[6].get_tag("XT") == "GENE1")
	assert(reads[6].get_tag("XF") == "UTR")
	assert(reads[7].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[7].has_tag("XF") is False)
	assert(reads[8].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[8].has_tag("XF") is False)
	assert(reads[9].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[9].has_tag("XF") is False)
	assert(reads[10].get_tag("XT") == "Unassigned_Ambiguous")
	assert(reads[10].has_tag("XF") is False)
	assert(reads[11].get_tag("XT") == "GENE1")
	assert(reads[11].get_tag("XF") == "UTR")
	assert(reads[12].get_tag("XT") == "GENE1")
	assert(reads[12].get_tag("XF") == "EXON")
