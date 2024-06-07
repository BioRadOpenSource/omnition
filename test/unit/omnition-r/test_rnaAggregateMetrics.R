# unit tests to run for rnaAggregateMetrics.R

# script is an executable, read lines and include only dependency and functions
lines <- scan("../../../bin/rnaAggregateMetrics.R",
    what = character(), sep = "\n")

# find start of library calls
library_start <- grep("# Loading dependencies", lines)

# find start of executable
exe_start <- grep("# Analysis", lines)

# trim down sourced file to include only loading dependencies and functions
lines <- lines[c(library_start:exe_start)]

# source the trimmed file
source(textConnection(paste(lines, collapse = "\n")))

# begin testing

metric_files_dir <- "files/metric-files/"

#test aggregate read counts
read_counts_df <- aggregate_read_counts(metric_files_dir)

fastq_counts_df <- read_counts_df %>%
    filter(process == "fastqc")

test_species <- "DemoRnaMerge_S1"
test_value <- 97000

test_that("read count aggregation is correct",{
    #test read counts were created correctly
     fastq_counts_df %>%
        is.data.frame(.) %>%
        expect_true(.)  
    #read counts are appended by read number
    expect_true("read_count_r2" %in% fastq_counts_df$metric)
    #read is not a column in the df
    expect_false("read" %in% colnames(fastq_counts_df))
    #test data pulled is expected and correct type
    expect_is(fastq_counts_df$sample,'character')
    expect_is(fastq_counts_df$process,'character')
    expect_is(fastq_counts_df$metric,'character')
    expect_is(fastq_counts_df$value,'numeric')
    expect_equal(fastq_counts_df$value[1],test_value)
    expect_equal(fastq_counts_df$value[2],test_value)
    expect_equal(fastq_counts_df$sample[1],test_species)
    expect_equal(fastq_counts_df$sample[2],test_species)
})

downsampled_DOs_df <- read_csv("files/metric-files/DemoRnaDropletMixed_S1_downsample_deconvolution_oligos_read_counts.csv")
test_value_species <- "DemoRnaDropletMixed_S1"
debarcoder_pct <- 80
cutadapt_trim_pct <- 80
test_df <- data.frame(matrix(ncol = 5, nrow = 0))

test_that("estimation of downsampled DOs inputs is correct", {
    downsampled_DOs_df %>%
    is.data.frame(.) %>%
    expect_true(.)

    #test back estimated counts
    test_df <- estimate_downsampled_DOS(downsampled_DOs_df, "DemoRnaDropletMixed_S1",debarcoder_pct, test_df)
    print(test_df$value[1])
    expect_equal(as.integer(test_df$value[1]), 9910)

    #test data pulled is expected and correct type
    expect_is(downsampled_DOs_df$sample,'character')
    expect_is(downsampled_DOs_df$process,'character')
    expect_is(downsampled_DOs_df$metric,'character')
    expect_is(downsampled_DOs_df$value,'numeric')
})

#test estimating downsampled and total data
downsampled_cDNA_df <- read_csv("files/metric-files/DemoRnaDropletMixed_S1_downsample_cDNA_read_counts.csv")
test_value_species <- "DemoRnaDropletMixed_S1"
debarcoder_pct <- 80
cutadapt_trim_pct <- 80
test_df <- data.frame(matrix(ncol = 5, nrow = 0))

test_that("estimation of downsampled cDNA inputs is correct", {
    downsampled_cDNA_df %>%
    is.data.frame(.) %>%
    expect_true(.)

    print(downsampled_cDNA_df$process)

    #test back estimated counts
    test_df <- estimate_downsampled_cDNA(downsampled_cDNA_df, "DemoRnaDropletMixed_S1",debarcoder_pct,cutadapt_trim_pct,test_df)
    expect_equal(as.integer(test_df$value[1]), 126464)  

    #test data pulled is expected and correct type
    expect_is(downsampled_cDNA_df$sample,'character')
    expect_is(downsampled_cDNA_df$process,'character')
    expect_is(downsampled_cDNA_df$metric,'character')
    expect_is(downsampled_cDNA_df$value,'numeric')  
})

#test mixed species counts
mixed_species_counts_df <- read_mixed_species_files(metric_files_dir)

test_value_mixed_cells <- 647.0
test_value_observed_mixed_species_multiplets <- 51.9
test_value_estimate_total_multiples <- 104.0
test_cell_purity <- 98.4

test_that("mixed species count is correct", {
    #test mixed species counts were created correctly
    mixed_species_counts_df %>%
        is.data.frame(.) %>%
        expect_true(.)
    #test crosstalk is process name
    expect_true("crosstalk" %in% mixed_species_counts_df$process)
    #test data pulled is expected and correct type
    expect_is(mixed_species_counts_df$sample,'factor')
    expect_is(mixed_species_counts_df$process,'factor')
    expect_is(mixed_species_counts_df$metric,'factor')
    expect_is(mixed_species_counts_df$value,'numeric')  
    expect_equal(mixed_species_counts_df$value[1],test_value_mixed_cells)
    expect_equal(mixed_species_counts_df$value[2],test_value_observed_mixed_species_multiplets)
    expect_equal(mixed_species_counts_df$value[3],test_value_estimate_total_multiples)
    expect_equal(mixed_species_counts_df$value[4],test_cell_purity)
})

#test read summarize expression

expression_files_dir <- "files/test-summarize-expression-files"

expr_df <- read_summarize_expression_files(expression_files_dir)

s1_expr <- read_csv("files/test-summarize-expression-files/Sample1.final.scrnaseq_counts.csv")
s2_expr <- read_csv("files/test-summarize-expression-files/Sample2.final.scrnaseq_counts.csv")
s1_allow <- read_csv("files/test-summarize-expression-files/Sample1_barcode_allowlist.csv", col_names = F)
s2_allow <- read_csv("files/test-summarize-expression-files/Sample2_barcode_allowlist.csv", col_names = F)

test_that("expression file aggregation is correct", {
    # test that output is correct
    expr_df %>%
        is_tibble(.) %>%
        expect_true(.)
    # confirm that counts are correct; should include all cells
    s1_summary <- expr_df %>% filter(sample=="Sample1")
    expect_true(s1_summary %>%
        filter(metric == "deduplicated_reads") %>%
        select(value) == sum(s1_expr$input_reads))
    expect_true(s1_summary %>%
        filter(metric == "umi") %>%
        select(value) == sum(s1_expr$umi))
    expect_true(s1_summary %>%
        filter(metric == "genic_reads") %>%
        select(value) == sum(s1_expr$genic_reads))
    # not repeating this for every single feature
    s2_summary <- expr_df %>% filter(sample=="Sample2")
    expect_true(s2_summary %>%
        filter(metric == "deduplicated_reads") %>%
        select(value) == sum(s2_expr$input_reads))
    expect_true(s2_summary %>%
        filter(metric == "umi") %>%
        select(value) == sum(s2_expr$umi))
    expect_true(s2_summary %>%
        filter(metric == "genic_reads") %>%
        select(value) == sum(s2_expr$genic_reads))
    # check to ensure that reads_in_cells has excluded the cells not on allowlist
    s1_in_cells <- s1_expr %>%
        filter(barcode %in% s1_allow$X1) %>%
        summarize(frac_in_cells = sum(input_reads))
    s1_in_cells <- (s1_in_cells / sum(s1_expr$input_reads)) * 100
    expect_true(
        s1_in_cells == expr_df %>%
            filter(sample=="Sample1" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
    s2_in_cells <- s2_expr %>%
        filter(barcode %in% s2_allow$X1) %>%
        summarize(frac_in_cells = sum(input_reads))
    s2_in_cells <- (s2_in_cells / sum(s2_expr$input_reads)) * 100
    expect_true(
        s2_in_cells == expr_df %>%
            filter(sample=="Sample2" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
    # check species specific counts-per-cell
    # naming
    expect_equal(c("Sample1_DemoReferenceRnaMixedHomoSapiens_S1",
     "Sample2_DemoReferenceRnaMixedHomoSapiens_S1",
     "Sample1_DemoReferenceRnaMixedMusMusculus_S2",
     "Sample2_DemoReferenceRnaMixedMusMusculus_S2")
     %in% unique(expr_df$sample), c(TRUE,TRUE,TRUE,TRUE))
    s1_species1 <- s1_expr %>%
                filter(barcode %in% s1_allow$X1) %>%
                summarize(frac_human = sum(input_DemoReferenceRnaMixedHomoSapiens_S1_reads))
    s1_species1_tot <- s1_expr %>%
                summarize(all_human = sum(input_DemoReferenceRnaMixedHomoSapiens_S1_reads))
    s1_species2 <- s1_expr %>%
                filter(barcode %in% s1_allow$X1) %>%
                summarize(frac_mouse = sum(input_DemoReferenceRnaMixedMusMusculus_S2_reads))
    s1_species2_tot <- s1_expr %>%
                summarize(all_mouse = sum(input_DemoReferenceRnaMixedMusMusculus_S2_reads))
    expect_true(
        (s1_species1 / s1_species1_tot)*100 == expr_df %>%
            filter(sample=="Sample1_DemoReferenceRnaMixedHomoSapiens_S1" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
    expect_true(
        (s1_species2 / s1_species2_tot)*100 == expr_df %>%
            filter(sample=="Sample1_DemoReferenceRnaMixedMusMusculus_S2" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
    s2_species1 <- s2_expr %>%
                filter(barcode %in% s2_allow$X1) %>%
                summarize(frac_human = sum(input_DemoReferenceRnaMixedHomoSapiens_S1_reads))
    s2_species1_tot <- s2_expr %>%
                summarize(all_human = sum(input_DemoReferenceRnaMixedHomoSapiens_S1_reads))
    s2_species2 <- s2_expr %>%
                filter(barcode %in% s2_allow$X1) %>%
                summarize(frac_mouse = sum(input_DemoReferenceRnaMixedMusMusculus_S2_reads))
    s2_species2_tot <- s2_expr %>%
                summarize(all_mouse = sum(input_DemoReferenceRnaMixedMusMusculus_S2_reads))
    expect_true(
        (s2_species1 / s2_species1_tot)*100 == expr_df %>%
            filter(sample=="Sample2_DemoReferenceRnaMixedHomoSapiens_S1" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
    expect_true(
        (s2_species2 / s2_species2_tot)*100 == expr_df %>%
            filter(sample=="Sample2_DemoReferenceRnaMixedMusMusculus_S2" & metric=="fraction_reads_in_cells") %>%
            select(value)
    )
})
