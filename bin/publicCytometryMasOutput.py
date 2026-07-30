#!/usr/bin/env python3

import click
import gzip
from pathlib import Path
import scipy.sparse as sp
import scipy.io as sio


def load_barcodes(barcode_file):
    """Load barcodes from a file into a list."""
    with open(barcode_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def load_genes(genes_file):
    """Load genes/features from a file into a list of tuples."""
    with open(genes_file, 'r') as f:
        return [tuple(line.strip().split('\t')) for line in f if line.strip()]


def load_sparse_matrix(matrix_file):
    """Load a sparse matrix from an MTX file (gzipped or not)."""
    open_func = gzip.open if str(matrix_file).endswith('.gz') else open
    with open_func(matrix_file, 'rb') as f:
        return sp.csr_matrix(sio.mmread(f))


def format_features(genes, feature_type):
    """
    Format genes into 3-column feature format.

    Args:
        genes: List of tuples (gene data)
        feature_type: Feature type to add if needed

    Returns:
        List of tuples (gene_id, gene_name, feature_type)
    """
    features = []
    for gene_parts in genes:
        if len(gene_parts) == 1:
            features.append((gene_parts[0], gene_parts[0], feature_type))
        elif len(gene_parts) == 2:
            features.append((gene_parts[0], gene_parts[1], feature_type))
        else:
            features.append((gene_parts[0], gene_parts[1], gene_parts[2]))
    return features


def write_output_files(output_dir, barcodes, features, matrix):
    """Write the three output files (Barcodes, Features, Matrix) compressed."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Write barcodes directly to gzip
    out_barcodes = output_path / "Barcodes.tsv.gz"
    with gzip.open(out_barcodes, 'wt') as f:
        f.write('\n'.join(barcodes) + '\n')
    print(f"  Saved: {out_barcodes}")

    # Write features directly to gzip
    out_features = output_path / "Features.tsv.gz"
    with gzip.open(out_features, 'wt') as f:
        for feat in features:
            f.write(f"{feat[0]}\t{feat[1]}\t{feat[2]}\n")
    print(f"  Saved: {out_features}")

    # Write matrix directly to gzip
    out_matrix = output_path / "Matrix.mtx.gz"
    with gzip.open(out_matrix, 'wb') as f:
        sio.mmwrite(f, matrix, comment='', field='integer')
    print(f"  Saved: {out_matrix}")


def get_file_paths(base_dir, prefix):
    """Get the three file paths for a given prefix."""
    base_path = Path(base_dir)
    return {
        'barcodes': base_path / f"{prefix}.filtered.barcodes.tsv",
        'genes': base_path / f"{prefix}.filtered.genes.tsv",
        'matrix': base_path / f"{prefix}.filtered.mtx.gz"
    }


def merge_matrices(cyto_dir, rna_dir, output_dir, prefix_pattern="*"):
    """
    Merge RNA and cytometry sparse matrices into a single output.
    """
    cyto_path = Path(cyto_dir)
    cyto_barcode_files = list(cyto_path.glob(f"{prefix_pattern}.filtered.barcodes.tsv"))

    if not cyto_barcode_files:
        print(f"No cytometry files matching pattern '{prefix_pattern}.filtered.barcodes.tsv' found in {cyto_dir}")
        return

    for cyto_barcode_file in cyto_barcode_files:
        prefix = cyto_barcode_file.name.replace(".filtered.barcodes.tsv", "")
        print(f"\n{'='*60}")
        print(f"Merging matrices for prefix: {prefix}")
        print(f"{'='*60}")

        # Get file paths
        cyto_files = get_file_paths(cyto_dir, prefix)
        rna_files = get_file_paths(rna_dir, prefix)

        # Check if all required files exist
        required_files = [
            (cyto_files['genes'], "cytometry genes"),
            (cyto_files['matrix'], "cytometry matrix"),
            (rna_files['barcodes'], "RNA barcodes"),
            (rna_files['genes'], "RNA genes"),
            (rna_files['matrix'], "RNA matrix")
        ]

        missing_files = [f"{name}: {file}" for file, name in required_files if not file.exists()]

        if missing_files:
            print(f"  Warning: Missing required files:")
            for mf in missing_files:
                print(f"    - {mf}")
            continue

        try:
            # Load barcodes
            print(f"Loading barcodes...")
            cyto_barcodes = load_barcodes(cyto_files['barcodes'])
            rna_barcodes = load_barcodes(rna_files['barcodes'])

            # Find common barcodes
            common_barcodes = sorted(set(cyto_barcodes) & set(rna_barcodes))

            print(f"  Cytometry barcodes: {len(cyto_barcodes)}")
            print(f"  RNA barcodes: {len(rna_barcodes)}")
            print(f"  Common barcodes: {len(common_barcodes)}")

            if not common_barcodes:
                print(f"  Warning: No common barcodes found!")
                continue

            # Create barcode index mappings
            cyto_bc_to_idx = {bc: i for i, bc in enumerate(cyto_barcodes)}
            rna_bc_to_idx = {bc: i for i, bc in enumerate(rna_barcodes)}

            # Load genes/features
            print(f"Loading features...")
            cyto_genes = load_genes(cyto_files['genes'])
            rna_genes = load_genes(rna_files['genes'])

            print(f"  Cytometry features: {len(cyto_genes)}")
            print(f"  RNA features: {len(rna_genes)}")

            # Load matrices
            print(f"Loading matrices...")
            cyto_matrix = load_sparse_matrix(cyto_files['matrix'])
            rna_matrix = load_sparse_matrix(rna_files['matrix'])

            print(f"  Cytometry matrix shape: {cyto_matrix.shape}")
            print(f"  RNA matrix shape: {rna_matrix.shape}")

            # Subset matrices to common barcodes
            print(f"Subsetting matrices to common barcodes...")
            cyto_indices = [cyto_bc_to_idx[bc] for bc in common_barcodes]
            rna_indices = [rna_bc_to_idx[bc] for bc in common_barcodes]

            cyto_matrix_subset = cyto_matrix[:, cyto_indices]
            rna_matrix_subset = rna_matrix[:, rna_indices]

            # Merge matrices vertically
            print(f"Merging matrices...")
            merged_matrix = sp.vstack([rna_matrix_subset, cyto_matrix_subset])
            print(f"  Merged matrix shape: {merged_matrix.shape}")

            # Format and merge features
            rna_features = format_features(rna_genes, "Gene Expression")
            cyto_features = format_features(cyto_genes, "Antibody Capture")
            merged_features = rna_features + cyto_features

            # Write output files
            print(f"Writing output files...")
            write_output_files(output_dir, common_barcodes, merged_features, merged_matrix)

            print(f"  ✓ Successfully merged {prefix}")

        except Exception as e:
            print(f"  ✗ Error merging {prefix}: {e}")
            import traceback
            traceback.print_exc()


def find_and_convert(input_dir, output_dir, prefix_pattern="*", feature_type="Antibody Capture"):
    """
    Find all matching sparse matrix files and convert them.
    """
    input_path = Path(input_dir)
    barcode_files = list(input_path.glob(f"{prefix_pattern}.filtered.barcodes.tsv"))

    if not barcode_files:
        print(f"No files matching pattern '{prefix_pattern}.filtered.barcodes.tsv' found in {input_dir}")
        return

    for barcode_file in barcode_files:
        prefix = barcode_file.name.replace(".filtered.barcodes.tsv", "")
        print(f"\nProcessing file set with prefix: {prefix}")

        # Get file paths
        files = get_file_paths(input_dir, prefix)

        # Check if all required files exist
        if not files['genes'].exists():
            print(f"  Warning: Missing genes file: {files['genes']}")
            continue
        if not files['matrix'].exists():
            print(f"  Warning: Missing matrix file: {files['matrix']}")
            continue

        try:
            # Load data
            barcodes = load_barcodes(files['barcodes'])
            genes = load_genes(files['genes'])
            matrix = load_sparse_matrix(files['matrix'])

            # Format features
            features = format_features(genes, feature_type)

            # Write output files
            write_output_files(output_dir, barcodes, features, matrix)

            print(f"  ✓ Successfully converted {prefix}")
        except Exception as e:
            print(f"  ✗ Error converting {prefix}: {e}")


@click.command()
@click.option('--input-dir', '-i', default='.',
              help='Input directory containing filtered matrix files')
@click.option('--output-dir', '-o', default='.',
              help='Output directory for converted files')
@click.option('--prefix', '-p', default='*',
              help='File prefix pattern to match (default: * for all files)')
@click.option('--feature-type', '-f', default='Antibody Capture',
              help='Feature type string to add to Features file')
@click.option('--rna-dir', '-r', default=None,
              help='Directory containing RNA matrix files to merge with cytometry data')
def main(input_dir, output_dir, prefix, feature_type, rna_dir):
    """
    Convert sparse matrix files to MAS format.

    Can also merge RNA and cytometry matrices when --rna-dir is specified.
    """
    if rna_dir:
        print("MODE: Merge RNA and Cytometry matrices")
        merge_matrices(
            cyto_dir=input_dir,
            rna_dir=rna_dir,
            output_dir=output_dir,
            prefix_pattern=prefix
        )
    else:
        print("MODE: Convert cytometry matrix only")
        find_and_convert(
            input_dir=input_dir,
            output_dir=output_dir,
            prefix_pattern=prefix,
            feature_type=feature_type
        )


if __name__ == '__main__':
    main()
