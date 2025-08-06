#!/usr/bin/env python3
"""
Script to find evidence of Gly177Asp mutation (rs2305089) in TBXT gene
Position: chr6:166165782 (GRCh38)
Mutation: G>A (GGT>GAT, Gly>Asp)
"""

import pysam
import argparse
from collections import defaultdict

def find_gly177asp_mutation(bam_file, min_mapq=20, min_baseq=20):
    """
    Find evidence of Gly177Asp mutation in BAM file
    
    Args:
        bam_file: Path to BAM file
        min_mapq: Minimum mapping quality
        min_baseq: Minimum base quality
    
    Returns:
        Dictionary with mutation evidence
    """
    
    # rs2305089 coordinates (GRCh38)
    chrom = "chr6"
    pos = 166165782  # 0-based for pysam: 166165781
    
    # Open BAM file
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Count bases at the position
    base_counts = defaultdict(int)
    total_reads = 0
    
    # Iterate over reads covering the position
    for pileupcolumn in bamfile.pileup(chrom, pos-1, pos, stepper="samtools"):
        if pileupcolumn.pos == pos-1:  # pysam uses 0-based coordinates
            
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    
                    read = pileupread.alignment
                    
                    # Filter by mapping quality
                    if read.mapping_quality < min_mapq:
                        continue
                    
                    # Get base and quality at position
                    query_pos = pileupread.query_position
                    if query_pos is not None:
                        base = read.query_sequence[query_pos]
                        qual = read.query_qualities[query_pos]
                        
                        # Filter by base quality
                        if qual < min_baseq:
                            continue
                        
                        # Note: pysam already handles strand orientation for us
                        # The base we get is already relative to the reference strand
                        
                        base_counts[base] += 1
                        total_reads += 1
    
    bamfile.close()
    
    # Calculate results
    # Since TBXT is on negative strand, we look for C>T (complement of G>A)
    ref_count = base_counts.get('C', 0)  # Reference C (complement of G)
    alt_count = base_counts.get('T', 0)  # Alternative T (complement of A)
    
    if total_reads > 0:
        alt_freq = alt_count / total_reads
        ref_freq = ref_count / total_reads
    else:
        alt_freq = 0
        ref_freq = 0
    
    return {
        'position': f"{chrom}:{pos}",
        'total_reads': total_reads,
        'ref_C_count': ref_count,
        'alt_T_count': alt_count,
        'other_bases': {k: v for k, v in base_counts.items() if k not in ['C', 'T']},
        'alt_frequency': alt_freq,
        'ref_frequency': ref_freq,
        'has_mutation': alt_count > 0,
        'likely_genotype': get_genotype(ref_count, alt_count, total_reads)
    }

def get_genotype(ref_count, alt_count, total_reads, het_threshold=0.2):
    """Estimate genotype based on allele frequencies"""
    if total_reads == 0:
        return "No coverage"
    
    alt_freq = alt_count / total_reads
    
    if alt_freq == 0:
        return "0/0 (C/C - wild type)"
    elif alt_freq >= (1 - het_threshold):
        return "1/1 (T/T - homozygous mutation)"
    elif alt_freq >= het_threshold:
        return "0/1 (C/T - heterozygous)"
    else:
        return "Likely 0/0 (low alt frequency)"

def main():
    parser = argparse.ArgumentParser(description="Find Gly177Asp mutation evidence")
    parser.add_argument("bam_files", nargs='+', help="BAM files to analyze")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality")
    parser.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality")
    parser.add_argument("--output", "-o", help="Output TSV file (default: stdout)")
    
    args = parser.parse_args()
    
    import os
    import sys
    
    # Prepare output
    if args.output:
        outfile = open(args.output, 'w')
    else:
        outfile = sys.stdout
    
    # Print header
    header = [
        "sample_name",
        "position", 
        "total_reads",
        "ref_C_count",
        "alt_T_count", 
        "alt_frequency",
        "ref_frequency",
        "genotype",
        "has_mutation",
        "other_bases"
    ]
    outfile.write('\t'.join(header) + '\n')
    
    # Process each BAM file
    for bam_file in args.bam_files:
        try:
            # Get basename without extension
            sample_name = os.path.splitext(os.path.basename(bam_file))[0]
            
            # Analyze the file
            result = find_gly177asp_mutation(bam_file, args.min_mapq, args.min_baseq)
            
            # Format other bases as string
            other_bases_str = ','.join([f"{k}:{v}" for k, v in result['other_bases'].items()]) if result['other_bases'] else "none"
            
            # Write results
            row = [
                sample_name,
                result['position'],
                str(result['total_reads']),
                str(result['ref_C_count']),
                str(result['alt_T_count']),
                f"{result['alt_frequency']:.4f}",
                f"{result['ref_frequency']:.4f}",
                result['likely_genotype'],
                str(result['has_mutation']),
                other_bases_str
            ]
            outfile.write('\t'.join(row) + '\n')
            
            # Progress to stderr
            print(f"Processed: {sample_name}", file=sys.stderr)
            
        except Exception as e:
            print(f"Error processing {bam_file}: {e}", file=sys.stderr)
            # Write error row
            row = [
                os.path.splitext(os.path.basename(bam_file))[0],
                "chr6:166165782",
                "ERROR",
                "ERROR", 
                "ERROR",
                "ERROR",
                "ERROR",
                f"Error: {str(e)}",
                "ERROR",
                "ERROR"
            ]
            outfile.write('\t'.join(row) + '\n')
    
    if args.output:
        outfile.close()
        print(f"Results written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()