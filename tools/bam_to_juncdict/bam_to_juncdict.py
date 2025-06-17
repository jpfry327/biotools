#!/usr/bin/env python3

import pysam
import pickle
import gzip
import argparse
from datetime import datetime

def timestamp():
   return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def main():
   parser = argparse.ArgumentParser(description='Extract splice junctions from BAM file')
   parser.add_argument('-i', '--input', required=True, help='Input BAM file path')
   parser.add_argument('-o', '--output', required=True, help='Output pickle file path')
   parser.add_argument('--mapq', type=int, default=50, help='Minimum mapping quality (default: 50)')
   
   args = parser.parse_args()
   
   print(f"[{timestamp()}] Processing BAM file: {args.input}")
   print(f"[{timestamp()}] Output file: {args.output}")
   print(f"[{timestamp()}] Minimum MAPQ: {args.mapq}")
   
   bam_file = pysam.AlignmentFile(args.input, 'r')
   
   introns = dict()
   
   chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
   
   for chrm in chroms:
       print(f'[{timestamp()}] Processing {chrm}')
       introns[chrm] = bam_file.find_introns(read for read in bam_file.fetch(chrm) if read.mapping_quality >= args.mapq)
   
   bam_file.close()
   
   print(f"[{timestamp()}] Saving to {args.output}")
   with gzip.open(args.output, 'wb') as f:
       pickle.dump(introns, f, protocol=pickle.HIGHEST_PROTOCOL)
   
   print(f"[{timestamp()}] Done!")

if __name__ == "__main__":
   main()