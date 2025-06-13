#!/bin/bash

# FitHiChIP to UCSC Interact Format Converter
# Converts FitHiChIP interaction BED files to UCSC Interact format

# Usage function
usage() {
    echo "Usage: $0 -i <input_file> [-o <output_file>] [-e <experiment_name>] [-s <score_column>] [-v <value_column>] [-h]"
    echo ""
    echo "Options:"
    echo "  -i <input_file>      Input FitHiChIP interaction BED file (required)"
    echo "  -o <output_file>     Output UCSC Interact file (default: input_file.interact.bed)"
    echo "  -e <experiment_name> Experiment name (default: FitHiChIP)"
    echo "  -s <score_column>    Column number for scoring (default: 7, which is contact count 'cc')"
    echo "  -v <value_column>    Column number for interaction value (default: 7)"
    echo "  -h                   Show this help message"
    echo ""
    echo "Expected FitHiChIP format:"
    echo "chr1 s1 e1 chr2 s2 e2 cc Coverage1 isPeak1 Bias1 Mapp1 GCContent1 RESites1 Coverage2 isPeak2 Bias2 Mapp2 GCContent2 RESites2 p exp_cc_Bias p_Bias dbinom_Bias P-Value_Bias Q-Value_Bias"
    echo ""
    echo "Output UCSC Interact format (18 columns):"
    echo "chrom chromStart chromEnd name score value exp color sourceChrom sourceStart sourceEnd sourceName sourceStrand targetChrom targetStart targetEnd targetName targetStrand"
}

# Default values
INPUT_FILE=""
OUTPUT_FILE=""
EXPERIMENT_NAME="FitHiChIP"
SCORE_COLUMN=7
VALUE_COLUMN=7

# Parse command line arguments
while getopts "i:o:e:s:v:h" opt; do
    case $opt in
        i)
            INPUT_FILE="$OPTARG"
            ;;
        o)
            OUTPUT_FILE="$OPTARG"
            ;;
        e)
            EXPERIMENT_NAME="$OPTARG"
            ;;
        s)
            SCORE_COLUMN="$OPTARG"
            ;;
        v)
            VALUE_COLUMN="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done

# Check if input file is provided
if [ -z "$INPUT_FILE" ]; then
    echo "Error: Input file is required"
    usage
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' does not exist"
    exit 1
fi

# Set default output file if not provided
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="${INPUT_FILE%.bed}.interact.bed"
fi

echo "=== FitHiChIP to UCSC Interact Format Converter ==="
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Experiment name: $EXPERIMENT_NAME"
echo "Score column: $SCORE_COLUMN"
echo "Value column: $VALUE_COLUMN"
echo ""

# First, get the maximum value for score scaling
echo "Calculating score scaling factor..."
MAX_SCORE=$(tail -n +2 "$INPUT_FILE" | awk -v col="$SCORE_COLUMN" '{if ($col > max) max = $col} END {print max}')

if [ -z "$MAX_SCORE" ] || [ "$MAX_SCORE" = "0" ]; then
    echo "Warning: Could not determine maximum score or maximum is 0. Using 1000 as default."
    MAX_SCORE=1000
fi

echo "Maximum score value: $MAX_SCORE"

# Count input lines
INPUT_LINES=$(tail -n +2 "$INPUT_FILE" | wc -l)
echo "Processing $INPUT_LINES interactions..."

# Convert FitHiChIP to UCSC Interact format
# Skip header line and process data
tail -n +2 "$INPUT_FILE" | awk -v max_score="$MAX_SCORE" \
                              -v exp_name="$EXPERIMENT_NAME" \
                              -v score_col="$SCORE_COLUMN" \
                              -v value_col="$VALUE_COLUMN" '
BEGIN {
    OFS = "\t"
}
{
    # Input columns: chr1 s1 e1 chr2 s2 e2 cc ...
    chr1 = $1
    s1 = $2
    e1 = $3
    chr2 = $4
    s2 = $5
    e2 = $6
    
    # Get score and value from specified columns
    score_val = $score_col
    value_val = $value_col
    
    # Check if intrachromosomal
    is_intrachrom = (chr1 == chr2)
    
    # UCSC Interact format fields:
    
    # 1. chrom - For intrachromosomal use the chromosome, for interchromosomal use first chrom
    chrom = chr1
    
    # 2. chromStart - Start of lower region (for intrachrom) or this region (for interchrom)
    if (is_intrachrom) {
        chromStart = (s1 < s2) ? s1 : s2
    } else {
        chromStart = s1
    }
    
    # 3. chromEnd - End of upper region (for intrachrom) or this region (for interchrom)  
    if (is_intrachrom) {
        chromEnd = (e1 > e2) ? e1 : e2
    } else {
        chromEnd = e1
    }
    
    # 4. name - Interaction name
    name = chr1 ":" s1 "-" e1 "/" chr2 ":" s2 "-" e2
    
    # 5. score - Scale to 0-1000 range
    score = int((score_val / max_score) * 1000)
    if (score > 1000) score = 1000
    if (score < 0) score = 0
    
    # 6. value - Interaction strength
    value = value_val
    
    # 7. exp - Experiment name
    experiment = exp_name
    
    # 8. color - Default color
    color = "0"
    
    # 9. sourceChrom - First region chromosome
    sourceChrom = chr1
    
    # 10. sourceStart - First region start
    sourceStart = s1
    
    # 11. sourceEnd - First region end
    sourceEnd = e1
    
    # 12. sourceName - First region name
    sourceName = chr1 ":" s1 "-" e1
    
    # 13. sourceStrand - Not available in FitHiChIP
    sourceStrand = "."
    
    # 14. targetChrom - Second region chromosome
    targetChrom = chr2
    
    # 15. targetStart - Second region start
    targetStart = s2
    
    # 16. targetEnd - Second region end
    targetEnd = e2
    
    # 17. targetName - Second region name
    targetName = chr2 ":" s2 "-" e2
    
    # 18. targetStrand - Not available in FitHiChIP
    targetStrand = "."
    
    # Output UCSC Interact format
    print chrom, chromStart, chromEnd, name, score, value, experiment, color, \
          sourceChrom, sourceStart, sourceEnd, sourceName, sourceStrand, \
          targetChrom, targetStart, targetEnd, targetName, targetStrand
}' > "$OUTPUT_FILE"

# Check if conversion was successful
if [ $? -eq 0 ]; then
    OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
    echo ""
    echo "Conversion completed successfully!"
    echo "Output lines: $OUTPUT_LINES"
    echo "Output file: $OUTPUT_FILE"
    
    # Show first few lines as example
    echo ""
    echo "First 3 lines of output:"
    head -n 3 "$OUTPUT_FILE" | cat -n
    
    # Show some statistics
    echo ""
    echo "Statistics:"
    echo "Score range: $(awk '{print $5}' "$OUTPUT_FILE" | sort -n | head -n 1) - $(awk '{print $5}' "$OUTPUT_FILE" | sort -n | tail -n 1)"
    echo "Value range: $(awk '{print $6}' "$OUTPUT_FILE" | sort -n | head -n 1) - $(awk '{print $6}' "$OUTPUT_FILE" | sort -n | tail -n 1)"
    
    # Count chromosomes
    echo "Chromosomes involved:"
    awk '{print $1}' "$OUTPUT_FILE" | sort | uniq -c | head -n 10
    
else
    echo "Error: Conversion failed"
    exit 1
fi

echo ""
echo "=== Conversion Complete ==="
