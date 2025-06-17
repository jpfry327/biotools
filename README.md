# biotools
Collection of bioinformatics scripts and utilities

## Tools

### Tractatus
- R package with commonly used functions and tools. Install with `devtools::install_github("jpfry327/biotools", subdir = "tractatus")`. 

### Format Converters
- [`fithichip_to_interact`](format-converters/fithichip_to_interact/) - Convert FitHiChIP BED to UCSC Interact format

### Tools
- [`bam_to_juncdict`](tools/bam_to_juncdict/) - convert BAM file to python dictionary of form `junctions[chrm][(X, Y)]` for junctions at chrm:X-Y
