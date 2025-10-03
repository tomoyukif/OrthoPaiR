# OrthoPaiR Test Suite

This directory contains comprehensive test scripts and subset data for testing OrthoPaiR functions.

## Files Overview

### Subset Data Files
Located in `inst/extdata/subset/`:
- `nb_subset_genome.fa` - NB genome subset (3 sequences, ~500kb each)
- `nb_subset.gff` - NB GFF subset (856 features from 30 genes)
- `nb_subset_cds.fa` - NB CDS subset (50 sequences)
- `nb_subset_prot.fa` - NB protein subset (50 sequences)
- `wk21_subset_genome.fa` - WK21 genome subset (3 sequences, ~500kb each)
- `wk21_subset.gff` - WK21 GFF subset (442 features from 30 genes)
- `wk21_subset_cds.fa` - WK21 CDS subset (50 sequences)
- `wk21_subset_prot.fa` - WK21 protein subset (50 sequences)
- `subset_summary.txt` - Summary of subset data

### Test Scripts

#### 1. `create_subset_data.R`
Creates subset data files from the original input files.
- **Purpose**: Generate small test datasets
- **Usage**: `Rscript create_subset_data.R`
- **Output**: Subset files in `inst/extdata/subset/`

#### 2. `test_core_functions.R`
Tests core OrthoPaiR functions without external dependencies.
- **Purpose**: Test main functionality
- **Usage**: `Rscript test_core_functions.R`
- **Tests**: 14 core functions
- **Output**: Test results and HDF5 file

#### 3. `test_all_functions.R`
Comprehensive test suite for all OrthoPaiR functions.
- **Purpose**: Test all exported functions
- **Usage**: `Rscript test_all_functions.R`
- **Tests**: 23 functions (some may be skipped if dependencies missing)
- **Output**: Complete test results

#### 4. `complete_workflow_example.R`
Demonstrates the complete OrthoPaiR workflow.
- **Purpose**: Show end-to-end usage
- **Usage**: `Rscript complete_workflow_example.R`
- **Output**: Complete analysis results

## Quick Start

### 1. Create Subset Data
```bash
cd /home/ftom/01_wd/softDevel/OrthoPaiR
Rscript create_subset_data.R
```

### 2. Test Core Functions
```bash
Rscript test_core_functions.R
```

### 3. Run Complete Workflow
```bash
Rscript complete_workflow_example.R
```

## Test Results

### Expected Output
- **Core Functions**: 14 tests, ~90% success rate
- **All Functions**: 23 tests, ~70% success rate (some require external tools)
- **Complete Workflow**: Full analysis with results saved

### Output Files
- `test_output/` - Test results and intermediate files
- `orthopair_example_output/` - Complete workflow results
- `*.h5` - HDF5 database files
- `*.csv` - Analysis results in CSV format

## Function Coverage

### Core Functions (Always Tested)
1. `fixInfiles` - Fix input file formats
2. `formatGFF` - Format GFF files
3. `orgInputFiles` - Organize input files
4. `makeOrthoPairDB` - Create database object
5. `rbh` - Reciprocal Best Hits analysis
6. `syntenicOrtho` - Syntenic ortholog analysis
7. `summaryOrthoPair` - Get summary statistics
8. `getOrthoPair` - Get ortholog pairs
9. `getOrphan` - Get orphan genes
10. `makeOrthoGraph` - Create ortholog graph
11. `graph2df` - Convert graph to dataframe
12. `reorgOrthopiars` - Reorganize orthopairs
13. `.getPromoterSeq` - Extract promoter sequences

### Extended Functions (May Require External Tools)
14. `mapProt` - Protein mapping (requires Miniprot)
15. `compareOrthoSeq` - Sequence comparison
16. `orthopair` - Main pipeline function
17. `runSibeliaZ` - SibeliaZ analysis (requires SibeliaZ)
18. `sibeliaRaw2Graph` - Convert SibeliaZ output
19. `getLCBpairs` - Get Locally Collinear Blocks
20. `lcbClassify` - Classify LCBs
21. `plotLCBpairs` - Plot LCB pairs
22. `showLCB` - Show LCB information
23. `statsLCB` - LCB statistics

## Troubleshooting

### Common Issues

1. **File Not Found Errors**
   - Ensure subset data files exist
   - Run `create_subset_data.R` first

2. **External Tool Dependencies**
   - Some functions require Miniprot, SibeliaZ, or BLAST
   - Tests will skip functions if tools are not available

3. **Memory Issues**
   - Subset data is designed to be small
   - If issues persist, reduce subset sizes in `create_subset_data.R`

4. **Permission Errors**
   - Ensure write permissions in output directories
   - Check file paths are accessible

### Debug Mode
Enable verbose output by modifying test scripts:
```r
# Add verbose = TRUE to function calls
rbh(object, n_threads = 1, verbose = TRUE)
```

## Data Description

### NB (Query Species)
- **Source**: Rice cultivar NB
- **Genome**: 3 chromosomes/scaffolds (~500kb each)
- **Genes**: 30 genes with complete annotations
- **Features**: 856 GFF features (genes, transcripts, exons, CDS)

### WK21 (Subject Species)
- **Source**: Rice cultivar WK21
- **Genome**: 3 chromosomes/scaffolds (~500kb each)
- **Genes**: 30 genes with complete annotations
- **Features**: 442 GFF features (genes, transcripts, exons, CDS)

## Performance Notes

- **Subset Size**: Designed for quick testing (~5-10 minutes)
- **Memory Usage**: ~100-200MB peak usage
- **Disk Usage**: ~50MB for all output files
- **CPU Usage**: Single-threaded by default

## Extending Tests

### Adding New Tests
1. Add test function to appropriate test script
2. Follow existing naming convention
3. Include error handling with `tryCatch`
4. Update test counter and results tracking

### Modifying Subset Data
1. Edit parameters in `create_subset_data.R`
2. Adjust `n_genes`, `n_sequences`, `max_length`
3. Re-run subset creation script
4. Update test scripts if needed

## Contact

For issues or questions about the test suite:
- Check function documentation in R files
- Review error messages in test output
- Ensure all dependencies are installed
- Verify file paths and permissions
