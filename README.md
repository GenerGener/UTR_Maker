# UTR_maker
Make 5' and 3' UTRs from annotated cDNA mRNA transcript models.

# Tutorial

This tutorial demonstrates how to use the UTR Maker script "utr_maker.py" to analyze GenBank files and extract 5' and 3' UTR regions from mRNA transcript models. We'll use two example files: MZ242719.txt and MZ242720.txt.

## Prerequisites

Before using the script, ensure you have the following installed:
- Python 3.6 or later
- BioPython library (`pip install biopython`)

## Script Overview

The UTR Maker script provides a class-based approach to:
1. Read GenBank files
2. Identify locus segments
3. Find coding sequence boundaries
4. Extract 5' and 3' UTRs
5. Save the results

## Installation

1. Save the script as `utr_maker.py`
2. Place your GenBank files in the same directory or provide the full path

## Basic Usage

```python
from utr_maker import UTRMaker

# Create UTRMaker instance
maker = UTRMaker("MZ242719.txt")

# Save UTRs to files
maker.save_utrs("MZ242719")

# Get detailed information
details = maker.get_utr_details()
```

## Example Analysis

Let's analyze our two example files:

### Example 1: MZ242719.txt

This file contains a complete HIV-1 genome. The UTRs are structured as follows:

5' UTR components:
- Locus segment 1 (positions 1-290): Minimal 5' UTR
- Part of segment 2 before gag starts (positions 291-336)
- Total length: 336 bp

3' UTR components:
- Sequences after nef gene (positions 8955-9173)
- Total length: 219 bp

```python
maker = UTRMaker("MZ242719.txt")
details = maker.get_utr_details()
print(f"5' UTR length: {details['5_utr_length']} bp")
print(f"3' UTR length: {details['3_utr_length']} bp")
```

### Example 2: MZ242720.txt

This file represents a different transcript model. The UTRs are structured as:

5' UTR components:
- Locus segment 1 (positions 1-290): Minimal 5' UTR
- Locus segment 9 (positions 291-359): tat/rev segment
- Part of segment 10 before coding starts (positions 360-375)
- Total length: 375 bp

3' UTR components:
- Sequences after the last coding sequence (positions 3723-3941)
- Total length: 219 bp

```python
maker = UTRMaker("MZ242720.txt")
details = maker.get_utr_details()
print(f"5' UTR length: {details['5_utr_length']} bp")
print(f"3' UTR length: {details['3_utr_length']} bp")
```

## Output Files

For each input file, the script generates two FASTA files:
- `{filename}_5UTR.fasta`: Contains the 5' UTR sequence
- `{filename}_3UTR.fasta`: Contains the 3' UTR sequence

## Advanced Features

The script provides several utility methods:

1. `find_segment_boundaries()`: Get all locus segment positions
2. `find_coding_boundaries()`: Find first CDS start and last CDS stop
3. `get_utr_details()`: Get comprehensive information about UTRs and segments

## Error Handling

The script includes basic error handling for:
- Missing files
- Invalid GenBank format
- Missing coding sequences
- Missing segment annotations

## Notes

- The script assumes GenBank format input files
- UTR regions are determined based on coding sequence boundaries
- Locus segments are identified from misc_feature annotations
- All positions are 0-based following Python convention

## Troubleshooting

If you encounter issues:
1. Verify input file format
2. Check file permissions
3. Ensure BioPython is properly installed
4. Verify CDS annotations exist in the GenBank file

For more complex analyses or custom modifications, refer to the script's source code and BioPython documentation.
