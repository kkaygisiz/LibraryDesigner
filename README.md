# LibraryDesigner
Design focused combinatorial peptide libraries from a set of starting sequences

The notebook recombines residues at corresponding positions and selects libraries that maximize starting-sequence coverage under explicit limits on library number and size. All analysis code is contained in one notebook.

## Quick start

Use Python 3.10 or newer. Download this repository and open a terminal in its folder:
```bash
python -m pip install -r requirements.txt
python -m jupyter lab LibraryDesigner.ipynb
```

Set the input path in the settings cell, then run the notebook from top to bottom:

## Input

The simplest CSV contains a sequence column and optional unique identifiers:

```csv
id,sequence
pep_001,ALGS
pep_002,VLAS
```

FASTA, TSV, and one-sequence-per-line text are supported. A `.log` extension is accepted only when the file contains one of those structured formats. 

A CSV may contain `sequence` plus any additional columns. All original columns are preserved in `input_records.csv`; an extra score column does not automatically influence the design. With `WEIGHT_COLUMN = None`, accepted unique sequences have equal priority. Weighting must be selected explicitly with the correct score.

In canonical letter mode, `X` is unknown and is not a concrete residue choice. The default `INVALID_SEQUENCE_POLICY = "skip"` excludes records containing unknown or other invalid symbols, prints a warning and exclusion table, and writes `excluded_sequences.csv`. Set the policy to `"error"` to stop at the first invalid sequence instead. Unknown residues are never guessed or expanded to all amino acids. Coverage is calculated over the accepted unique sequence/group entries, with the original and excluded record counts separately recorded in the metadata.

Optional settings:

| Setting | Meaning |
| --- | --- |
| `WEIGHT_COLUMN` | Nonnegative priorities, such as appropriately transformed enrichment or confidence scores. Default: equal weights. |
| `GROUP_COLUMN` | Keep different scaffolds, registers, targets, or user-defined families separate. |
| `TOKEN_MODE` | `letters` for standard amino acids, `braces` for `{A}{NMe-Phe}`, or `space` for `A NMe-Phe`. |
| `CANDIDATE_MODE` | `auto`, `exact`, or `heuristic`; see below. |
| `OUTPUT_DIR` | Output folder. Re-running replaces files in that folder. |
| `RUN_LABEL` | Descriptive title for figures and run metadata. |
| `INVALID_SEQUENCE_POLICY` | `skip` invalid sequence records with a report, or `error` to stop. |
| `LOGO_MODE` | `probability` for residue frequencies (default), or `information` for bits. |
| `LOGO_POSITIONS_PER_ROW` | Wrap long sequences across rows in each library logo (default: 30). |

Different lengths are always separated. Equal length does **not** establish a common binding register: prepare positionally comparable inputs or supply a group column. No sequence alignment or hierarchical clustering is performed. Duplicate sequence/group entries are collapsed, their identifiers retained, and the maximum weight used. Negative or missing weights are rejected rather than transformed silently.

The library budget is **global across all groups and lengths**, not per group.

## How it works

1. Pool residues at each position within candidate groups of starting sequences.
2. Define each library as all combinations of these positional options. For example,
   `[A,V] × [L] × [G,A] × [S]` produces four peptides.
3. Use binary integer linear programming to maximize weighted unique-sequence coverage
   within the library-count limit. Candidates exceeding the size limit are excluded.
5. Enumerate all selected library members and verify sizes and coverage independently.

### Search scope matters

Enumerating all subsets of 50 sequences would require considering up to 1,125,899,906,842,623 nonempty subsets. The implementation instead deduplicates feasible position-set combinations, prunes combinations exceeding the size cap, and bounds the work performed.

- **`auto`** attempts exhaustive construction within the configured budgets. If a budget
  is reached, it retains a bounded pool and adds multi-start greedy candidates.
- **`exact`** requires exhaustive construction and stops with an error if a budget is reached.
- **`heuristic`** constructs a bounded pool directly using multi-start greedy expansion.

The integer solver optimizes selection over the generated pool. A global coverage claim is made only when candidate construction is exhaustive and the primary solve is certified optimal. Otherwise, results are explicitly labeled as pool-optimal or feasible but uncertified. No hierarchical clustering is required.

## Outputs

| File | Contents |
| --- | --- |
| `design.csv` | Selected libraries, sizes, covered parent IDs, and positional options. |
| `composition.csv` | One row per position, with allowed monomers and fixed/variable status. |
| `members.csv` | Every peptide in each selected library; labels input versus new combinations. |
| `coverage.csv` | Every unique starting sequence, its weight, and the libraries covering it. |
| `input_records.csv` | All original input records and additional columns, including excluded records. |
| `normalized_sequences.csv` | Accepted normalized records before duplicate sequences are merged. |
| `excluded_sequences.csv` | Excluded sequences, original record numbers, identifiers, and reasons. |
| `run_metadata.json` | Input hash, settings, versions, search completeness, and solver status. |
| `coverage.png/.svg` | Input coverage map and library sizes. |
| `composition.png/.svg` | Number of residue options at each position. |
| `logos/LIBxx_logo.png/.svg/.pdf` | Sequence logo for each selected library, also displayed in the notebook. |
| `logos/all_library_logos.pdf` | One PDF containing every library logo, one library per page. |
| `logos/logo_frequencies.csv` | Exact residue counts and frequencies for each library and position. |
| `logos/logo_heights.csv` | Plotted letter heights and the selected scale. |

JSON-valued CSV fields preserve monomer names containing commas. Synthesizing the same peptide in two libraries counts twice in the summed-library-size objective; coverage
of a starting sequence is counted only once.

### Sequence logos

Every designed library receives a WebLogo-style sequence logo. Letter height is the fraction of **all designed members** carrying that residue at that position.
Each library member is counted once, irrespective of parent scores. For a complete combinatorial library, allowed residues at a position therefore have equal frequency.
A fixed position shows one full-height letter. These are design frequencies, not measured molar abundances or frequencies among the input hits.

Optional `LOGO_MODE = "information"` scales letters by information relative to a uniform 20-amino-acid background. A column's total height is log₂(20) minus its
Shannon entropy; each letter receives its frequency times that height. No finite-sample correction is applied to the fully enumerated design. Information mode supports
canonical letters; probability mode also renders explicit custom monomer labels. The plotting implementation uses the existing Matplotlib dependency and needs no WebLogo installation or online service.

## Example and validation

`example` contains input and outputs from the example data file. The example data file was generated through a PepMLM-650M run[1] targeting MDM2. The run used a PDB-derived target sequence to generate 50 candidate peptides of length 10 with a top-K setting of 3. For the generated peptide sequences `X` denotes an unspecified residue, so sequences containing `X` are excluded from library design and reported separately. The column Pseudo Perplexity is a model score for the generated peptide. Lower values indicate a better score. The model score is not considered for library design in the example output. It can be used if explicitly configured as a weight.

[1] Chen, T., Dumas, M., Watson, R., et al. “PepMLM: Target Sequence-Conditioned Generation of Therapeutic Peptide Binders via Span Masked Language Modeling.” *arXiv* (2024). https://doi.org/10.48550/arXiv.2310.03842
