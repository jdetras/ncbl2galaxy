# Architecture

## Data flow
1. Query PubMed with a rice-focused search term.
2. Link PubMed IDs to SRA IDs through NCBI E-utilities.
3. Resolve run accessions and ENA FASTQ metadata.
4. Group runs by sample accession.
5. Upload reads and reference FASTA to Galaxy.
6. Invoke workflow per run with sample-aware history naming and resumable state.

## Reliability
- HTTP retries with exponential backoff for NCBI, ENA, and Galaxy API calls.
- Resume support via `.ncbi_to_galaxy_state.json` to avoid reprocessing successful runs.

## Security model
- No secrets in source; Galaxy API key provided at runtime.
- CI security checks: Bandit, dependency audit, and secret scanning.
- Repository security policy in `SECURITY.md`.

## Galaxy compatibility matrix
Validated workflow templates target these Galaxy tool wrappers:
- `bwa_mem2` `2.2.1+galaxy0`
- `samtools_sort` `2.0.5`
- `freebayes` `1.3.6+galaxy0`

If your Galaxy server has different wrapper versions, import workflow templates and remap missing tools in the workflow editor.
