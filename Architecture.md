# FlaGs2 — Architecture Overview

A developer-facing map of how the pipeline fits together. For usage, see
`USER_GUIDE.md`.

---

## Shape of the program

A linear pipeline. Each stage takes the previous stage's output and adds to it;
nothing loops back.

```
input list
    ↓  AccessionListReader
protein accessions (+ optional genome)
    ↓  LocalGenomeResolver → ProteinAssemblyMapper
protein → [genome ids]
    ↓  AssemblyDownloader / MgnifyGenomeDownloader
genome id → GenomeFiles(gff, faa, rna, genome)
    ↓  NeighborhoodExtractor
[FlankingGene] + protein/RNA sequence tables
    ↓  NeighborhoodClusterer / RnaClusterer
families: [[accession, ...], ...]
    ↓  OperonView / NeighborhoodVisualizer / ReportWriter
SVG figures + TSV tables + FASTA sequences
```

Optional stages (`--tmhmm`, `--signalp`, `--sismis`) hang off the side after
extraction and feed extra layers into the figures.

---

## Modules

| Module | Responsibility | Imported |
|---|---|---|
| `FlaGs2.py` | data pipeline and CLI | always |
| `flags2_view.py` | shared styling, `OperonView` renderer | always |
| `flags2_tree.py` | MAFFT + trimming + VeryFastTree or IQ-TREE, tree figure | on `--tree`/`--iqtree` |
| `flags2_domains.py` | pyhmmer domain scan | on `--domains` |
| `flags2_features.py` | DeepTMHMM / SignalP via BioLib | on `--tmhmm`/`--signalp` |
| `flags2_secretion.py` | Sismis secretion-system scan | on `--sismis` |

Optional modules are imported lazily inside the branch that needs them, so a
default run never touches `mafft`, `pybiolib` or `sismis` — and a missing
dependency degrades to a warning rather than an import error at startup.

---

## Core types

Three small carriers move data between stages:

```python
GenomeFiles(gff, faa, rna, genome)   # paths; any may be None
FlankingGene(accession, strand, start, end, product, offset, query, is_rna, contig)
families: List[List[str]]            # each inner list is one family's accessions
```

`GenomeFiles` is what every genome source returns, whether downloaded from NCBI,
downloaded from MGnify, or found on disk. Downstream code never learns where a
genome came from.

`FlankingGene.offset` is signed distance from the query (`0` = query), already
flipped when the query is on the minus strand, so renderers can lay rows out
without knowing about strands. `FlankingGene.query` holds a **row id**
(`protein|genome`), not a bare protein — one protein in several genomes makes
several independent rows.

---

## Resolution and download

Genomes are resolved in priority order: local directory, then NCBI. Anything
matching `MGYG\d+` is routed to MGnify instead.

**`ProteinAssemblyMapper`** turns bare protein accessions into genome ids via
NCBI's IPG database in one batched request. `XP_` proteins are not in IPG, so
they take a separate BioProject → assembly path. Under `--no_cross_db` it drops
assemblies from the other database (RefSeq protein → `GCF_`, INSDC protein →
`GCA_`) before truncating to `-m`. Filtering happens after the fetch rather than
in the query, which leaves the IPG round trip unchanged and lets
`dropped_cross_db` record what was excluded for the `-vb` summary. Assemblies
supplied in the input file bypass the mapper entirely and are never filtered.

**`_GenomeDownloader`** holds everything the two remote sources share: HTTP
session, retry policy, worker pool, rate limiter, per-file streaming. Subclasses
supply only their URL scheme:

- `AssemblyDownloader` — NCBI's partitioned FTP layout, so it must first list a
  directory to discover the versioned name, then fetch gzipped files from it.
- `MgnifyGenomeDownloader` — one API call returns direct URLs; files are plain
  text. MGnify has no RNA-only FASTA, so the RNA slot stays empty and the
  extractor cuts RNA sequences out of the genome FASTA instead.

Each assembly's files are fetched concurrently, and assemblies are fetched
concurrently with each other.

**`RateLimiter`** caps requests per second across all threads by handing out
timed slots. This is deliberately separate from the worker count: more threads
or a faster machine would otherwise mean a faster request rate and rejected
downloads. The cap is 5/s, or 10/s when an NCBI API key is supplied.

---

## Extraction

`NeighborhoodExtractor` parses GFF and FASTA into gene records, then takes a
window of `±flank` genes around the query.

The GFF parser handles two annotation styles without being told which it has:

- **NCBI PGAP** — a `gene` feature followed by a `CDS` carrying the accession in
  `protein_id=`. The CDS attaches itself to the preceding gene record.
- **Prokka / Prodigal** (MGnify) — bare `CDS` features with no parent gene, where
  the locus tag *is* the accession. With no gene record to attach to, the CDS
  becomes its own record.

`_cds_accession` resolves the accession by trying `protein_id=`, then `ID=`
(stripping NCBI's `cds-` prefix), then the locus tag, and only then `Name=`.
`Name=` is last because Prodigal often sets it to a gene symbol such as `hisZ_2`,
which matches no FASTA record.

Genes are sorted by `(contig, start)` after parsing. Neighbours are selected by
list index, so file order must equal genomic order — true for NCBI, not for
MGnify, which writes all CDS features before all ncRNA features.

Parsed GFF and FASTA tables are cached per assembly, since one genome is usually
queried by several proteins. Genome FASTAs are the exception: they are loaded
only when an RNA actually needs slicing, and only one is held at a time, because
they are orders of magnitude larger than the other tables.

---

## Clustering

`NeighborhoodClusterer` runs jackhmmer with every flanking protein as a query
against all of them, building an adjacency map from the included hits, then takes
connected components as families. Clustering is symmetric by construction: if A
hits B, they end up in the same component regardless of direction. The adjacency
map is kept on the clusterer after the run and handed to `ReportWriter`, which
writes it as `_jackhits.tsv` — the only record of *why* two proteins share a
family.

`_outdesc.txt` orders its family blocks by total occurrences descending, which is
how main's file read. That ordering is local to the file: family *labels* come
from `family_numbers()` and stay consistent with the figures and `_clusters.tsv`,
so the labels in `_outdesc.txt` are not ascending. Main had both because it
derived its numbering from the same sort.

`RnaClusterer` mirrors this with nhmmer for RNA genes, falling back to grouping
by normalised product name when no nucleotide sequence is available.

---

## Tree building

`TreeBuilder` runs three stages: MAFFT `--auto`, a gap-threshold column trim, then
an inference engine. The trim keeps columns where at least `gap_threshold` of the
sequences carry a residue — the same rule as `trimal -gt`, which is what ete3's
`trimal01` ran in the old pipeline. It is implemented in `_trim` rather than
shelled out because `-gt` is a deterministic column filter and adding a binary
dependency for ten lines is not worth it; the heuristic modes (`-automated1`,
`-gappyout`) would need the real trimal. If trimming would empty the alignment,
the untrimmed one is used instead.

The engine is `veryfasttree` by default and `iqtree` under `--iqtree`
(ModelFinder, plus 1000 ultrafast bootstrap replicates when there are at least
four taxa — below that IQ-TREE refuses to bootstrap). The binary is looked up as
`iqtree3`, `iqtree2`, then `iqtree`, since distributions disagree on the name.

The trimmed alignment stays on the builder as `.alignment` so `main()` can write
it as `_tree.aln` without re-running anything.

`write_run_info()` runs immediately after the output path is resolved, before any
network or file work, so a crashed run still leaves `_runinfo.txt` and a copy of
the input list. It reads defaults off the parser rather than hardcoding them, so
the "set explicitly" split stays correct as options are added. Anything named in
`SECRET_ARGS`/`SECRET_FLAGS` is masked in both the option dump and the recorded
command line — `_redact_command_line` handles `--flag value` and `--flag=value`,
since only covering the first form leaks the key in the second.

---

## Rendering and reporting

`OperonView` draws all three figures — neighbourhoods, domains, secretion — from
one code path with a `mode` switch, so layout stays identical between them.
`flags2_tree.NeighborhoodVisualizer` is separate because it must align rows to
tree leaves.

Text width is measured from a per-character Arial metrics table rather than
estimated, and the font stack is pinned to metrically identical faces (Arial,
Liberation Sans, Helvetica). A generic `sans-serif` would resolve differently per
platform and render text at a width the layout never reserved.

`ReportWriter` owns every TSV, text and FASTA output. It takes the leaf order under
`--tree_order` and re-keys its per-row map, so the tables and the figures agree
on row order rather than only the figures being sorted. The FASTA outputs come
from the extractor's sequence tables; queries are excluded from
`_flankgene.fasta` because `extractor.sequences` holds every gene in the window,
query included.

---

## Concurrency

| Where | Model | Bound by |
|---|---|---|
| Downloads | thread pool, ≤10 assemblies, files within each in parallel | network, rate limiter |
| Clustering | thread pool over queries, `cpus=1` each | CPU |
| Domain scan | pyhmmer internal threads | CPU |
| tmhmm / signalp / sismis | 3-thread pool, run together | cloud / subprocess |

The last group runs concurrently because each spends its time waiting on a remote
service or a subprocess rather than on local CPU. Each measures its own elapsed
time, so the `-vb` timing table reports true per-tool cost even though the three
overlap in wall time. `TOTAL` is measured wall clock, not the sum of stages —
the stages overlap, so summing them would over-count.

### BioLib submission constraints

Three things in `flags2_features.py` look removable and are not:

- **`_SUBMIT_LOCK` around `load()` + `cli()`.** pybiolib's `attempt_sign_in()` is
  a check-then-act on a shared singleton, and its `UserState` lock file is
  created with `fail_fast_on_lock_acquire=True`, so a second concurrent caller
  gets no retry. Worse, the loser sets `_is_in_memory_only` on the shared object,
  which makes the winner's `__exit__` skip releasing the lock file — leaving a
  stale lock that degrades *every later run on that machine*. `warm_up()` burns
  the one-shot sign-in on the main thread before the pool starts; the lock covers
  the rest. `job.wait()` stays outside it so the cloud jobs still overlap.
- **The `chdir` into the scratch directory.** The fasta argument must reach
  `app.cli()` as a bare relative name. Relative names are mounted at
  `/query.fasta` and passed through unchanged; an absolute path is mounted under
  a hashed directory and rewritten *without* its leading slash, so the app
  receives `hash/query.fasta` and never finds its input. `chdir` is process-wide,
  which is safe only because it is held inside `_SUBMIT_LOCK` and because
  `main()` absolutises `args.output` and `args.temporary` before any thread
  starts — the concurrent sismis task holds paths under both.
- **SignalP's `--output_dir output`.** The app's own `generate_output.py` reads
  `output/output.json`. Any other directory name makes the remote step fail.

---

## Extending it

**Adding a genome source.** Subclass `_GenomeDownloader`, implement `_fetch_one`
returning `GenomeFiles`, define `SUFFIX`, and add a routing test alongside
`is_mgnify_accession` in `main()`. Nothing downstream needs to change. If the
source uses an unusual GFF dialect, extend `_cds_accession` rather than branching
on source anywhere else.

**Adding a per-protein annotation layer.** Follow `flags2_features.py`: return
`{accession: [(kind, start, end)]}`, add a lazy import and a `_run_*` function in
`main()`, and register it in the parallel task group. If it goes through BioLib,
subclass `_BioLibScanner` rather than calling `biolib` directly, so it inherits
`_SUBMIT_LOCK` and the scratch-directory handling described above.