# FlaGs2 — User Guide

Predicting protein functional association by analysis of conservation of genomic context (Flanking Genes).

---

## Installing

With conda:

```bash
bash build.sh
conda activate FlaGs2
```

`build.sh` creates the environment, verifies it, and offers the two large
optional extras (sismis and the Pfam-A database). It is safe to re-run — it
detects what is already installed and skips it.

Or install the core dependencies by hand:

```bash
pip install biopython "pyhmmer>=0.12,<0.13" requests
```

Optional features need extra pieces, each only required if you use the matching
flag:

| Flag | Needs |
|---|---|
| `--tree`, `--tree_order` | `mafft` and `VeryFastTree` on `PATH` |
| `--domains` | an HMM database — run `pfamA_loader.sh` to fetch Pfam-A |
| `--tmhmm`, `--signalp` | `pip install pybiolib` and a network connection |
| `--sismis` | `pip install sismis` |

If a tool is missing, FlaGs2 prints a warning, skips that feature, and finishes
the rest of the run.

Keep `pyhmmer` in the 0.12 series. FlaGs2 is developed against it, and sismis
(via gecco) requires it — an unpinned install can leave the two in conflict.

---

## The input list

One query per line. Two accepted forms, which can be mixed in the same file:

```
WP_047256880.1                          # protein only — FlaGs2 finds the genome
WP_047256880.1    GCF_000001765.3       # protein + genome, tab-separated
MGYG000454827_00001   MGYG000454827     # MGnify genome
```

**Protein only.** FlaGs2 resolves the genome through NCBI. This works for RefSeq
and GenBank proteins; `-m` controls how many genomes a protein may expand to when
it appears in several.

**Protein + genome.** Skips the lookup. The genome may be:

- an NCBI assembly — `GCF_...` (RefSeq) or `GCA_...` (GenBank)
- an MGnify Genomes accession — `MGYG...`

For MGnify genomes the protein accession must be the **exact locus tag** from
that genome's annotation (`MGYG000454827_00001`, not an `MGYP...` protein ID).
MGnify has no equivalent of NCBI's lookup service, so a bare MGnify protein
accession cannot be resolved to its genome — always supply the pair.

To find real locus tags for a genome:

```bash
curl -s "https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000454827/downloads/MGYG000454827.faa" \
  | grep '^>' | head | sed 's/^>//' | cut -d' ' -f1
```

---

## Running it

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O results
```

`-u` is required by NCBI on any Entrez request. `-O` names the output directory
*and* becomes the prefix on every file inside it.

That is the whole default pipeline: it resolves each query to a genome, pulls the
flanking genes, clusters them, and writes the figure and tables. Everything else
is optional. See [All options](#all-options) for the full list and
[Examples](#examples) for common combinations.

## All options

### Required

| Option | Description |
|---|---|
| `-i`, `--input_list FILE` | Query list, one per line. See [The input list](#the-input-list). |
| `-u`, `--user_email ADDR` | Your email. NCBI requires it on every Entrez request. Not used for anything else. |

### Where genomes come from

| Option | Default | Description |
|---|---|---|
| `--use_local DIR` | — | Search a directory of local genomes before going to NCBI. A genome is a `.gff` and `.faa` sharing a basename; `.fna` and RNA FASTAs are picked up if present. Files may be gzipped. Anything not found falls back to NCBI. |
| `-m`, `--max_assemblies N` | `1` | How many genomes one protein may expand to when it occurs in several. Each genome becomes its own row. Raise to compare strains. Above `1`, row labels in the figures become `protein\|genome` so the rows stay distinguishable. |
| `-api`, `--api_key KEY` | — | NCBI API key. Also raises the download rate cap from 5/s to 10/s. |
| `-tmp`, `--temporary DIR` | `./genomes` | Where downloads are stored. Deleted at the end unless `-k`. |
| `-k`, `--keep` | off | Keep downloaded genomes instead of deleting them. Useful for reruns — the directory can be fed straight back in via `--use_local`. |

### What gets analysed

| Option | Default | Description |
|---|---|---|
| `-g`, `--gene N` | `4` | Flanking genes to take each side of the query. |
| `-e`, `--ethreshold X` | `1e-3` | Inclusion E-value for clustering. Lower is stricter, giving more and smaller families. |
| `-n`, `--number N` | `3` | Jackhmmer iterations. More iterations find remoter homology but blur family boundaries. |
| `--cluster_rna` | off | Also cluster flanking RNA genes into families. Uses nhmmer on RNA sequences where available, otherwise groups by product name. |

### Output and progress

| Option | Default | Description |
|---|---|---|
| `-O`, `--output DIR` | `output` | Result directory. A `_YYYYMMDD_HHMMSS` stamp of the run start is appended so repeated runs do not overwrite each other, and the stamped name is also the prefix on every file inside, so `-O myrun` produces `myrun_20260810_093134/myrun_20260810_093134_neighbors.svg`. |
| `--no_timestamp` | off | Use `-O` verbatim, without the stamp. Repeated runs then overwrite each other; use it when a pipeline needs a fixed path. |
| `-vb`, `--verbose` | off | Per-stage progress and a timing breakdown. Worth using on any long run. |
| `-v`, `--version` | — | Print the version and exit. |
| `-h`, `--help` | — | Print all options and exit. |

### Optional figures and annotation

Each needs an extra dependency. If it is missing, FlaGs2 warns, skips that
feature, and completes the rest of the run.

| Option | Needs | Description |
|---|---|---|
| `--tree` | mafft, VeryFastTree | Also build a phylogenetic tree with the neighbourhoods aligned to its leaves (`_tree.svg`, `_tree.nwk`). Does not change the main figure. |
| `--tree_order` | mafft, VeryFastTree | Order rows by tree leaf order, in the main figure and in `_operon.tsv`. Implies `--tree`. |
| `--domains` | `--hmmdb` | Scan flanking proteins for domains and write `_domains.svg` plus `_domains.tsv`. |
| `--hmmdb FILE` | — | HMM database for `--domains`, e.g. `pfam_db/Pfam-A.hmm`. Run `pfamA_loader.sh` to fetch it. |
| `--clans FILE` | — | `Pfam-A.clans.tsv.gz`. Colours domains by clan rather than family, which groups related domains together. |
| `--tmhmm` | pybiolib, network | Predict transmembrane regions with DeepTMHMM, drawn as double red dotted lines on the domain figure. Uploads your sequences to the BioLib cloud. |
| `--signalp` | pybiolib, network | Predict signal peptides with SignalP-6, drawn as black triangles on the domain figure. Also uploads sequences. |
| `--sismis` | sismis | Scan each genome for secretion systems and write `_secretion.tsv` plus `_secretion.svg`, noting which neighbourhoods each hit overlaps. Downloads the genomic FASTA per genome. |

`--tmhmm`, `--signalp` and `--sismis` run concurrently with each other, since
each spends a lot of time waiting on a remote service or a subprocess.

### Tree figure appearance

| Option | Default | Description |
|---|---|---|
| `-ts`, `--tshape N` | `20` | Size of the gene triangles in the tree figure. |
| `-tf`, `--tfontsize N` | `13` | Font size inside those triangles. |

### Performance

| Option | Default | Description |
|---|---|---|
| `-c`, `--cpu N` | auto | Worker cap for clustering, domain scanning and downloads. |

Download rate is capped independently of `-c` at 5 requests/second, or 10 with
`-api`. This is deliberate: raising the worker count on a fast machine would
otherwise raise the request rate and get the run throttled or rejected.

### Legacy arguments

Original FlaGs argument names still work and are translated automatically:

| Legacy | Becomes |
|---|---|
| `-p`, `--proteinList` | `-i` |
| `-a`, `--assemblyList` | `-i` |
| `-l`, `--localGenomeList` | `-i` |
| `-ld`, `--localGenomeDirectory` | `--use_local` |
| `-o`, `--out_prefix` | `-O` |
| `-db` | `--hmmdb` |
| `-t` | `--tree` |
| `-to` | `--tree_order` |
| `-r`, `--redundant N` | `-m N` (`-r a` becomes `-m 100000`) |

Passing several legacy input files merges them into `_flags2_merged_input.txt`
in the working directory, reordering two-column files to protein-first.

---

## Examples

Default run — neighbours figure and data tables:

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O myrun
```

Wider neighbourhood across several strains, ordered by a tree:

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O myrun \
  -g 6 -m 5 --tree --tree_order -vb
```

Domain annotation with clan colouring:

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O myrun \
  --domains --hmmdb pfam_db/Pfam-A.hmm --clans pfam_db/Pfam-A.clans.tsv.gz
```

Everything on, with an API key for faster downloads:

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O myrun -api YOUR_KEY \
  --tree --tree_order --domains --hmmdb pfam_db/Pfam-A.hmm \
  --tmhmm --signalp --sismis --cluster_rna -vb
```

Reuse genomes from a previous run instead of downloading again:

```bash
python3 FlaGs2.py -i input.txt -u you@example.com -O run1 -k
python3 FlaGs2.py -i input.txt -u you@example.com -O run2 --use_local ./genomes
```

---

## Output files

Every file is prefixed with the output directory's name, stamp included — a run
with `-O results` writes `results_20260810_093134/results_20260810_093134_operon.tsv`.
The tables below drop the stamp and write `results_...` for readability.

**Figures**

| File | Contents |
|---|---|
| `results_neighbors.svg` | the main diagram: one row per query, genes as arrows coloured by family |
| `results_tree.svg` | same rows aligned to a phylogenetic tree (`--tree`) |
| `results_tree.nwk` | the tree in Newick format |
| `results_domains.svg` | neighbourhoods with domains, TM regions and signal peptides drawn on |
| `results_secretion.svg` | neighbourhoods with predicted secretion systems marked (`--sismis`) |

**Tables**

| File | Contents |
|---|---|
| `results_operon.tsv` | one row per flanking gene: query, genome, family, strand, offset, coordinates, length, contig, product |
| `results_clusters.tsv` | each family and its members |
| `results_outdesc.txt` | families as readable blocks: `family(occurrences)`, accession, product description |
| `results_domains.tsv` | one row per domain hit: protein, family, domain, clan, coordinates, E-value (`--domains`) |
| `results_jackhits.tsv` | per-protein jackhmmer inclusion lists — the audit trail behind the families |
| `results_speciesInfo.txt` | genome and organism per query row |
| `results_QueryStatus.txt` | which genomes each query resolved to, and whether it produced a row |
| `results_accessionIssues.txt` | queries that produced nothing, and why |
| `results_flankgene_Report.log` | each neighbourhood as a compact family chain |
| `results_secretion.tsv` | Sismis hits and which neighbourhoods they overlap (`--sismis`) |
| `results_sismis_diagnostics.txt` | per-genome Sismis status: scanned, nothing found, or skipped and why (`--sismis`) |

### Reading the figure

Genes are arrows pointing in their direction of transcription, normalised so the
query always points right. Arrows sharing a colour and number are one family.
The query itself is outlined in black at the centre of each row.

Grey means the gene had no family (a singleton). RNA genes keep a green outline,
pseudogenes navy.

### Reading `results_operon.tsv`

The `query` column is the query protein and `assembly` is the genome it was found
in. They are separate columns because one protein appearing in several genomes
produces one row per genome, so the pair — not the protein alone — identifies a
row. The `accession` column is the flanking gene itself. `offset` is the position
relative to the query: `0` is the query, negative upstream, positive downstream.
`contig` is the sequence the gene lies on, which is what Sismis hits are matched
against.

With `--tree_order`, rows appear in tree leaf order rather than input order.

---

## When something goes wrong

**"No flanking neighbourhoods could be extracted for any query."**
Nothing matched. Check `results_accessionIssues.txt` — it separates *no genome
resolved* from *genome found but the protein was not in it*. The second usually
means the accession does not appear in that genome's annotation, which for
MGnify genomes normally means the locus tag is wrong.

**Queries silently missing from the figure.**
`results_QueryStatus.txt` lists every query and whether it produced a row.

**Downloads look stuck.**
Run with `-vb`; downloads report as they complete. Requests are rate-limited
(5/s, or 10/s with `-api`) so a large list takes a while by design — this keeps
NCBI and EBI from rejecting the run.

**The diagram is unexpectedly wide.**
A neighbour lying very far from the query stretches the canvas. This normally
means a fragmented assembly where the query sits near a contig edge.

**A tool was skipped.**
Warnings name the missing dependency and the install command. The run continues
without that feature.