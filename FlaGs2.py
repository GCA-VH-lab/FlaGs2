import argparse
import gzip
import os
import re
import shutil
import sys
import tempfile
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, NamedTuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock

from flags2_view import OperonView, family_numbers

def plural(n, word, plural_form=None):
	return "{} {}".format(n, word if n == 1 else (plural_form or word + "s"))


def translate_legacy_args(argv):
	rename = {
		"-a": "-i", "--assemblyList": "-i",
		"-p": "-i", "--proteinList": "-i",
		"-l": "-i", "--localGenomeList": "-i",
		"-ld": "--use_local", "--localGenomeDirectory": "--use_local",
		"-o": "-O", "--out_prefix": "-O",
		"-db": "--hmmdb",
	}
	flag_rename = {
		"-t": "--tree", "-to": "--tree_order",
	}
	out = []
	input_files = []
	paired_files = []
	i = 0
	while i < len(argv):
		tok = argv[i]
		#
		if tok.startswith("-") and "=" in tok:
			name, _, val = tok.partition("=")
			argv = argv[:i] + [name, val] + argv[i + 1:]
			tok = argv[i]

		if tok in ("-r", "--redundant"):
			val = argv[i + 1] if i + 1 < len(argv) else ""
			new_val = "100000" if val.lower() == "a" else val
			out += ["-m", new_val]
			i += 2
			continue
		if tok in ("-p", "--proteinList"):
			if i + 1 < len(argv):
				input_files.append(argv[i + 1])
			i += 2
			continue
		if tok in ("-a", "--assemblyList", "-l", "--localGenomeList"):
			if i + 1 < len(argv):
				paired_files.append(argv[i + 1])
			i += 2
			continue
		if tok in rename:
			out += [rename[tok], argv[i + 1]] if i + 1 < len(argv) else [rename[tok]]
			i += 2 if i + 1 < len(argv) else 1
			continue
		if tok in flag_rename:
			out.append(flag_rename[tok])
			i += 1
			continue
		out.append(tok)
		i += 1
	legacy_inputs = input_files + paired_files
	if legacy_inputs:
		if len(input_files) == 1 and not paired_files:
			out += ["-i", input_files[0]]      
		else:
			merged = "_flags2_merged_input.txt"
			with open(merged, "w") as mf:
				for path in input_files:        
					with open(path) as f:
						for line in f:
							if line.strip():
								mf.write(line.strip() + "\n")
				for path in paired_files:        
					with open(path) as f:
						for line in f:
							line = line.strip()
							if not line:
								continue
							cols = line.replace(" ", "").split("\t")
							if len(cols) >= 2:
								mf.write("{}\t{}\n".format(cols[1], cols[0]))
							else:
								mf.write(line + "\n")
			if len(legacy_inputs) > 1:
				print("Note: combined {} legacy input file(s) into {}".format(
					len(legacy_inputs), merged))
			out += ["-i", merged]
	return out
class AccessionListReader: 
	def __init__(self, path):
		self.path = path

	def read(self):
		proteins_assembly = []
		proteins_only = []
		with open(self.path, 'r') as f:
			for n, line in enumerate(f, 1):
				line = line.strip()
				if not line:
					continue
				if '\t' in line:
					fields = line.split('\t')
					if len(fields) >= 2 and fields[0] and fields[1]:
						proteins_assembly.append([fields[0], fields[1]])
					else:
						print("Warning: line {} of the input list is malformed, skipping it: {!r}".format(n, line))
				else:
					proteins_only.append(line)
		return proteins_assembly, proteins_only


class ProteinAssemblyMapper: 
	def __init__(self, email: str, max_assemblies: int = 1,
				 api_key: Optional[str] = None, ncbi_time: float = 0.4):
		self.max_assemblies = max_assemblies
		self.ncbi_time = ncbi_time
		self.accessions_in: Dict[str, Dict[str, set]] = {}
		Entrez.email = email
		if api_key:
			Entrez.api_key = api_key

	def map(self, proteins: List[str]) -> Dict[str, List[str]]:
		if not proteins:
			return {}
		xp = [p for p in proteins if p.startswith("XP_")]
		ipg = [p for p in proteins if not p.startswith("XP_")]
		found: Dict[str, set] = {}
		if ipg:
			found.update(self._map_ipg(ipg))
		for acc in xp:
			found[acc] = self._map_xp(acc)
		gcf_first = lambda a: (0 if a[:3] == "GCF" else 1, a)
		return {acc: sorted(asm, key=gcf_first)[:self.max_assemblies]
				for acc, asm in found.items()}

	def _map_ipg(self, proteins: List[str]) -> Dict[str, set]:
		queries = set(proteins)
		time.sleep(self.ncbi_time)
		handle = Entrez.efetch(db="ipg", id=",".join(proteins),
							   rettype="ipg", retmode="text")
		data = handle.read()
		handle.close()
		if isinstance(data, bytes):
			data = data.decode("utf-8", errors="replace")
		groups: Dict[str, dict] = {}
		for line in data.splitlines():
			if line[0:2] == "Id" or not re.search(r"GC._\d*\.\d", line):
				continue
			fields = line.rstrip().split("\t")
			ipg_id, acc, assembly = fields[0], fields[6], fields[-1]
			grp = groups.setdefault(ipg_id, {"queries": set(), "by_asm": {}})
			grp["by_asm"].setdefault(assembly, set()).add(acc)
			if acc in queries:
				grp["queries"].add(acc)

		found = {acc: set() for acc in proteins}
		for grp in groups.values():
			for q in grp["queries"]:
				self.accessions_in.setdefault(q, {})
				for assembly, accs in grp["by_asm"].items():
					found[q].add(assembly)
					self.accessions_in[q].setdefault(assembly, set()).update(accs)
		return found

	def _map_xp(self, accession: str) -> set:
		try:
			time.sleep(self.ncbi_time)
			handle = Entrez.efetch(db="protein", id=accession,
								   rettype="gbwithparts", retmode="text")
			record = SeqIO.read(handle, "genbank")
			handle.close()
		except Exception:
			return set()
		bioprojects = [x.split(":", 1)[1] for x in record.dbxrefs
					   if x.split(":", 1)[0] == "BioProject"]
		assemblies: set = set()
		for bp in bioprojects:
			assemblies.update(self._assemblies_for_bioproject(bp))
		self.accessions_in.setdefault(accession, {})
		for asm in assemblies:
			self.accessions_in[accession].setdefault(asm, set()).add(accession)
		return assemblies

	def _assemblies_for_bioproject(self, bioproject: str) -> set:
		try:
			time.sleep(self.ncbi_time)
			search = Entrez.read(Entrez.esearch(db="bioproject", term=bioproject))
			ids = search.get("IdList", [])
			if not ids:
				return set()
			time.sleep(self.ncbi_time)
			links = Entrez.read(Entrez.elink(dbfrom="bioproject", db="assembly",
											 id=",".join(ids)))
			asm_ids = []
			for linkset in links:
				for db in linkset.get("LinkSetDb", []):
					asm_ids.extend(link["Id"] for link in db.get("Link", []))
			if not asm_ids:
				return set()
			time.sleep(self.ncbi_time)
			summary = Entrez.read(Entrez.esummary(db="assembly",
												  id=",".join(asm_ids)))
			docs = summary["DocumentSummarySet"]["DocumentSummary"]
			out = set()
			for d in docs:
				acc = d.get("AssemblyAccession", "")
				if re.match(r"GC._\d+\.\d", acc):
					out.add(acc)
			return out
		except Exception:
			return set()


class RateLimiter: 

	def __init__(self, rate: float):
		self.min_interval = 1.0 / rate if rate > 0 else 0.0
		self._lock = threading.Lock()
		self._next_slot = time.monotonic()

	def wait(self):
		if self.min_interval <= 0:
			return
		with self._lock:
			now = time.monotonic()
			start = max(now, self._next_slot)
			self._next_slot = start + self.min_interval
		delay = start - now
		if delay > 0:
			time.sleep(delay)


class GenomeFiles(NamedTuple):
	gff: Optional[str] = None
	faa: Optional[str] = None
	rna: Optional[str] = None       
	genome: Optional[str] = None    


class _GenomeDownloader: 
	SUFFIX: Dict[str, str] = {}

	def __init__(self, out_dir: Optional[str] = None, workers: int = 8,
				 rate: float = 5.0, want_rna: bool = False, want_genome: bool = False):
		self.out_dir = out_dir or tempfile.gettempdir()
		self.workers = workers
		self.want_rna = want_rna
		self.want_genome = want_genome
		self.limiter = RateLimiter(rate)
		self.failures: Dict[str, str] = {}  
		os.makedirs(self.out_dir, exist_ok=True)
		self.session = requests.Session()
		retry = Retry(total=5, backoff_factor=0.5, respect_retry_after_header=True,
					  status_forcelist=[429, 500, 502, 503, 504],
					  allowed_methods=frozenset(["GET"]))
		size = workers * 4
		self.session.mount("https://", HTTPAdapter(
			max_retries=retry, pool_connections=size, pool_maxsize=size))

	def download_many(self, assemblies: List[str], progress_cb=None) -> Dict[str, GenomeFiles]:
		if not assemblies:
			return {}
		results = {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			futures = {pool.submit(self._fetch_one, a): a for a in assemblies}
			for done, fut in enumerate(as_completed(futures), 1):
				results[futures[fut]] = fut.result()
				if progress_cb:
					progress_cb(done, len(assemblies))
		return results

	def _slots(self) -> List[str]:
		slots = ["gff", "faa"]
		if self.want_rna and "rna" in self.SUFFIX:
			slots.append("rna")
		if self.want_genome:
			slots.append("genome")
		return slots

	def _fetch_files(self, jobs: Dict[str, Tuple[Optional[str], str]]) -> GenomeFiles:

		with ThreadPoolExecutor(max_workers=max(len(jobs), 1)) as pool:
			done = {slot: pool.submit(self._stream, url, local)
					for slot, (url, local) in jobs.items()}
			return GenomeFiles(**{slot: jobs[slot][1]
								  for slot, fut in done.items() if fut.result()})

	def _stream(self, url: Optional[str], local: str) -> bool:
		if not url:
			return False
		for _ in range(3):
			try:
				self.limiter.wait()
				with self.session.get(url, stream=True, timeout=120) as r:
					if r.status_code != 200:
						return False  
					with open(local, "wb") as fout:
						for chunk in r.iter_content(chunk_size=1 << 16):
							fout.write(chunk)
				if os.path.getsize(local) > 0:
					self.failures.pop(url, None)
					return True
			except Exception as e:
				self.failures[url] = "{}: {}".format(type(e).__name__, e)
				time.sleep(0.5)
		return False


class AssemblyDownloader(_GenomeDownloader): 
	BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
	SUFFIX = {"gff": "_genomic.gff.gz", "faa": "_protein.faa.gz",
			  "rna": "_rna_from_genomic.fna.gz", "genome": "_genomic.fna.gz"}

	def _partition_url(self, assembly: str) -> str:
		prefix, digits = assembly.split("_")
		digits = digits.split(".")[0]
		return "{}/{}/{}/{}/{}".format(
			self.BASE, prefix, digits[0:3], digits[3:6], digits[6:9])

	def _versioned_dir(self, assembly: str) -> Optional[str]:
		for _ in range(3):
			try:
				self.limiter.wait()
				r = self.session.get(self._partition_url(assembly) + "/", timeout=30)
				r.raise_for_status()
				self.failures.pop(assembly, None)
				for name in re.findall(r'href="([^"/]+)/"', r.text):
					if name.startswith(assembly):
						return name
				return None  
			except Exception as e:
				self.failures[assembly] = "{}: {}".format(type(e).__name__, e)
				time.sleep(0.5)
		return None

	def _fetch_one(self, assembly: str) -> GenomeFiles:
		vdir = self._versioned_dir(assembly)
		if vdir is None:
			return GenomeFiles()
		base = "{}/{}/{}".format(self._partition_url(assembly), vdir, vdir)
		return self._fetch_files({
			slot: (base + self.SUFFIX[slot],
				   os.path.join(self.out_dir, assembly + self.SUFFIX[slot]))
			for slot in self._slots()})


class MgnifyGenomeDownloader(_GenomeDownloader): 
	API = "https://www.ebi.ac.uk/metagenomics/api/v1/genomes/{}/downloads"
	ACCESSION_RE = re.compile(r"^MGYG\d+$")
	SUFFIX = {"gff": ".gff", "faa": ".faa", "genome": ".fna"}

	@classmethod
	def is_mgnify_accession(cls, assembly: str) -> bool:
		return bool(cls.ACCESSION_RE.match(assembly))

	def _fetch_one(self, assembly: str) -> GenomeFiles:
		try:
			self.limiter.wait()
			r = self.session.get(self.API.format(assembly), timeout=30)
			r.raise_for_status()
			urls = {f["id"]: f["links"]["self"] for f in r.json().get("data", [])
					if f.get("links", {}).get("self")}
			self.failures.pop(assembly, None)
		except Exception as e:
			self.failures[assembly] = "{}: {}".format(type(e).__name__, e)
			return GenomeFiles()
		return self._fetch_files({
			slot: (urls.get(assembly + self.SUFFIX[slot]),
				   os.path.join(self.out_dir, assembly + self.SUFFIX[slot]))
			for slot in self._slots()})


class LocalGenomeResolver: 
	GFF_EXT = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
	FAA_EXT = (".faa", ".faa.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz")
	RNA_EXT = (".rna.fna", ".rna.fna.gz", ".rna.fa", ".rna.fa.gz",
			   "_rna_from_genomic.fna", "_rna_from_genomic.fna.gz")
	GENOME_EXT = (".fna", ".fna.gz")

	def __init__(self, directory: str):
		self.genomes = self._pair_files(directory)
		self._protein_index = None

	def _pair_files(self, directory):
		found = {"gff": {}, "faa": {}, "rna": {}, "genome": {}}
		try:
			names = sorted(os.listdir(directory))
		except OSError:
			return {}
		claimed = set()

		for slot, exts in (("rna", self.RNA_EXT), ("gff", self.GFF_EXT),
						   ("faa", self.FAA_EXT), ("genome", self.GENOME_EXT)):
			for name in names:
				path = os.path.join(directory, name)
				if name in claimed or not os.path.isfile(path):
					continue
				base = self._basename(name, exts)
				if base is not None:
					found[slot][base] = path
					claimed.add(name)

		return {b: GenomeFiles(found["gff"][b], found["faa"][b],
							   found["rna"].get(b), found["genome"].get(b))
				for b in found["gff"] if b in found["faa"]}

	INFIX = ("_genomic", "_protein", "_rna_from_genomic", "_cds_from_genomic")

	@staticmethod
	def _basename(name, exts):
		for ext in sorted(exts, key=len, reverse=True):
			if name.endswith(ext):
				stem = name[:-len(ext)]
				for infix in LocalGenomeResolver.INFIX:
					if stem.endswith(infix):
						stem = stem[:-len(infix)]
						break
				return stem
		return None

	def _build_protein_index(self):
		index = {}
		for base, files in self.genomes.items():
			with NeighborhoodExtractor._open(files.faa) as fh:
				for rec in SeqIO.parse(fh, "fasta"):
					index.setdefault(rec.id, base)
		self._protein_index = index

	def resolve_pair(self, assembly: str) -> Optional[Tuple[str, GenomeFiles]]:
		if assembly in self.genomes:
			return assembly, self.genomes[assembly]
		for base in self.genomes:
			if base.startswith(assembly) or assembly.startswith(base):
				return base, self.genomes[base]
		return None

	def resolve_protein(self, protein: str) -> Optional[Tuple[str, GenomeFiles]]:
		if self._protein_index is None:
			self._build_protein_index()
		base = self._protein_index.get(protein)
		return (base, self.genomes[base]) if base else None


class FlankingGene(NamedTuple):
	accession: str   # protein/RNA accession, or '<biotype>*' for unidentified non-coding
	strand: str      # '+'/'-', normalized relative to the query gene
	start: int
	end: int
	product: str
	offset: int      # 0 = query, negative = upstream, positive = downstream
	query: str       # the query protein this neighbor belongs to
	is_rna: bool = False   # True for tRNA/rRNA/ncRNA genes (keep the RNA outline)
	contig: str = ""       # contig/sequence id the gene lies on (for locus matching)


class NeighborhoodExtractor: 
	def __init__(self, flank: int = 4, label_assembly: bool = False):
		self.flank = flank
		self.label_assembly = label_assembly
		self.sequences: Dict[str, str] = {}        # all flanking proteins: accession -> sequence
		self.query_sequences: Dict[str, str] = {}  # query proteins only: accession -> sequence
		self.row_sequences: Dict[str, str] = {}    # row id -> query sequence (one leaf per assembly-row)
		self.rna_sequences: Dict[str, str] = {}    # flanking RNAs: accession -> nucleotide sequence
		self.rna_products: Dict[str, str] = {}     # flanking RNAs: accession -> product name
		self.species: Dict[str, str] = {}          # row id -> organism name
		self.row_label: Dict[str, str] = {}        # row id -> display label
		self._gff_cache: Dict[str, List[dict]] = {}
		self._faa_cache: Dict[str, Dict[str, Tuple[str, str]]] = {}
		self._rna_cache: Dict[str, Dict[str, str]] = {}
		self._genome_cache: Tuple[Optional[str], Dict[str, str]] = (None, {})

	def extract(self, assembly: str, gff_path: str, faa_path: str,
				query: str, acceptable: Optional[set] = None,
				rna_path: Optional[str] = None,
				genome_path: Optional[str] = None) -> List[FlankingGene]:
		acceptable = acceptable or {query}
		genes = self._genes(assembly, gff_path)
		idx = next((i for i, g in enumerate(genes) if g["accession"] in acceptable), None)
		if idx is None:
			return []

		faa = self._faa(assembly, faa_path)
		rna = self._rna(assembly, rna_path) if rna_path else {}
		contig = genes[idx]["contig"]
		qstrand = genes[idx]["strand"]
		lo, hi = max(0, idx - self.flank), min(len(genes), idx + self.flank + 1)
		q_acc = genes[idx]["accession"]
		q_organism = faa.get(q_acc, (None, ""))[1]
		row_id = "{}|{}".format(query, assembly)
		# With -m > 1 a query contributes one row per assembly, so the accession alone
		# no longer identifies a row in the figures.
		name = "{}|{}".format(query, assembly) if self.label_assembly else query
		self.row_label[row_id] = "{}  {}".format(name, q_organism) if q_organism else name

		neighborhood = []
		for j in range(lo, hi):
			g = genes[j]
			if g["contig"] != contig:
				continue
			acc = g["accession"]
			if g["is_rna"]:
				self.rna_products[acc] = g["product"]
				rseq = rna.get(g["locus_tag"]) or rna.get(acc)
				if rseq is None and genome_path:

					rseq = self._slice_genome(
						self._genome(assembly, genome_path),
						g["contig"], g["start"], g["end"], g["strand"])
				if rseq is not None:
					self.rna_sequences[acc] = rseq
			else:
				seq, organism = faa.get(acc, (None, ""))
				if seq is not None:
					self.sequences[acc] = seq
					if j == idx:
						self.query_sequences[acc] = seq
						self.row_sequences[row_id] = seq   
						self.species[row_id] = organism
			offset = j - idx
			if qstrand == "-":
				offset = -offset
			neighborhood.append(FlankingGene(
				accession=acc,
				strand=self._norm_strand(qstrand, g["strand"]),
				start=g["start"], end=g["end"],
				product=g["product"],
				offset=offset,
				query=row_id,
				is_rna=g["is_rna"],
				contig=g["contig"],
			))
		return neighborhood

	@staticmethod
	def _open(path: str):
		if path.endswith(".gz"):
			return gzip.open(path, "rt", encoding="utf-8", errors="replace")
		return open(path, "rt", encoding="utf-8", errors="replace")

	def _genes(self, assembly: str, gff_path: str) -> List[dict]:
		if assembly in self._gff_cache:
			return self._gff_cache[assembly]
		genes = []
		with self._open(gff_path) as fh:
			for raw in fh:
				if raw.startswith("#"):
					continue
				col = raw.rstrip("\n").split("\t")
				if len(col) < 9:
					continue
				feature, attrs = col[2], col[8]
				if feature.endswith("gene"):
					genes.append(self._record(col, None, "",
						biotype=self._attr(attrs, "gene_biotype") or "",
						locus_tag=self._attr(attrs, "locus_tag") or "", is_rna=False))
				elif feature == "CDS":
					locus_tag = self._attr(attrs, "locus_tag") or ""
					accession = self._cds_accession(attrs, locus_tag)
					product = self._attr(attrs, "product") or ""
					if genes and genes[-1]["accession"] is None:
						genes[-1]["accession"] = accession
						genes[-1]["product"] = product
					else:

						genes.append(self._record(col, accession, product,
							biotype="protein_coding", locus_tag=locus_tag, is_rna=False))
				elif feature.endswith("RNA"):

					locus_tag = self._attr(attrs, "locus_tag") or ""
					accession = (self._attr(attrs, "Name") or self._attr(attrs, "transcript_id")
								 or self._strip_id_prefix(self._attr(attrs, "ID")) or locus_tag)
					product = self._attr(attrs, "product") or ""
					if genes and genes[-1]["accession"] is None:
						genes[-1]["accession"] = accession
						genes[-1]["product"] = product
						genes[-1]["is_rna"] = True
						if not genes[-1]["locus_tag"]:
							genes[-1]["locus_tag"] = locus_tag
					else:
						genes.append(self._record(col, accession, product,
							biotype=feature, locus_tag=locus_tag, is_rna=True))
		for g in genes:
			if not g["accession"]:
				g["accession"] = (g["biotype"] or "noProtein") + "*"

		genes.sort(key=lambda g: (g["contig"], g["start"]))
		self._gff_cache[assembly] = genes
		return genes

	@staticmethod
	def _record(col: List[str], accession: Optional[str], product: str,
				biotype: str, locus_tag: str, is_rna: bool) -> dict:
		return {
			"contig": col[0], "start": int(col[3]), "end": int(col[4]),
			"strand": col[6], "accession": accession, "product": product,
			"biotype": biotype, "locus_tag": locus_tag, "is_rna": is_rna,
		}

	def _faa(self, assembly: str, faa_path: str) -> Dict[str, Tuple[str, str]]:
		if assembly in self._faa_cache:
			return self._faa_cache[assembly]
		table = {}
		with self._open(faa_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				m = re.search(r"\[([^\]]+)\]\s*$", rec.description)
				organism = m.group(1) if m else ""
				table[rec.id] = (str(rec.seq), organism)
		self._faa_cache[assembly] = table
		return table

	def _rna(self, assembly: str, rna_path: str) -> Dict[str, str]:

		if assembly in self._rna_cache:
			return self._rna_cache[assembly]
		table = {}
		with self._open(rna_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				seq = str(rec.seq)
				m = re.search(r"\[locus_tag=([^\]]+)\]", rec.description)
				if m:
					table[m.group(1)] = seq
				table[rec.id] = seq
		self._rna_cache[assembly] = table
		return table

	def _genome(self, assembly: str, genome_path: str) -> Dict[str, str]:

		cached_for, table = self._genome_cache
		if cached_for == assembly:
			return table
		table = {}
		with self._open(genome_path) as fh:
			for rec in SeqIO.parse(fh, "fasta"):
				table[rec.id] = str(rec.seq)
		self._genome_cache = (assembly, table)
		return table

	@staticmethod
	def _slice_genome(genome_seqs: Dict[str, str], contig: str,
					   start: int, end: int, strand: str) -> Optional[str]:
		seq = genome_seqs.get(contig)
		if seq is None:
			return None
		sub = seq[start - 1:end] 
		if strand == "-":
			sub = str(Seq(sub).reverse_complement())
		return sub or None

	@staticmethod
	def _attr(attributes: str, key: str) -> Optional[str]:
		m = re.search(r"(?:^|;){}=([^;]*)".format(re.escape(key)), attributes)
		return m.group(1) if m else None

	@classmethod
	def _cds_accession(cls, attrs: str, locus_tag: str) -> Optional[str]:
		return (cls._attr(attrs, "protein_id")
				or cls._strip_id_prefix(cls._attr(attrs, "ID"))
				or locus_tag
				or cls._attr(attrs, "Name"))

	@staticmethod
	def _strip_id_prefix(id_attr: Optional[str]) -> Optional[str]:
		if not id_attr:
			return None
		for prefix in ("cds-", "rna-"):
			if id_attr.startswith(prefix):
				return id_attr[len(prefix):]
		return id_attr

	@staticmethod
	def _norm_strand(query_strand: str, gene_strand: str) -> str:
		if query_strand == "+":
			return gene_strand
		return "-" if gene_strand == "+" else "+"


class NeighborhoodClusterer:
	def __init__(self, iterations: int = 3, incE: float = 1e-3,
				 workers: Optional[int] = None):
		self.iterations = iterations
		self.incE = incE
		self.workers = workers
		self.alphabet = Alphabet.amino()
		self.adjacency: Dict[str, set] = {}

	def cluster(self, sequences: Dict[str, str]) -> List[List[str]]:
		if not sequences:
			return []

		digital = {name: TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
				   for name, seq in sequences.items()}
		block = DigitalSequenceBlock(self.alphabet, list(digital.values()))

		def search_one(item):
			name, query = item
			result = list(pyhmmer.hmmer.jackhmmer(
				[query], block,
				max_iterations=self.iterations,
				incE=self.incE,
				cpus=1,
			))[0]
			return name, {self._name(h) for h in result.hits if h.included}

		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			adjacency = dict(pool.map(search_one, digital.items()))

		self.adjacency = adjacency
		return self._connected_components(adjacency)

	@staticmethod
	def _name(hit) -> str:
		n = hit.name
		return n.decode() if isinstance(n, bytes) else n

	@staticmethod
	def _connected_components(adjacency: Dict[str, set]) -> List[List[str]]:
		seen, families = set(), []
		for node in adjacency:
			if node in seen:
				continue
			stack, component = [node], set()
			while stack:
				x = stack.pop()
				if x in component:
					continue
				component.add(x)
				seen.add(x)
				stack.extend(adjacency.get(x, set()) - component)
			families.append(sorted(component))
		families.sort(key=len, reverse=True)
		return families


class RnaClusterer: 
	def __init__(self, incE: float = 1e-3, workers: Optional[int] = None):
		self.incE = incE
		self.workers = workers
		self.alphabet = Alphabet.dna()
		self.adjacency: Dict[str, set] = {}

	def cluster(self, sequences: Dict[str, str]) -> List[List[str]]:
		if not sequences:
			return []
		digital = {name: TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
				   for name, seq in sequences.items()}
		block = DigitalSequenceBlock(self.alphabet, list(digital.values()))

		def search_one(item):
			name, query = item
			hits = list(pyhmmer.hmmer.nhmmer([query], block, incE=self.incE, cpus=1))[0]
			return name, {NeighborhoodClusterer._name(h) for h in hits if h.included}

		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			adjacency = dict(pool.map(search_one, digital.items()))
		self.adjacency = adjacency
		return NeighborhoodClusterer._connected_components(adjacency)

	@staticmethod
	def cluster_by_name(products: Dict[str, str]) -> List[List[str]]:
		groups: Dict[str, List[str]] = {}
		for acc, product in products.items():
			key = RnaClusterer._normalise(product) or acc
			groups.setdefault(key, []).append(acc)
		families = [sorted(v) for v in groups.values()]
		families.sort(key=len, reverse=True)
		return families

	@staticmethod
	def _normalise(product: str) -> str:
		return " ".join(product.lower().split())


class ReportWriter: 
	def __init__(self, neighborhoods, families, species,
				 queries, protein_to_assemblies, matched,
				 order=None, adjacency=None):
		self.neighborhoods = neighborhoods
		self.families = families
		self.species = species
		self.queries = queries
		self.protein_to_assemblies = protein_to_assemblies
		self.matched = matched
		self.adjacency = adjacency or {}
		rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
		self.fam_of = family_numbers(families, rna_accessions)
		self.by_query = {}
		for g in neighborhoods:
			self.by_query.setdefault(g.query, []).append(g)
		if order:
			rank = {row: i for i, row in enumerate(order)}
			self.by_query = {row: self.by_query[row] for row in
							 sorted(self.by_query, key=lambda r: rank.get(r, len(rank)))}
		self.occurrences = Counter(g.accession for g in neighborhoods)
		self.products = {}
		for g in neighborhoods:
			self.products.setdefault(g.accession, g.product)

	def write_all(self, out_path):
		self.operon_tsv(out_path("_operon.tsv"))
		self.clusters_tsv(out_path("_clusters.tsv"))
		self.outdesc_txt(out_path("_outdesc.txt"))
		self.species_info(out_path("_speciesInfo.txt"))
		self.query_status(out_path("_QueryStatus.txt"))
		self.flankgene_report(out_path("_flankgene_Report.log"))
		if self.adjacency:
			self.jackhits_tsv(out_path("_jackhits.tsv"))
		return self.accession_issues(out_path("_accessionIssues.txt"))

	@staticmethod
	def _split_row(row_id):
		query, _, assembly = row_id.partition("|")
		return query, assembly

	def operon_tsv(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\tfamily\tstrand\toffset\t"
					  "start\tend\tlength\tcontig\taccession\tproduct\n")
			for row_id in self.by_query:
				query, assembly = self._split_row(row_id)
				sp = self.species.get(row_id, "")
				for g in sorted(self.by_query[row_id], key=lambda x: x.offset):
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
						query, assembly or "-", sp, self.fam_of.get(g.accession, "-"),
						g.strand, g.offset, g.start, g.end, g.end - g.start + 1,
						g.contig or "-", g.accession, g.product))

	def clusters_tsv(self, path):
		with open(path, "w") as out:
			out.write("#family\tsize\tmembers\n")
			for fam in self.families:
				label = self.fam_of.get(fam[0], "-") if len(fam) > 1 else "-"
				out.write("{}\t{}\t{}\n".format(label, len(fam), ",".join(fam)))

	def outdesc_txt(self, path):
		with open(path, "w") as out:
			for fam in self.families:
				if len(fam) < 2:
					continue
				label = self.fam_of.get(fam[0], "-")
				for acc in sorted(fam, key=lambda a: -self.occurrences.get(a, 0)):
					out.write("{}({})\t{}\t{}\n".format(
						label, self.occurrences.get(acc, 0), acc,
						self.products.get(acc, "")))
				out.write("\n\n")

	def jackhits_tsv(self, path):
		with open(path, "w") as out:
			out.write("#accession\tfamily\tn_hits\thits\n")
			for acc in sorted(self.adjacency):
				hits = sorted(self.adjacency[acc])
				out.write("{}\t{}\t{}\t{}\n".format(
					acc, self.fam_of.get(acc, "-"), len(hits), ";".join(hits)))

	def species_info(self, path):
		with open(path, "w") as out:
			out.write("#query\tassembly\tspecies\n")
			for row_id in sorted(self.species):
				query, assembly = self._split_row(row_id)
				out.write("{}\t{}\t{}\n".format(query, assembly or "-",
											    self.species[row_id]))

	def query_status(self, path):
		with open(path, "w") as out:
			out.write("#query\tassemblies\tflanking_genes_found\n")
			for q in self.queries:
				asms = self.protein_to_assemblies.get(q, [])
				status = "Yes" if q in self.matched else "No"
				out.write("{}\t{}\t{}\n".format(q, ";".join(asms) if asms else "-", status))

	def flankgene_report(self, path):
		with open(path, "w") as out:
			for query in self.by_query:
				genes = sorted(self.by_query[query], key=lambda x: x.offset)
				chain = " ".join("{}({})".format(g.accession, self.fam_of.get(g.accession, "-"))
								 for g in genes)
				out.write("{}\t{}\n".format(query, chain))

	def accession_issues(self, path):
		lines = []
		for q in self.queries:
			if not self.protein_to_assemblies.get(q):
				lines.append("{}\tno assembly resolved (not found locally or via NCBI)".format(q))
			elif q not in self.matched:
				lines.append("{}\tassembly resolved but no flanking neighborhood extracted".format(q))
		with open(path, "w") as out:
			out.write("#query\tissue\n")
			for line in lines:
				out.write(line + "\n")
		return len(lines)

VERSION = "2.2.1"

def main():
	usage = ''' Description: Identify flanking genes and cluster them based on similarity. Requirement= Python3, BioPython, pyhmmer, requests. '''
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument("-i", "--input_list", required=True, help=" Protein Accession eg. WP_047256880.1, optionally tab-separated with an assembly Identifier eg. GCF_000001765.3 (NCBI RefSeq/GenBank) or MGYG000454827 (MGnify Genomes catalogue -- the protein accession must then be the exact locus tag used in that genome's annotation). One per line. ")
	parser.add_argument("-u", "--user_email", required=True, help=" User Email Address (required by NCBI Entrez). ")
	parser.add_argument("-api", "--api_key", help=" NCBI API Key. ")
	parser.add_argument("-g", "--gene", type=int, default=4, help=" Number of flanking genes up/downstream. Default = 4 ")
	parser.add_argument("-m", "--max_assemblies", type=int, default=1, help=" Max assemblies per protein. Default = 1 ")
	parser.add_argument("-e", "--ethreshold", type=float, default=1e-3, help=" Jackhmmer inclusion E-value threshold for clustering flanking genes. Default = 1e-3 ")
	parser.add_argument("-n", "--number", type=int, default=3, help=" Number of jackhmmer iterations for clustering. Default = 3 ")
	parser.add_argument("-c", "--cpu", type=int, help=" Max parallel CPU workers (default: auto-detect). ")
	parser.add_argument("-ts", "--tshape", type=int, default=20, help=" Size of the flanking-gene triangles in the tree view. Default = 20 ")
	parser.add_argument("-tf", "--tfontsize", type=int, default=13, help=" Font size inside the tree-view triangles. Default = 13 ")
	parser.add_argument("-tmp", "--temporary", default="./genomes", help=" Temporary directory for downloaded assemblies; deleted at the end. Default = ./genomes ")
	parser.add_argument("-O", "--output", default="output", help=" Directory for result files; its name is also the file prefix. A YYYYMMDD_HHMMSS stamp of the run start is appended, so repeated runs do not overwrite each other. Default = output ")
	parser.add_argument("--no_timestamp", action="store_true", help=" Use -O/--output verbatim instead of appending a date-time stamp. Repeated runs then overwrite each other; useful for scripted pipelines that need a fixed path. ")
	parser.add_argument("--tree", action="store_true", help=" Also build a phylogenetic tree with aligned neighbourhood triangles (<dir>_tree.svg). Does not affect the neighbours output. ")
	parser.add_argument("-to", "--tree_order", action="store_true", help=" Order the neighbours output by tree leaf order (implies --tree). ")
	parser.add_argument("--domains", action="store_true", help=" Scan flanking proteins for domains and write <dir>_domains.svg (requires --hmmdb). ")
	parser.add_argument("--hmmdb", default="./pfam_db/Pfam-A.hmm", help=" HMM database file (e.g. Pfam) for domain scanning. ")
	parser.add_argument("--clans", help=" Pfam-A.clans.tsv(.gz): colour domains by clan instead of family. ")
	parser.add_argument("--tmhmm", action="store_true", help=" Predict transmembrane regions (DeepTMHMM, via BioLib cloud) and draw them as double red dotted lines in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("--signalp", action="store_true", help=" Predict signal peptides (SignalP-6, via BioLib cloud) and draw them as black triangles in the domain figure. Off by default; needs pybiolib and network. ")
	parser.add_argument("--sismis", action="store_true", help=" Scan each assembly's genomic FASTA for secretion systems using Sismis (github.com/lmc297/Sismis) and write <dir>_secretion.tsv, noting which query neighborhoods (if any) each hit overlaps. Downloads the genomic FASTA per assembly. Off by default; needs sismis installed (pip install sismis). ")
	parser.add_argument("-k", "--keep", action="store_true", help=" Keep the downloaded assemblies instead of deleting the temporary directory at the end. ")
	parser.add_argument("--use_local", metavar="DIR", help=" Directory of local .gff/.faa genome files to search before falling back to NCBI. Files may be gzipped; a genome is a .gff and .faa sharing a basename. ")
	parser.add_argument("--cluster_rna", action="store_true", help=" Cluster flanking RNA genes into families too (off by default). Uses RNA sequences (nhmmer) when available, otherwise groups by product name. ")
	parser.add_argument("-vb", "--verbose", action="store_true", help=" Print progress and a stage-by-stage funnel. ")
	parser.add_argument("-v", "--version", action="version", version="FlaGs2 " + VERSION)
	args = parser.parse_args(translate_legacy_args(sys.argv[1:]))

	if not os.path.isfile(args.input_list):
		sys.exit("Error: input list not found: {}".format(args.input_list))
	if args.use_local and not os.path.isdir(args.use_local):
		sys.exit("Error: --use_local directory not found: {}".format(args.use_local))
	if args.domains and args.hmmdb and not os.path.isfile(args.hmmdb):
		sys.exit("Error: --hmmdb file not found: {}".format(args.hmmdb))
	if args.clans and not os.path.isfile(args.clans):
		sys.exit("Error: --clans file not found: {}".format(args.clans))

	if not args.no_timestamp:
		args.output = "{}_{}".format(os.path.normpath(args.output),
									 time.strftime("%Y%m%d_%H%M%S"))

	args.output = os.path.abspath(args.output)
	args.temporary = os.path.abspath(args.temporary)

	timings = {}
	t_start = time.perf_counter()
	t0 = t_start

	proteins_assembly, proteins_only = AccessionListReader(args.input_list).read()
	timings["1_read_input"] = time.perf_counter() - t0; t0 = time.perf_counter()
	all_queries = [p for p, _ in proteins_assembly] + list(proteins_only)
	if args.verbose:
		print(">> read {} queries ({} paired, {} protein-only)".format(
			len(all_queries), len(proteins_assembly), len(proteins_only)), flush=True)

	local = LocalGenomeResolver(args.use_local) if args.use_local else None
	local_files: Dict[str, GenomeFiles] = {}   # source id -> files
	protein_to_assemblies = {}                 # protein -> [source ids]
	local_acceptable = {}                      # protein -> {source id: {accession}}
	pending_pairs = []                         # (protein, assembly) still needing NCBI
	pending_only = []                          # bare proteins still needing NCBI

	if local:
		for protein, assembly in proteins_assembly:
			hit = local.resolve_pair(assembly)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_pairs.append([protein, assembly])
		for protein in proteins_only:
			hit = local.resolve_protein(protein)
			if hit:
				base, files = hit
				local_files[base] = files
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_only.append(protein)
		if args.verbose:
			print(">> local: {} genomes indexed; resolved {} of {} queries locally".format(
				len(local.genomes), len(protein_to_assemblies),
				len(proteins_assembly) + len(proteins_only)), flush=True)
	else:
		pending_pairs = proteins_assembly
		pending_only = proteins_only

	mapper = ProteinAssemblyMapper(email=args.user_email, api_key=args.api_key,
								   max_assemblies=args.max_assemblies)
	ncbi_map = mapper.map(pending_only)
	for protein, asms in ncbi_map.items():
		protein_to_assemblies.setdefault(protein, []).extend(asms)
	for protein, assembly in pending_pairs:
		protein_to_assemblies.setdefault(protein, []).append(assembly)
	timings["2_ipg_mapping"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose and pending_only:
		print(">> NCBI IPG: mapped {} of {} remaining proteins to assemblies".format(
			len(ncbi_map), len(pending_only)), flush=True)

	assemblies = sorted({asm for asms in protein_to_assemblies.values()
						 for asm in asms if asm not in local_files})
	mgnify_assemblies = [a for a in assemblies if MgnifyGenomeDownloader.is_mgnify_accession(a)]
	ncbi_assemblies = [a for a in assemblies if a not in set(mgnify_assemblies)]
	dl_workers = args.cpu if args.cpu else min(max(len(assemblies), 1), 10)
	dl_rate = 10.0 if args.api_key else 5.0

	def progress(label):
		if not args.verbose:
			return None
		return lambda done, total: print(
			">> {} download: {}/{}".format(label, done, total), flush=True)

	if args.verbose and assemblies:
		print(">> downloading {}...".format(plural(len(assemblies), "genome")), flush=True)
	downloaded: Dict[str, GenomeFiles] = {}
	failures: Dict[str, str] = {}
	if ncbi_assemblies:
		dl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
								want_rna=args.cluster_rna, want_genome=args.sismis)
		downloaded.update(dl.download_many(ncbi_assemblies, progress("NCBI")))
		failures.update(dl.failures)
		timings["3a_download_ncbi"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if mgnify_assemblies:
		mg = MgnifyGenomeDownloader(out_dir=args.temporary, workers=dl_workers, rate=dl_rate,
									want_genome=args.sismis or args.cluster_rna)
		downloaded.update(mg.download_many(mgnify_assemblies, progress("MGnify")))
		failures.update(mg.failures)
		timings["3b_download_mgnify"] = time.perf_counter() - t0; t0 = time.perf_counter()
	downloaded.update(local_files)
	if args.verbose:
		ready = sum(1 for f in downloaded.values() if f.gff and f.faa)
		print(">> genomes ready: {} of {} usable ({} NCBI, {} MGnify, {} local)".format(
			ready, len(downloaded), len(ncbi_assemblies), len(mgnify_assemblies),
			len(local_files)), flush=True)
		if failures:

			print(">> {}, first few:".format(plural(len(failures), "download error")), flush=True)
			for key, msg in list(failures.items())[:5]:
				print("     {}: {}".format(key.rsplit("/", 1)[-1] or key, msg), flush=True)

	extractor = NeighborhoodExtractor(flank=args.gene,
									  label_assembly=args.max_assemblies > 1)
	all_neighborhoods = []
	matched = set()  
	for protein, asms in protein_to_assemblies.items():
		for asm in asms:
			files = downloaded.get(asm, GenomeFiles())
			if not (files.gff and files.faa):
				continue
			if asm in local_files:
				acceptable = local_acceptable.get(protein, {}).get(asm)
			else:
				acceptable = mapper.accessions_in.get(protein, {}).get(asm)
			rows = extractor.extract(
				asm, files.gff, files.faa, protein, acceptable,
				rna_path=files.rna if args.cluster_rna else None,
				genome_path=files.genome if args.cluster_rna else None)
			if rows:
				matched.add(protein)
				all_neighborhoods.extend(rows)
	timings["4_extract_neighbors"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		print(">> extracted {} flanking-gene records; {} of {} queries matched".format(
			len(all_neighborhoods), len(matched), len(all_queries)), flush=True)

	sismis_mod = None
	sismis_hits: List = []
	sismis_rows: Dict[str, Tuple[str, str, int, int]] = {}
	sismis_statuses: Dict[str, str] = {}
	features = {}

	def _run_sismis():
		nonlocal sismis_mod
		t = time.perf_counter()
		rows, hits, statuses = {}, [], {}
		if not (args.sismis and all_neighborhoods):
			return hits, rows, statuses, 0.0
		try:
			import flags2_secretion as mod
		except ImportError:
			print("Warning: --sismis needs the sismis package (pip install sismis); skipping secretion-system detection.")
			return hits, rows, statuses, time.perf_counter() - t
		for g in all_neighborhoods:
			assembly = g.query.rsplit("|", 1)[-1]
			if g.query in rows:
				_, contig, lo, hi = rows[g.query]
				rows[g.query] = (assembly, contig, min(lo, g.start), max(hi, g.end))
			else:
				rows[g.query] = (assembly, g.contig, g.start, g.end)
		scanner = mod.SismisScanner(out_dir=os.path.join(args.output, "sismis"))
		for assembly in sorted({asm for asm, _, _, _ in rows.values()}):
			genome_path = downloaded.get(assembly, GenomeFiles()).genome
			if not genome_path:
				statuses[assembly] = "skipped: no genomic FASTA downloaded"
				continue
			try:
				found = scanner.scan_assembly(assembly, genome_path)
				hits.extend(found)
				statuses[assembly] = ("{} secretion system(s) predicted".format(len(found))
									   if found else "no secretion system predicted")
			except Exception as e:
				statuses[assembly] = "error: {}".format(e)
		sismis_mod = mod
		return hits, rows, statuses, time.perf_counter() - t

	def _run_tmhmm():
		t = time.perf_counter()
		if not args.tmhmm:
			return {}, 0.0
		try:
			import flags2_features as feat_mod
			tm = feat_mod.TMScanner().scan(extractor.sequences, want_signal=not args.signalp)
		except ImportError:
			print("Warning: --tmhmm needs pybiolib (pip install pybiolib); skipping transmembrane prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: DeepTMHMM did not finish, skipping transmembrane regions ({}).".format(e))
			return {}, time.perf_counter() - t
		return tm, time.perf_counter() - t

	def _run_signalp():
		t = time.perf_counter()
		if not args.signalp:
			return {}, 0.0
		try:
			import flags2_features as feat_mod
			sp = feat_mod.SignalPScanner().scan(extractor.sequences)
		except ImportError:
			print("Warning: --signalp needs pybiolib (pip install pybiolib); skipping signal-peptide prediction.")
			return {}, time.perf_counter() - t
		except Exception as e:
			print("Warning: SignalP did not finish, skipping signal peptides ({}).".format(e))
			return {}, time.perf_counter() - t
		return sp, time.perf_counter() - t

	tm, sp = {}, {}
	active = [n for n, flag in (("sismis", args.sismis), ("tmhmm", args.tmhmm),
								 ("signalp", args.signalp)) if flag]
	if active:
		if args.verbose:
			print(">> running {} in the background (cloud/subprocess, not local CPU)...".format(
				", ".join(active)), flush=True)
		if args.tmhmm or args.signalp:
			import flags2_features as feat_mod
			feat_mod.warm_up()
		task = {"sismis": _run_sismis, "tmhmm": _run_tmhmm, "signalp": _run_signalp}
		with ThreadPoolExecutor(max_workers=3) as pool:
			futures = {pool.submit(task[name]): name for name in active}
			for fut in as_completed(futures):
				name = futures[fut]
				if name == "sismis":
					sismis_hits, sismis_rows, sismis_statuses, elapsed = fut.result()
					timings["sismis_scan"] = elapsed
					if sismis_mod and args.verbose:
						print(">> sismis: scanned {}, found {}".format(
							plural(len(sismis_statuses), "assembly", "assemblies"),
							plural(len(sismis_hits), "secretion system")), flush=True)
				elif name == "tmhmm":
					tm, elapsed = fut.result()
					timings["tmhmm_scan"] = elapsed
					if args.verbose:
						print(">> DeepTMHMM: features on {} proteins".format(len(tm)), flush=True)
				elif name == "signalp":
					sp, elapsed = fut.result()
					timings["signalp_scan"] = elapsed
					if args.verbose:
						print(">> SignalP: signal peptides on {} proteins".format(len(sp)), flush=True)
		for acc, regs in tm.items():
			features.setdefault(acc, []).extend(regs)
		for acc, regs in sp.items():
			features.setdefault(acc, []).extend(regs)

	if not args.keep:
		shutil.rmtree(args.temporary, ignore_errors=True)

	if not all_neighborhoods:
		os.makedirs(args.output, exist_ok=True)
		prefix = os.path.basename(os.path.normpath(args.output))
		reporter = ReportWriter(all_neighborhoods, [], extractor.species,
								all_queries, protein_to_assemblies, matched)
		issues_path = os.path.join(args.output, prefix + "_accessionIssues.txt")
		reporter.accession_issues(issues_path)
		reporter.query_status(os.path.join(args.output, prefix + "_QueryStatus.txt"))
		sys.exit("No flanking neighbourhoods could be extracted for any query. "
				 "See {} for per-query details.".format(issues_path))

	t0 = time.perf_counter()
	if args.verbose:
		print(">> clustering {} flanking proteins...".format(len(extractor.sequences)), flush=True)
	clusterer = NeighborhoodClusterer(iterations=args.number, incE=args.ethreshold,
									  workers=args.cpu)
	families = clusterer.cluster(extractor.sequences)

	rna_families = []
	if args.cluster_rna:
		rna = RnaClusterer(incE=args.ethreshold, workers=args.cpu)
		have_seq = set(extractor.rna_sequences)
		all_rna = set(extractor.rna_products)
		missing = all_rna - have_seq
		rna_families = rna.cluster(extractor.rna_sequences)
		if missing:
			fallback = {acc: extractor.rna_products[acc] for acc in missing}
			rna_families += rna.cluster_by_name(fallback)
			print("Warning: {} of {} flanking RNAs had no nucleotide sequence available, "
				  "so they were grouped by product name rather than by sequence."
				  .format(len(missing), len(all_rna)))
		families = families + rna_families
	timings["5_clustering"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		msg = ">> clustered {} into {}".format(
			plural(len(extractor.sequences), "flanking protein"),
			plural(len(families) - len(rna_families), "family", "families"))
		if args.cluster_rna:
			msg += "; {} into {}".format(
				plural(len(extractor.rna_products), "RNA"),
				plural(len(rna_families), "family", "families"))
		print(msg, flush=True)

	want_tree = args.tree or args.tree_order
	newick = ""
	leaf_order = None
	tree_mod = None
	if want_tree:
		if args.verbose:
			print(">> building tree ({} leaves)...".format(len(extractor.row_sequences)), flush=True)
		import flags2_tree as tree_mod
		builder = tree_mod.TreeBuilder(threads=0)
		newick, _ = builder.build(extractor.row_sequences)
		if newick:
			leaf_order = tree_mod.ladderized_leaf_order(newick)
	if want_tree:
		timings["6_tree"] = time.perf_counter() - t0
	t0 = time.perf_counter()
	if args.verbose and want_tree:
		print(">> tree: {}".format("built ({} leaves)".format(len(leaf_order))
			  if newick else "skipped (fewer than 3 query sequences)"), flush=True)

	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	def out_path(suffix):
		return os.path.join(args.output, prefix + suffix)
	order = leaf_order if args.tree_order else None
	operon = OperonView(mode="families")
	with open(out_path("_neighbors.svg"), "w") as out:
		out.write(operon.render(all_neighborhoods, families, extractor.species, order,
								labels=extractor.row_label))

	tree_written = False
	if newick:
		viz = tree_mod.NeighborhoodVisualizer(gene_h=args.tshape, font=args.tfontsize,
											  row_h=max(args.tshape + 4, args.tfontsize + 4))
		with open(out_path("_tree.svg"), "w") as out:
			out.write(viz.render(newick, all_neighborhoods, families,
								 extractor.species, labels=extractor.row_label))
		with open(out_path("_tree.nwk"), "w") as out:
			out.write(newick + "\n")
		tree_written = True
	timings["7_visualize"] = time.perf_counter() - t0; t0 = time.perf_counter()
	domains_written = False
	domain_table_written = False

	want_domain_fig = args.domains or features
	if want_domain_fig:
		domains = {}
		clan_map = None
		if args.domains and not args.hmmdb:
			print("Warning: --domains needs --hmmdb; drawing the figure without domains.")
		elif args.domains:
			try:
				import flags2_domains as dom_mod
				if args.verbose:
					print(">> scanning {} proteins for domains...".format(
						len(extractor.sequences)), flush=True)
				scanner = dom_mod.DomainScanner(args.hmmdb, cpus=args.cpu or 0)
				domains = scanner.scan(extractor.sequences)
				clan_map = dom_mod.DomainScanner.load_clans(args.clans) if args.clans else None
			except Exception as e:
				print("Warning: could not read the HMM database, drawing the figure without domains ({}).".format(e))
		if domains:
			try:
				dom_mod.DomainScanner.write_report(
					domains, out_path("_domains.tsv"), clans=clan_map,
					families=family_numbers(families,
											{g.accession for g in all_neighborhoods if g.is_rna}))
				domain_table_written = True
			except Exception as e:
				print("Warning: could not write the domain table ({}).".format(e))
		try:
			dom_viz = OperonView(mode="domains")
			with open(out_path("_domains.svg"), "w") as out:
				out.write(dom_viz.render(all_neighborhoods, families,
										 extractor.species, order,
										 domains=domains, clans=clan_map,
										 features=features, labels=extractor.row_label))
			domains_written = True
		except Exception as e:
			print("Warning: could not write the domain figure ({}).".format(e))
		timings["8_domains"] = time.perf_counter() - t0
	t0 = time.perf_counter()

	secretion_written = False
	if args.sismis and sismis_mod:
		t0 = time.perf_counter()
		matches = sismis_mod.match_rows(sismis_hits, sismis_rows)
		sismis_mod.write_report(sismis_hits, matches, out_path("_secretion.tsv"))
		sismis_mod.write_diagnostics(sismis_statuses, out_path("_sismis_diagnostics.txt"))
		try:
			sis_viz = OperonView(mode="secretion")
			with open(out_path("_secretion.svg"), "w") as out:
				out.write(sis_viz.render(all_neighborhoods, families,
										 extractor.species, order,
										 secretion=sismis_hits, labels=extractor.row_label))
			secretion_written = True
		except Exception as e:
			print("Warning: could not write the secretion figure ({}).".format(e))
		timings["9_sismis_report"] = time.perf_counter() - t0

	adjacency = dict(clusterer.adjacency)
	if args.cluster_rna:
		adjacency.update(rna.adjacency)
	reporter = ReportWriter(all_neighborhoods, families, extractor.species,
							all_queries, protein_to_assemblies, matched,
							order=order, adjacency=adjacency)
	n_issues = reporter.write_all(out_path)
	if args.verbose:
		print(">> wrote data tables and reports ({} queries with issues)".format(n_issues),
			  flush=True)

	print("\n{} -> {}".format(
		plural(len(extractor.sequences), "flanking protein"),
		plural(len(families) - len(rna_families), "family", "families")))
	print("\noutputs in {}/".format(os.path.relpath(args.output)
									if args.output.startswith(os.getcwd() + os.sep)
									else args.output))
	print("  {}_neighbors.svg".format(prefix))
	if tree_written:
		print("  {}_tree.svg / {}_tree.nwk".format(prefix, prefix))
	elif want_tree:
		print("  (tree skipped: fewer than 3 query sequences)")
	if domains_written:
		print("  {}_domains.svg".format(prefix))
	if domain_table_written:
		print("  {}_domains.tsv".format(prefix))
	if secretion_written:
		print("  {}_secretion.svg".format(prefix))
	if args.sismis and sismis_mod:
		print("  {}_secretion.tsv / {}_sismis_diagnostics.txt".format(prefix, prefix))
	for suffix in ("_operon.tsv", "_clusters.tsv", "_outdesc.txt", "_speciesInfo.txt",
				   "_QueryStatus.txt", "_flankgene_Report.log", "_jackhits.tsv",
				   "_accessionIssues.txt"):
		print("  {}{}".format(prefix, suffix))

	if args.verbose:
		print("\n--- timing (seconds) ---")
		for stage in sorted(timings):
			print("  {:24s} {:8.2f}".format(stage, timings[stage]))
		print("  {:24s} {:8.2f}".format("TOTAL", time.perf_counter() - t_start))


if __name__ == '__main__':
	try:
		main()
	except FileNotFoundError as e:
		sys.exit("Error: file not found - {}".format(e))
	except KeyboardInterrupt:
		sys.exit("\nInterrupted.")
	except Exception as e:
		sys.exit("Error: {}".format(e))