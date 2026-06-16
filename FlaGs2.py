import argparse
import colorsys
import gzip
import os
import re
import shutil
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Tuple, NamedTuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import Entrez, SeqIO
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock

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
			out += ["-i", input_files[0]]      # no rewrite needed
		else:
			merged = "_flags2_merged_input.txt"
			with open(merged, "w") as mf:
				for path in input_files:        # bare proteins: copy as-is
					with open(path) as f:
						for line in f:
							if line.strip():
								mf.write(line.strip() + "\n")
				for path in paired_files:        # two-column: swap to protein-first
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
def family_numbers(families):
	number = {}
	for i, fam in enumerate([f for f in families if len(f) > 1], 1):
		for acc in fam:
			number[acc] = i
	return number


class AccessionListReader: #reads the input file and classifies entries
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
						print("Skipping malformed line {}: {!r}".format(n, line))
				else:
					proteins_only.append(line)
		return proteins_assembly, proteins_only


class ProteinAssemblyMapper: #maps protein accessions -> assemblies via NCBI IPG or BioProject
	def __init__(self, email: str, max_assemblies: int = 5,
				 api_key: Optional[str] = None, ncbi_time: float = 0.4):
		self.max_assemblies = max_assemblies
		self.ncbi_time = ncbi_time
		# {query: {assembly: {accession, ...}}}
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
		# prefer GCF (RefSeq) over GCA (GenBank), then alphabetical, before capping
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


class AssemblyDownloader: #downloads GFF + FASTA from NCBI over HTTPS
	BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
	WANT = ("_genomic.gff.gz", "_protein.faa.gz")

	def __init__(self, out_dir: Optional[str] = None, workers: int = 16):
		self.out_dir = out_dir or tempfile.gettempdir()
		self.workers = workers
		os.makedirs(self.out_dir, exist_ok=True)
		self.session = requests.Session()
		retry = Retry(total=5, backoff_factor=0.5,
					  status_forcelist=[429, 500, 502, 503, 504],
					  allowed_methods=frozenset(["GET"]))
		adapter = HTTPAdapter(max_retries=retry,
							  pool_connections=workers, pool_maxsize=workers)
		self.session.mount("https://", adapter)

	def download_many(self, assemblies: List[str]) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
		if not assemblies:
			return {}
		with ThreadPoolExecutor(max_workers=self.workers) as pool:
			results = pool.map(self._fetch_one, assemblies)
		return dict(zip(assemblies, results))

	def _partition_url(self, assembly: str) -> str:
		prefix, digits = assembly.split("_")
		digits = digits.split(".")[0]
		return "{}/{}/{}/{}/{}".format(
			self.BASE, prefix, digits[0:3], digits[3:6], digits[6:9])

	def _versioned_dir(self, assembly: str) -> Optional[str]:
		for _ in range(3):
			try:
				r = self.session.get(self._partition_url(assembly) + "/", timeout=30)
				r.raise_for_status()
				for name in re.findall(r'href="([^"/]+)/"', r.text):
					if name.startswith(assembly):
						return name
				return None  # listing succeeded but no matching dir: genuinely absent
			except Exception:
				time.sleep(0.5)
		return None

	def _fetch_one(self, assembly: str) -> Tuple[Optional[str], Optional[str]]:
		vdir = self._versioned_dir(assembly)
		if vdir is None:
			return None, None
		base = "{}/{}/{}".format(self._partition_url(assembly), vdir, vdir)
		paths = [None, None]
		for i, suffix in enumerate(self.WANT):
			local = os.path.join(self.out_dir, assembly + suffix[8:])
			if self._stream(base + suffix, local):
				paths[i] = local
		return paths[0], paths[1]

	def _stream(self, url: str, local: str) -> bool:
		for _ in range(3):
			try:
				with self.session.get(url, stream=True, timeout=120) as r:
					if r.status_code != 200:
						return False   # 404 etc: file genuinely absent, don't retry
					with open(local, "wb") as fout:
						for chunk in r.iter_content(chunk_size=1 << 16):
							fout.write(chunk)
				if os.path.getsize(local) > 0:
					return True
			except Exception:
				time.sleep(0.5)
		return False


class LocalGenomeResolver: #resolves proteins to local .gff/.faa genomes before NCBI
	GFF_EXT = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
	FAA_EXT = (".faa", ".faa.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz")

	def __init__(self, directory: str):
		self.genomes = self._pair_files(directory)
		self._protein_index = None

	def _pair_files(self, directory):
		gff, faa = {}, {}
		try:
			names = os.listdir(directory)
		except OSError:
			return {}
		for name in names:
			path = os.path.join(directory, name)
			if not os.path.isfile(path):
				continue
			base = self._basename(name, self.GFF_EXT)
			if base is not None:
				gff[base] = path
				continue
			base = self._basename(name, self.FAA_EXT)
			if base is not None:
				faa[base] = path
		return {b: (gff[b], faa[b]) for b in gff if b in faa}

	@staticmethod
	def _basename(name, exts):
		for ext in sorted(exts, key=len, reverse=True):
			if name.endswith(ext):
				return name[:-len(ext)]
		return None

	def _build_protein_index(self):
		index = {}
		for base, (_, faa_path) in self.genomes.items():
			with NeighborhoodExtractor._open(faa_path) as fh:
				for rec in SeqIO.parse(fh, "fasta"):
					index.setdefault(rec.id, base)
		self._protein_index = index

	def resolve_pair(self, assembly: str):
		if assembly in self.genomes:
			return (assembly,) + self.genomes[assembly]
		for base in self.genomes:
			if base.startswith(assembly) or assembly.startswith(base):
				return (base,) + self.genomes[base]
		return None

	def resolve_protein(self, protein: str):
		if self._protein_index is None:
			self._build_protein_index()
		base = self._protein_index.get(protein)
		if base is None:
			return None
		return (base,) + self.genomes[base]


class FlankingGene(NamedTuple):
	accession: str   # protein accession, or '<biotype>*' for pseudo/non-coding
	strand: str      # '+'/'-', normalized relative to the query gene
	start: int
	end: int
	product: str
	offset: int      # 0 = query, negative = upstream, positive = downstream
	query: str       # the query protein this neighbor belongs to


class NeighborhoodExtractor: #parses GFF+FAA, extracts flanking genes around a query
	def __init__(self, flank: int = 4):
		self.flank = flank
		self.sequences: Dict[str, str] = {}        # all flanking proteins: accession -> sequence
		self.query_sequences: Dict[str, str] = {}  # query proteins only: accession -> sequence
		self.species: Dict[str, str] = {}          # query accession -> organism name
		self._gff_cache: Dict[str, List[dict]] = {}
		self._faa_cache: Dict[str, Dict[str, Tuple[str, str]]] = {}

	def extract(self, assembly: str, gff_path: str, faa_path: str,
				query: str, acceptable: Optional[set] = None) -> List[FlankingGene]:
		acceptable = acceptable or {query}
		genes = self._genes(assembly, gff_path)
		idx = next((i for i, g in enumerate(genes) if g["accession"] in acceptable), None)
		if idx is None:
			return []

		faa = self._faa(assembly, faa_path)
		contig = genes[idx]["contig"]
		qstrand = genes[idx]["strand"]
		lo, hi = max(0, idx - self.flank), min(len(genes), idx + self.flank + 1)

		neighborhood = []
		for j in range(lo, hi):
			g = genes[j]
			if g["contig"] != contig:
				continue
			acc = g["accession"]
			seq, organism = faa.get(acc, (None, ""))
			if seq is not None:
				self.sequences[acc] = seq
				if j == idx:
					self.query_sequences[acc] = seq
					self.species[query] = organism
			offset = j - idx
			if qstrand == "-":
				offset = -offset
			neighborhood.append(FlankingGene(
				accession=acc,
				strand=self._norm_strand(qstrand, g["strand"]),
				start=g["start"], end=g["end"],
				product=g["product"],
				offset=offset,
				query=query,
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
					genes.append({
						"contig": col[0], "start": int(col[3]), "end": int(col[4]),
						"strand": col[6], "accession": None, "product": "",
						"biotype": self._attr(attrs, "gene_biotype") or "",
					})
				elif feature == "CDS" and genes and genes[-1]["accession"] is None:
					genes[-1]["accession"] = self._attr(attrs, "Name")
					genes[-1]["product"] = self._attr(attrs, "product") or ""
		for g in genes:
			if not g["accession"]:
				g["accession"] = (g["biotype"] or "noProtein") + "*"
		self._gff_cache[assembly] = genes
		return genes

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

	@staticmethod
	def _attr(attributes: str, key: str) -> Optional[str]:
		m = re.search(r"(?:^|;){}=([^;]*)".format(re.escape(key)), attributes)
		return m.group(1) if m else None

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


class _FlaGsBase: #shared styling/colour helpers for the FlaGs visualizers
	GREY = "#d9d9d9"
	PSEUDO = ("#f2f2f3", "#000080")   # pseudogene: pale, navy outline
	RNA = ("#f2f2f2", "#008000")      # t/r/nc RNA: pale, green outline
	OTHER = ("#ffffff", "#bebebe")    # other non-coding: white, grey outline

	@staticmethod
	def _special_type(accession: str) -> Optional[str]:
		if not accession.endswith("*"):
			return None
		low = accession.lower()
		if low.startswith("pseudo") or low.startswith("ps"):
			return "pseudo"
		if "rna" in low:
			return "rna"
		return "other"

	def _style_for(self, accession: str, color: Dict[str, str]):
		special = self._special_type(accession)
		if special == "pseudo":
			return self.PSEUDO
		if special == "rna":
			return self.RNA
		if special == "other":
			return self.OTHER
		return color.get(accession, self.GREY), "#333"

	@staticmethod
	def _palette(n: int) -> List[str]:
		colors = []
		for i in range(n):
			r, g, b = colorsys.hsv_to_rgb(i / n if n else 0, 0.55, 0.85)
			colors.append("#{:02x}{:02x}{:02x}".format(
				int(r * 255), int(g * 255), int(b * 255)))
		return colors

	def _family_colors(self, families: List[List[str]]) -> Dict[str, str]:
		multi = [fam for fam in families if len(fam) > 1]
		palette = self._palette(len(multi))
		color = {}
		for fam, c in zip(multi, palette):
			for acc in fam:
				color[acc] = c
		for fam in families:
			if len(fam) == 1:
				color[fam[0]] = self.GREY
		return color

	def _family_numbers(self, families: List[List[str]]) -> Dict[str, int]:
		return family_numbers(families)

	@staticmethod
	def _escape(text: str) -> str:
		return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class OperonVisualizer(_FlaGsBase):
	def __init__(self, row_h: int = 22, gene_h: int = 14, bp_per_px: float = 24,
				 pad: int = 16, font: int = 13, show_numbers: bool = True):
		self.row_h = row_h; self.gene_h = gene_h
		self.bp_per_px = bp_per_px
		self.pad = pad; self.font = font
		self.show_numbers = show_numbers

	def render(self, neighborhoods: List["FlankingGene"],
			   families: List[List[str]],
			   species: Optional[Dict[str, str]] = None,
			   order: Optional[List[str]] = None) -> str:
		species = species or {}
		by_query: Dict[str, list] = {}
		for g in neighborhoods:
			by_query.setdefault(g.query, []).append(g)
		for q in by_query:
			by_query[q].sort(key=lambda x: x.start)

		color = self._family_colors(families)
		number = self._family_numbers(families)

		rows = [q for q in (order or list(by_query)) if q in by_query] or list(by_query)

		labels = {q: (q if not species.get(q) else "{}  {}".format(q, species[q]))
				  for q in rows}
		label_w = int(max((len(l) for l in labels.values()), default=0) * self.font * 0.55) + 6

		spans = {}
		max_left = max_right = 0
		for q in rows:
			genes = by_query[q]
			qg = next((g for g in genes if g.offset == 0), genes[len(genes) // 2])
			q_mid = (qg.start + qg.end) / 2
			left = (min(g.start for g in genes) - q_mid) / self.bp_per_px
			right = (max(g.end for g in genes) - q_mid) / self.bp_per_px
			spans[q] = q_mid
			max_left = min(max_left, left)
			max_right = max(max_right, right)

		track_x0 = self.pad + label_w
		center = track_x0 - max_left

		W = center + max_right + self.pad
		H = self.pad * 2 + len(rows) * self.row_h

		svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
			   'font-family="sans-serif" font-size="{}">'.format(W, H, self.font),
			   '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

		for i, q in enumerate(rows):
			y = self.pad + i * self.row_h + self.row_h / 2
			svg.append('<text x="{}" y="{}">{}</text>'.format(
				self.pad, y + self.font / 3, self._escape(labels[q])))
			q_mid = spans[q]
			reversed_row = self._row_reversed(by_query[q])
			for g in by_query[q]:
				if reversed_row:
					gx0 = center + (q_mid - g.end) / self.bp_per_px
					gx1 = center + (q_mid - g.start) / self.bp_per_px
				else:
					gx0 = center + (g.start - q_mid) / self.bp_per_px
					gx1 = center + (g.end - q_mid) / self.bp_per_px
				fill, outline = self._style_for(g.accession, color)
				stroke = "#000" if g.offset == 0 else outline
				sw = 2 if g.offset == 0 else 1
				num = number.get(g.accession) if self.show_numbers else None
				svg.append(self._block_arrow(gx0, gx1, y, g.strand, fill, stroke, sw, num))

		svg.append('</svg>')
		return "\n".join(svg)

	@staticmethod
	def _row_reversed(genes: List["FlankingGene"]) -> bool:
		q = next((g for g in genes if g.offset == 0), None)
		if q is None:
			return False
		pos = [g for g in genes if g.offset > 0]
		if not pos:
			neg = [g for g in genes if g.offset < 0]
			if not neg:
				return False
			ref = max(neg, key=lambda g: g.offset)
			return ref.start > q.start
		ref = max(pos, key=lambda g: g.offset)
		return ref.start < q.start

	def _block_arrow(self, x0, x1, cy, strand, fill, stroke, sw, number=None) -> str:
		h = self.gene_h
		length = max(x1 - x0, 6)              # floor so very short genes stay visible
		head = min(h, length * 0.5)
		top, bot = cy - h / 2, cy + h / 2
		if strand == "-":
			body_l = x0 + head
			pts = [(x0, cy), (body_l, top), (x0 + length, top),
				   (x0 + length, bot), (body_l, bot)]
		else:
			body_r = x0 + length - head
			pts = [(x0, top), (body_r, top), (x0 + length, cy),
				   (body_r, bot), (x0, bot)]
		points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
		out = ['<polygon points="{}" fill="{}" stroke="{}" stroke-width="{}"/>'.format(
			points, fill, stroke, sw)]
		if number is not None:
			out.append('<text x="{:.1f}" y="{:.1f}" font-size="{}" fill="black" '
					   'text-anchor="middle" dominant-baseline="central">{}</text>'.format(
						   x0 + length / 2, cy, self.font - 4, number))
		return "".join(out)


class ReportWriter: 
	def __init__(self, neighborhoods, families, species,
				 queries, protein_to_assemblies, matched):
		self.neighborhoods = neighborhoods
		self.families = families
		self.species = species
		self.queries = queries
		self.protein_to_assemblies = protein_to_assemblies
		self.matched = matched
		self.fam_of = family_numbers(families)
		self.by_query = {}
		for g in neighborhoods:
			self.by_query.setdefault(g.query, []).append(g)

	def write_all(self, out_path):
		self.operon_tsv(out_path("_operon.tsv"))
		self.clusters_tsv(out_path("_clusters.tsv"))
		self.species_info(out_path("_speciesInfo.txt"))
		self.query_status(out_path("_QueryStatus.txt"))
		self.flankgene_report(out_path("_flankgene_Report.log"))
		return self.accession_issues(out_path("_accessionIssues.txt"))

	def operon_tsv(self, path):
		with open(path, "w") as out:
			out.write("#query\tspecies\tfamily\tstrand\toffset\tstart\tend\taccession\tproduct\n")
			for query in self.by_query:
				sp = self.species.get(query, "")
				for g in sorted(self.by_query[query], key=lambda x: x.offset):
					out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
						query, sp, self.fam_of.get(g.accession, "-"), g.strand,
						g.offset, g.start, g.end, g.accession, g.product))

	def clusters_tsv(self, path):
		with open(path, "w") as out:
			out.write("#family\tsize\tmembers\n")
			for i, fam in enumerate(self.families, 1):
				out.write("{}\t{}\t{}\n".format(i, len(fam), ",".join(fam)))

	def species_info(self, path):
		with open(path, "w") as out:
			out.write("#query\tspecies\n")
			for query in sorted(self.species):
				out.write("{}\t{}\n".format(query, self.species[query]))

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

VERSION = "2.0.0"

def main():
	usage = ''' Description: Identify flanking genes and cluster them based on similarity. Requirement= Python3, BioPython, pyhmmer, requests. '''
	parser = argparse.ArgumentParser(description=usage)
	parser.add_argument("-i", "--input_list", required=True, help=" Protein Accession eg. WP_047256880.1, optionally tab-separated with an assembly Identifier eg. GCF_000001765.3, one per line. ")
	parser.add_argument("-u", "--user_email", required=True, help=" User Email Address (required by NCBI Entrez). ")
	parser.add_argument("-api", "--api_key", help=" NCBI API Key. ")
	parser.add_argument("-g", "--gene", type=int, default=4, help=" Number of flanking genes up/downstream. Default = 4 ")
	parser.add_argument("-m", "--max_assemblies", type=int, default=5, help=" Max assemblies per protein. Default = 5 ")
	parser.add_argument("-e", "--ethreshold", type=float, default=1e-3, help=" Jackhmmer inclusion E-value threshold for clustering flanking genes. Default = 1e-3 ")
	parser.add_argument("-n", "--number", type=int, default=3, help=" Number of jackhmmer iterations for clustering. Default = 3 ")
	parser.add_argument("-c", "--cpu", type=int, help=" Max parallel CPU workers (default: auto-detect). ")
	parser.add_argument("-ts", "--tshape", type=int, default=20, help=" Size of the flanking-gene triangles in the tree view. Default = 20 ")
	parser.add_argument("-tf", "--tfontsize", type=int, default=13, help=" Font size inside the tree-view triangles. Default = 13 ")
	parser.add_argument("-tmp", "--temporary", default="./genomes", help=" Temporary directory for downloaded assemblies; deleted at the end. Default = ./genomes ")
	parser.add_argument("-O", "--output", default="output", help=" Directory for result files; its name is also the file prefix. Default = output ")
	parser.add_argument("--tree", action="store_true", help=" Also build a phylogenetic tree with aligned neighbourhood triangles (<dir>_tree.svg). Does not affect the neighbours output. ")
	parser.add_argument("--tree_order", action="store_true", help=" Order the neighbours output by tree leaf order (implies --tree). ")
	parser.add_argument("--domains", action="store_true", help=" Scan flanking proteins for domains and write <dir>_domains.svg (requires --hmmdb). ")
	parser.add_argument("--hmmdb", help=" HMM database file (e.g. Pfam) for domain scanning. ")
	parser.add_argument("--clans", help=" Pfam-A.clans.tsv(.gz): colour domains by clan instead of family. ")
	parser.add_argument("-k", "--keep", action="store_true", help=" Keep the downloaded assemblies instead of deleting the temporary directory at the end. ")
	parser.add_argument("--use_local", metavar="DIR", help=" Directory of local .gff/.faa genome files to search before falling back to NCBI. Files may be gzipped; a genome is a .gff and .faa sharing a basename. ")
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

	timings = {}
	t0 = time.perf_counter()

	proteins_assembly, proteins_only = AccessionListReader(args.input_list).read()
	timings["1_read_input"] = time.perf_counter() - t0; t0 = time.perf_counter()
	all_queries = [p for p, _ in proteins_assembly] + list(proteins_only)
	if args.verbose:
		print(">> read {} queries ({} paired, {} protein-only)".format(
			len(all_queries), len(proteins_assembly), len(proteins_only)), flush=True)

	local = LocalGenomeResolver(args.use_local) if args.use_local else None
	local_files = {}                     # source id -> (gff_path, faa_path)
	protein_to_assemblies = {}           # protein -> [source ids]
	local_acceptable = {}                # protein -> {source id: {accession}}
	pending_pairs = []                   # (protein, assembly) still needing NCBI
	pending_only = []                    # bare proteins still needing NCBI

	if local:
		for protein, assembly in proteins_assembly:
			hit = local.resolve_pair(assembly)
			if hit:
				base, gff, faa = hit
				local_files[base] = (gff, faa)
				protein_to_assemblies.setdefault(protein, []).append(base)
				local_acceptable.setdefault(protein, {})[base] = {protein}
			else:
				pending_pairs.append([protein, assembly])
		for protein in proteins_only:
			hit = local.resolve_protein(protein)
			if hit:
				base, gff, faa = hit
				local_files[base] = (gff, faa)
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

	dl_workers = args.cpu if args.cpu else min(max(len(assemblies), 1), 16)
	dl = AssemblyDownloader(out_dir=args.temporary, workers=dl_workers)
	downloaded = dl.download_many(assemblies) if assemblies else {}
	downloaded.update(local_files)
	timings["3_download"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		dl_ok = sum(1 for g, f in downloaded.values() if g and f)
		print(">> genomes ready: {} of {} usable ({} downloaded, {} local)".format(
			dl_ok, len(downloaded), len(assemblies), len(local_files)), flush=True)

	extractor = NeighborhoodExtractor(flank=args.gene)
	all_neighborhoods = []
	for protein, asms in protein_to_assemblies.items():
		for asm in asms:
			gff, faa = downloaded.get(asm, (None, None))
			if gff and faa:
				if asm in local_files:
					acceptable = local_acceptable.get(protein, {}).get(asm)
				else:
					acceptable = mapper.accessions_in.get(protein, {}).get(asm)
				all_neighborhoods.extend(extractor.extract(asm, gff, faa, protein, acceptable))
	timings["4_extract_neighbors"] = time.perf_counter() - t0; t0 = time.perf_counter()
	matched = {n.query for n in all_neighborhoods}
	if args.verbose:
		print(">> extracted {} flanking-gene records; {} of {} queries matched".format(
			len(all_neighborhoods), len(matched), len(all_queries)), flush=True)

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

	clusterer = NeighborhoodClusterer(iterations=args.number, incE=args.ethreshold,
									  workers=args.cpu)
	families = clusterer.cluster(extractor.sequences)
	timings["5_clustering"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose:
		print(">> clustered {} flanking proteins into {} families".format(
			len(extractor.sequences), len(families)), flush=True)

	want_tree = args.tree or args.tree_order
	newick = ""
	leaf_order = None
	tree_mod = None
	if want_tree:
		import flags2_tree as tree_mod
		builder = tree_mod.TreeBuilder(threads=0)
		newick, _ = builder.build(extractor.query_sequences)
		if newick:
			leaf_order = tree_mod.ladderized_leaf_order(newick)
	timings["6_tree"] = time.perf_counter() - t0; t0 = time.perf_counter()
	if args.verbose and want_tree:
		print(">> tree: {}".format("built ({} leaves)".format(len(leaf_order))
			  if newick else "skipped (fewer than 3 query sequences)"), flush=True)

	os.makedirs(args.output, exist_ok=True)
	prefix = os.path.basename(os.path.normpath(args.output))
	def out_path(suffix):
		return os.path.join(args.output, prefix + suffix)
	order = leaf_order if args.tree_order else None
	operon = OperonVisualizer(show_numbers=True)
	with open(out_path("_neighbors.svg"), "w") as out:
		out.write(operon.render(all_neighborhoods, families, extractor.species, order))
	tree_written = False
	if newick:
		viz = tree_mod.NeighborhoodVisualizer(gene_h=args.tshape, font=args.tfontsize,
											  row_h=max(args.tshape + 4, args.tfontsize + 4))
		with open(out_path("_tree.svg"), "w") as out:
			out.write(viz.render(newick, all_neighborhoods, families, extractor.species))
		with open(out_path("_tree.nwk"), "w") as out:
			out.write(newick + "\n")
		tree_written = True
	timings["7_visualize"] = time.perf_counter() - t0
	domains_written = False
	if args.domains:
		if not args.hmmdb:
			print("--domains requires --hmmdb; skipping domain scan.")
		else:
			try:
				import flags2_domains as dom_mod
				scanner = dom_mod.DomainScanner(args.hmmdb, cpus=args.cpu or 0)
				domains = scanner.scan(extractor.sequences)
				clan_map = dom_mod.DomainScanner.load_clans(args.clans) if args.clans else None
				dom_viz = dom_mod.DomainVisualizer()
				with open(out_path("_domains.svg"), "w") as out:
					out.write(dom_viz.render(all_neighborhoods, families, domains,
											 extractor.species, order, clan_map))
				domains_written = True
			except Exception as e:
				print("   domain scan skipped: could not read/scan HMM database ({}).".format(e))
	timings["8_domains"] = time.perf_counter() - t0

	reporter = ReportWriter(all_neighborhoods, families, extractor.species,
							all_queries, protein_to_assemblies, matched)
	n_issues = reporter.write_all(out_path)
	if args.verbose:
		print(">> wrote data tables and reports ({} queries with issues)".format(n_issues),
			  flush=True)

	print("\n{} flanking proteins -> {} families".format(
		len(extractor.sequences), len(families)))
	print("\noutputs in {}/".format(args.output))
	print("  {}_neighbors.svg".format(prefix))
	if tree_written:
		print("  {}_tree.svg / {}_tree.nwk".format(prefix, prefix))
	elif want_tree:
		print("  (tree skipped: fewer than 3 query sequences)")
	if domains_written:
		print("  {}_domains.svg".format(prefix))
	for suffix in ("_operon.tsv", "_clusters.tsv", "_speciesInfo.txt",
				   "_QueryStatus.txt", "_flankgene_Report.log", "_accessionIssues.txt"):
		print("  {}{}".format(prefix, suffix))

	if args.verbose:
		print("\n--- timing (seconds) ---")
		for stage in sorted(timings):
			print("  {:24s} {:8.2f}".format(stage, timings[stage]))
		print("  {:24s} {:8.2f}".format("TOTAL", sum(timings.values())))


if __name__ == '__main__':
	try:
		main()
	except FileNotFoundError as e:
		sys.exit("Error: file not found - {}".format(e))
	except KeyboardInterrupt:
		sys.exit("\nInterrupted.")
	except Exception as e:
		sys.exit("Error: {}".format(e))