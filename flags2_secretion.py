import gzip
import os
import shutil
import subprocess
from typing import Dict, List, NamedTuple, Optional, Tuple


class SecretionHit(NamedTuple):
	assembly: str
	contig: str
	start: int          
	end: int
	type: str           
	probability: float  
	columns: Tuple[str, ...]   
	values: Tuple[str, ...]    


class SismisScanner: 
	CONTIG_COLS = ("sequence_id", "seq_id", "contig", "scaffold", "record_id")

	def __init__(self, out_dir: str):
		self.out_dir = out_dir
		os.makedirs(out_dir, exist_ok=True)
		self._cache: Dict[str, List[SecretionHit]] = {}   

	def scan_assembly(self, assembly: str, genome_path: str) -> List[SecretionHit]:
		if assembly in self._cache:
			return self._cache[assembly]

		asm_dir = os.path.join(self.out_dir, assembly)
		os.makedirs(asm_dir, exist_ok=True)
		fasta = self._decompressed(genome_path, asm_dir)

		sismis_out = os.path.join(asm_dir, "sismis")
		cmd = ["sismis", "run", "-g", fasta, "-o", sismis_out]
		try:
			proc = subprocess.run(cmd, capture_output=True, text=True)
		except FileNotFoundError:
			raise FileNotFoundError("sismis not found on PATH; install with 'pip install sismis'.")
		if proc.returncode != 0:
			raise RuntimeError("sismis failed on {}: {}".format(
				assembly, proc.stderr[-500:] if proc.stderr else "(no stderr)"))

		hits = self._parse_clusters(sismis_out, assembly)
		self._cache[assembly] = hits
		return hits

	@staticmethod
	def _decompressed(genome_path: str, asm_dir: str) -> str:
		if not genome_path.endswith(".gz"):
			return genome_path
		local = os.path.join(asm_dir, "genome.fna")
		with gzip.open(genome_path, "rt", encoding="utf-8", errors="replace") as fin, \
			 open(local, "w") as fout:
			shutil.copyfileobj(fin, fout)
		return local

	def _parse_clusters(self, sismis_out: str, assembly: str) -> List[SecretionHit]:
		clusters_tsv = self._find_clusters_tsv(sismis_out)
		if clusters_tsv is None:
			return []
		hits = []
		with open(clusters_tsv) as fh:
			header = fh.readline().rstrip("\n").split("\t")
			lower = [h.lower() for h in header]
			contig_i = next((lower.index(c) for c in self.CONTIG_COLS if c in lower), None)
			start_i = lower.index("start") if "start" in lower else None
			end_i = lower.index("end") if "end" in lower else None
			type_i = lower.index("type") if "type" in lower else None
			prob_i = lower.index("max_p") if "max_p" in lower else None
			missing = [name for name, i in (("contig", contig_i), ("start", start_i),
											 ("end", end_i), ("type", type_i),
											 ("max_p", prob_i)) if i is None]
			if missing:
				raise RuntimeError(
					"{} is missing column(s) {} needed to use its output; "
					"columns seen: {}".format(clusters_tsv, missing, header))
			for line in fh:
				if not line.strip():
					continue
				values = line.rstrip("\n").split("\t")
				hits.append(SecretionHit(
					assembly=assembly, contig=values[contig_i],
					start=int(values[start_i]), end=int(values[end_i]),
					type=values[type_i], probability=float(values[prob_i]),
					columns=tuple(header), values=tuple(values)))
		return hits

	@staticmethod
	def _find_clusters_tsv(sismis_out: str) -> Optional[str]:
		if not os.path.isdir(sismis_out):
			return None
		for name in os.listdir(sismis_out):
			if name.endswith(".clusters.tsv"):
				return os.path.join(sismis_out, name)
		return None


def match_rows(hits: List[SecretionHit], rows: Dict[str, Tuple[str, str, int, int]]
			   ) -> Dict[int, List[str]]:
	matches: Dict[int, List[str]] = {}
	for i, h in enumerate(hits):
		for row_id, (assembly, contig, lo, hi) in rows.items():
			if h.assembly == assembly and h.contig == contig and h.start <= hi and h.end >= lo:
				matches.setdefault(i, []).append(row_id)
	return matches


def write_report(hits: List[SecretionHit], matches: Dict[int, List[str]], path: str):
	with open(path, "w") as out:
		if not hits:
			out.write("#assembly\t(no secretion systems predicted)\n")
			return
		header = hits[0].columns
		out.write("#assembly\t{}\toverlapping_rows\n".format("\t".join(header)))
		for i, h in enumerate(hits):
			values = h.values
			if len(values) != len(header):
				values = (values + ("",) * len(header))[:len(header)]
			rows = ",".join(matches.get(i, [])) or "-"
			out.write("{}\t{}\t{}\n".format(h.assembly, "\t".join(values), rows))


def write_diagnostics(statuses: Dict[str, str], path: str):
	with open(path, "w") as out:
		out.write("#assembly\tstatus\n")
		for assembly in sorted(statuses):
			out.write("{}\t{}\n".format(assembly, statuses[assembly]))