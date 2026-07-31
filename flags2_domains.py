
import gzip
from typing import Dict, List, NamedTuple

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import HMMFile


class DomainHit(NamedTuple):
	protein: str
	name: str 
	start: int 
	end: int
	evalue: float


class DomainScanner: 
	def __init__(self, hmm_db: str, evalue: float = 1e-10, cpus: int = 0):
		self.hmm_db = hmm_db
		self.evalue = evalue
		self.cpus = cpus  
		self.alphabet = Alphabet.amino()

	def scan(self, sequences: Dict[str, str]) -> Dict[str, List[DomainHit]]:
		if not sequences:
			return {}
		block = DigitalSequenceBlock(self.alphabet, [
			TextSequence(name=name.encode(), sequence=seq).digitize(self.alphabet)
			for name, seq in sequences.items()])
		with HMMFile(self.hmm_db) as handle:
			if handle.is_pressed:
				profiles = list(handle.optimized_profiles())
			else:
				profiles = list(handle)

		hits: Dict[str, List[DomainHit]] = {name: [] for name in sequences}
		for top in pyhmmer.hmmer.hmmsearch(profiles, block,
										   E=self.evalue, cpus=self.cpus):
			domain_name = self._decode(top.query.name) 
			for hit in top:
				protein = self._decode(hit.name)
				for dom in hit.domains:
					if not dom.included:
						continue
					al = dom.alignment
					hits[protein].append(DomainHit(
						protein=protein, name=domain_name,
						start=al.target_from, end=al.target_to, evalue=dom.i_evalue))
		return hits

	@staticmethod
	def _decode(value) -> str:
		return value.decode() if isinstance(value, (bytes, bytearray)) else value

	@staticmethod
	def load_clans(path: str) -> Dict[str, str]:
		mapping: Dict[str, str] = {}
		opener = gzip.open if path.endswith(".gz") else open
		with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
			for line in fh:
				cols = line.rstrip("\n").split("\t")
				if len(cols) < 4:
					continue
				pfam_id, clan_id, clan_name, family_name = cols[0], cols[1], cols[2], cols[3]
				if not clan_id:
					continue 
				clan = clan_name or clan_id
				if pfam_id:
					mapping[pfam_id] = clan
				if family_name:
					mapping[family_name] = clan
		return mapping