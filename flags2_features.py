import os
import re
import tempfile
import threading
from typing import Dict, List, Tuple

FeatureRegion = Tuple[str, int, int]

_SUBMIT_LOCK = threading.Lock()
_WARMUP_LOCK = threading.Lock()
_warmed_up = False


def warm_up() -> bool:

	global _warmed_up
	with _WARMUP_LOCK:
		if _warmed_up:
			return True
		try:
			import biolib  
			from biolib.biolib_api_client import BiolibApiClient
		except ImportError:
			return False
		try:
			BiolibApiClient.get()
		except Exception:
			pass   
		_warmed_up = True
		return True


def _runs(topology: str, wanted: str, kind: str) -> List[FeatureRegion]:
	regions = []
	start = None
	for i, code in enumerate(topology):
		if code in wanted:
			if start is None:
				start = i
		elif start is not None:
			regions.append((kind, start + 1, i))
			start = None
	if start is not None:
		regions.append((kind, start + 1, len(topology)))
	return regions


def parse_deeptmhmm_3line(text: str, want_signal: bool = True) -> Dict[str, List[FeatureRegion]]:
	features: Dict[str, List[FeatureRegion]] = {}
	lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
	i = 0
	while i + 2 < len(lines) + 1 and i < len(lines):
		if not lines[i].startswith(">"):
			i += 1
			continue
		name = lines[i][1:].split("|")[0].split()[0].strip()
		topology = lines[i + 2] if i + 2 < len(lines) else ""
		regions = _runs(topology, "MB", "tm")
		if want_signal:
			regions += _runs(topology, "S", "signal")
		if regions:
			features[name] = regions
		i += 3
	return features


def parse_signalp_regions(text: str) -> Dict[str, List[FeatureRegion]]:
	features: Dict[str, List[FeatureRegion]] = {}
	lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
	i = 0
	while i < len(lines):
		if not lines[i].startswith(">"):
			i += 1
			continue
		name = lines[i][1:].split("|")[0].split()[0].strip()
		labels = lines[i + 2] if i + 2 < len(lines) else ""
		regions = _runs(labels, "STLP", "signal")
		if regions:
			features[name] = [regions[0]]
		i += 3
	return features


class _BioLibScanner:

	def __init__(self):
		import biolib 
		self._biolib = biolib

	def _write_fasta(self, sequences: Dict[str, str], path: str):
		with open(path, "w") as out:
			for name, seq in sequences.items():
				out.write(">{}\n{}\n".format(name, seq))

	def _run_batched(self, app_slug, sequences, result_suffix, parse,
					 args_template="--fasta {fasta}", batch_size=25):
		items = list(sequences.items())
		merged: Dict[str, List[FeatureRegion]] = {}
		for i in range(0, len(items), batch_size):
			chunk = dict(items[i:i + batch_size])
			text = self._run(app_slug, chunk, result_suffix, args_template=args_template)
			for name, regions in parse(text).items():
				merged.setdefault(name, []).extend(regions)
		return merged

	def _run(self, app_slug: str, sequences: Dict[str, str], result_suffix: str,
			 args_template: str = "--fasta {fasta}") -> str:
		tmp = tempfile.mkdtemp(prefix="flags2_biolib_")
		fasta_name = "query.fasta"
		self._write_fasta(sequences, os.path.join(tmp, fasta_name))
		with _SUBMIT_LOCK:
			prev_cwd = os.getcwd()
			os.chdir(tmp)
			try:
				app = self._biolib.load(app_slug)
				job = app.cli(args=args_template.format(fasta=fasta_name))
			finally:
				os.chdir(prev_cwd)
		if hasattr(job, "wait"):
			job.wait()  

		out_dir = os.path.join(tmp, "out")
		saved = []
		problems = []
		try:
			job.save_files(out_dir)
			for root, _, names in os.walk(out_dir):
				for name in names:
					full = os.path.join(root, name)
					saved.append(full)
					if full.replace("\\", "/").endswith(result_suffix):
						with open(full) as fh:
							return fh.read()
		except Exception as e:
			problems.append("save_files: {!r}".format(e))
		try:
			listed = job.list_output_files()
			paths = [f if isinstance(f, str) else getattr(f, "path", str(f)) for f in listed]
			match = next((p for p in paths if p.replace("\\", "/").endswith(result_suffix)), None)
			if match is not None:
				out_file = job.get_output_file(match)
				data = (out_file.get_data() if hasattr(out_file, "get_data")
						else out_file.get_file_handle().read())
				return data.decode() if isinstance(data, bytes) else data
			saved = saved or paths
		except Exception as e:
			problems.append("list_output_files: {!r}".format(e))
		stdout = ""
		try:
			stdout = job.get_stdout().decode(errors="replace")[-500:]
		except Exception as e:
			problems.append("get_stdout: {!r}".format(e))
		raise FileNotFoundError(
			"{} produced no {}.\n  files seen: {}\n  retrieval errors: {}\n"
			"  stdout tail: {}".format(
				app_slug, result_suffix, saved or "(none)",
				"; ".join(problems) or "(none)", stdout or "(none)"))


class TMScanner(_BioLibScanner):
	APP = "DTU/DeepTMHMM"
	RESULT = "predicted_topologies.3line"

	def scan(self, sequences: Dict[str, str], want_signal: bool = True
			 ) -> Dict[str, List[FeatureRegion]]:
		if not sequences:
			return {}
		parse = lambda text: parse_deeptmhmm_3line(text, want_signal=want_signal)
		return self._run_batched(self.APP, sequences, self.RESULT, parse)


class SignalPScanner(_BioLibScanner):

	APP = "DTU/SignalP-6"

	RESULT = "prediction_results.txt"
	ARGS = "--fastafile {fasta} --output_dir output --organism other --format txt --mode fast"

	def scan(self, sequences: Dict[str, str]) -> Dict[str, List[FeatureRegion]]:
		if not sequences:
			return {}
		return self._run_batched(self.APP, sequences, self.RESULT,
								 self._parse_prediction_results, args_template=self.ARGS)

	@staticmethod
	def _parse_prediction_results(text: str) -> Dict[str, List[FeatureRegion]]:
		features: Dict[str, List[FeatureRegion]] = {}
		cs_re = re.compile(r"CS pos:\s*(\d+)")
		for line in text.splitlines():
			if line.startswith("#") or not line.strip():
				continue
			col = line.split("\t")
			if len(col) < 2:
				continue
			name = col[0].split()[0].strip()
			pred = col[1].strip().upper()
			if pred in ("", "OTHER"):
				continue   
			end = None
			m = cs_re.search(line)
			if m:
				end = int(m.group(1))
			if end is None:
				continue
			features.setdefault(name, []).append(("signal", 1, end))
		return features