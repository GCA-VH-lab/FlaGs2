import colorsys
import gzip
from typing import Dict, List, Optional, NamedTuple

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from pyhmmer.plan7 import HMMFile

from FlaGs2 import _FlaGsBase, OperonVisualizer, FlankingGene


class DomainHit(NamedTuple):
	protein: str     # prot accession
	name: str        # dom name
	start: int       
	end: int
	evalue: float    # per-domain E-value


class DomainScanner: 

	def __init__(self, hmm_db: str, evalue: float = 1e-10, cpus: int = 0):
		self.hmm_db = hmm_db
		self.evalue = evalue
		self.cpus = cpus            # 0 -> all
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
					continue                       # family belongs to no clan
				clan = clan_name or clan_id
				if pfam_id:
					mapping[pfam_id] = clan
				if family_name:
					mapping[family_name] = clan
		return mapping


class DomainVisualizer(_FlaGsBase):

	LABEL_STEP = 12   

	def __init__(self, row_h: int = 22, gene_h: int = 8, bp_per_px: float = 18,
				 pad: int = 16, font: int = 13, domain_h: int = 6):
		self.row_h = row_h; self.gene_h = gene_h
		self.bp_per_px = bp_per_px
		self.pad = pad; self.font = font; self.domain_h = domain_h

	def render(self, neighborhoods: List["FlankingGene"],
			   families: List[List[str]],
			   domains: Dict[str, List[DomainHit]],
			   species: Optional[Dict[str, str]] = None,
			   order: Optional[List[str]] = None,
			   clans: Optional[Dict[str, str]] = None) -> str:
		species = species or {}
		clans = clans or {}
		by_query: Dict[str, list] = {}
		for g in neighborhoods:
			by_query.setdefault(g.query, []).append(g)
		for q in by_query:
			by_query[q].sort(key=lambda x: x.start)

		gene_color = self._family_colors(families)
		rows = [q for q in (order or list(by_query)) if q in by_query] or list(by_query)

		group_index: Dict[str, int] = {}   
		group_label: Dict[str, str] = {}    
		domain_group: Dict[str, str] = {}   
		for q in rows:
			for g in by_query[q]:
				for d in domains.get(g.accession, []):
					clan = clans.get(d.name)
					group = clan if clan else d.name
					if group not in group_index:
						group_index[group] = len(group_index)
						group_label[group] = group        
					domain_group[d.name] = group
		domain_index = {name: group_index[grp] for name, grp in domain_group.items()}
		domain_color = {name: self._contrast_color(group_index[grp])
						for name, grp in domain_group.items()}

		labels = {q: (q if not species.get(q) else "{}  {}".format(q, species[q]))
				  for q in rows}
		label_w = int(max((len(l) for l in labels.values()), default=0) * self.font * 0.55) + 6

		spans = {}
		max_left = max_right = 0
		for q in rows:
			genes = by_query[q]
			qg = next((g for g in genes if g.offset == 0), genes[len(genes) // 2])
			spans[q] = (qg.start + qg.end) / 2
			max_left = min(max_left, (min(g.start for g in genes) - spans[q]) / self.bp_per_px)
			max_right = max(max_right, (max(g.end for g in genes) - spans[q]) / self.bp_per_px)

		center = self.pad + label_w - max_left
		W = center + max_right + self.pad

		legend_items = sorted(
			((group_index[grp], group_label[grp], self._contrast_color(group_index[grp]))
			 for grp in group_index), key=lambda t: t[0])
		if legend_items:
			longest = max(len("{}. {}".format(num, lbl)) for num, lbl, _ in legend_items)
			col_w = int(longest * self.font * 0.55) + 26   
			avail = W - 2 * self.pad
			cols = max(1, min(len(legend_items), avail // col_w))
			legend_rows = -(-len(legend_items) // cols)    
			legend_h = (legend_rows + 1) * 16 + 8
		else:
			col_w = cols = legend_rows = 0
			legend_h = 0
		H = self.pad * 2 + len(rows) * self.row_h + legend_h

		svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
			   'font-family="sans-serif" font-size="{}">'.format(W, H, self.font),
			   '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

		clip_n = 0
		for i, q in enumerate(rows):
			y = self.pad + i * self.row_h + self.row_h / 2
			svg.append('<text x="{}" y="{}">{}</text>'.format(
				self.pad, y + self.font / 3, self._escape(labels[q])))
			q_mid = spans[q]
			reversed_row = OperonVisualizer._row_reversed(by_query[q])
			for g in by_query[q]:
				if reversed_row:
					gx0 = center + (q_mid - g.end) / self.bp_per_px
					gx1 = center + (q_mid - g.start) / self.bp_per_px
				else:
					gx0 = center + (g.start - q_mid) / self.bp_per_px
					gx1 = center + (g.end - q_mid) / self.bp_per_px
				drawn_strand = g.strand
				_, outline = self._style_for(g.accession, gene_color)
				stroke = "#000" if g.offset == 0 else outline
				sw = 2 if g.offset == 0 else 1
				wedges, labels_svg = self._domains_on_gene(
					g, gx0, gx1, y, drawn_strand, domains.get(g.accession, []),
					domain_index, domain_color)
				if wedges:
					clip_id = "clip{}".format(clip_n)
					clip_n += 1
					clip_pts = " ".join("{:.1f},{:.1f}".format(px, py)
										for px, py in self._arrow_points(gx0, gx1, y, drawn_strand))
					svg.append('<clipPath id="{}"><polygon points="{}"/></clipPath>'.format(
						clip_id, clip_pts))
					svg.append('<g clip-path="url(#{})">{}</g>'.format(clip_id, wedges))
				svg.append(self._gene_outline(gx0, gx1, y, drawn_strand, stroke, sw))
				svg.append(labels_svg)

		if legend_items:
			svg.append(self._domain_legend(legend_items,
										   self.pad, self.pad + len(rows) * self.row_h + 14,
										   cols, col_w))
		svg.append('</svg>')
		return "\n".join(svg)

	def _arrow_points(self, x0, x1, cy, strand):
		h = self.gene_h
		length = max(x1 - x0, 6)
		head = min(h, length * 0.5)
		top, bot = cy - h / 2, cy + h / 2
		if strand == "-":
			body_l = x0 + head
			return [(x0, cy), (body_l, top), (x0 + length, top),
					(x0 + length, bot), (body_l, bot)]
		body_r = x0 + length - head
		return [(x0, top), (body_r, top), (x0 + length, cy),
				(body_r, bot), (x0, bot)]

	def _gene_outline(self, x0, x1, cy, strand, stroke, sw) -> str:
		points = " ".join("{:.1f},{:.1f}".format(px, py)
						  for px, py in self._arrow_points(x0, x1, cy, strand))
		return ('<polygon points="{}" fill="rgba(0,0,0,0)" stroke="{}" '
				'stroke-width="{}"/>'.format(points, stroke, sw))

	def _domains_on_gene(self, gene, gx0, gx1, cy, drawn_strand, hits,
						 domain_index, domain_color):
		if not hits:
			return "", ""
		minus = (drawn_strand == "-")
		direction = -1 if minus else 1
		u = self.domain_h / 4.0
		step = self.LABEL_STEP

		prot_len = max((gene.end - gene.start) // 3, 1)
		span = gx1 - gx0

		def res_to_x(res):
			frac = min(max(res / prot_len, 0.0), 1.0)
			return gx1 - frac * span if minus else gx0 + frac * span

		wedges, labels = [], []
		prev_labels = []
		for d in hits:
			s, e = res_to_x(d.start), res_to_x(d.end)
			color = domain_color[d.name]
			pts = [(s, cy + 2 * u), (s, cy + 1 * u), (e, cy - 2 * u),
				   (e, cy + 2 * u), (s, cy + 2 * u)]
			points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
			wedges.append('<polygon points="{}" fill="{}"/>'.format(points, color))
			num = domain_index[d.name]
			label_x = (s + e) / 2
			if any(abs(label_x - p) < step for p in prev_labels):
				label_x = prev_labels[-1] + step * direction
			prev_labels.append(label_x)
			prev_labels.sort(reverse=minus)
			labels.append('<text x="{:.1f}" y="{:.1f}" font-size="{}" fill="{}" '
						  'text-anchor="middle">{}</text>'.format(
							  label_x, cy - 2 * u - 3, self.font - 4, color, num))
		return "".join(wedges), "".join(labels)

	def _domain_legend(self, legend_items, x, y, cols, col_w) -> str:
		parts = ['<text x="{}" y="{}" font-weight="bold">Domains</text>'.format(x, y)]
		for idx, (num, label, color) in enumerate(legend_items):
			col = idx % cols
			row = idx // cols
			ex = x + col * col_w
			ey = y + 16 + row * 16
			parts.append('<rect x="{}" y="{}" width="12" height="10" fill="{}"/>'.format(
				ex + 2, ey - 8, color))
			parts.append('<text x="{}" y="{}">{}. {}</text>'.format(
				ex + 20, ey, num, self._escape(label)))
		return "".join(parts)

	@staticmethod
	def _contrast_color(n: int) -> str:
		r, s = 3, 7
		c = (n % r / r) + n // r / (s + n // (r * s) * 2)
		if c < 0:
			c = 1 - c
		if c > 1:
			c = c - 1
		red, green, blue = colorsys.hls_to_rgb(c, 0.5, 0.5)
		return "#{:02x}{:02x}{:02x}".format(int(red * 255), int(green * 255), int(blue * 255))