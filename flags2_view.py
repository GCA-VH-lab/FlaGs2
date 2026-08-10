import colorsys
from typing import Dict, List, Optional


def family_numbers(families, rna_accessions=None):
	rna_accessions = rna_accessions or set()
	number = {}
	prot_n, rna_n = 0, 0
	for fam in families:
		if len(fam) < 2:
			continue
		if fam[0] in rna_accessions:
			rna_n += 1
			label = "R{}".format(rna_n)
		else:
			prot_n += 1
			label = str(prot_n)
		for acc in fam:
			number[acc] = label
	return number


class _FlaGsBase: 
	GREY = "#d9d9d9"
	PSEUDO = ("#f2f2f3", "#000080")   
	RNA = ("#f2f2f2", "#008000")      
	OTHER = ("#ffffff", "#bebebe")    

	FONT_FAMILY = "Arial, 'Liberation Sans', Helvetica, sans-serif"

	_CHAR_WIDTH = {
		' ':278,'!':278,'"':355,'#':556,'$':556,'%':889,'&':667,"'":191,
		'(':333,')':333,'*':389,'+':584,',':278,'-':333,'.':278,'/':278,
		'0':556,'1':556,'2':556,'3':556,'4':556,'5':556,'6':556,'7':556,'8':556,'9':556,
		':':278,';':278,'<':584,'=':584,'>':584,'?':556,'@':1015,
		'A':667,'B':667,'C':722,'D':722,'E':667,'F':611,'G':778,'H':722,'I':278,
		'J':500,'K':667,'L':556,'M':833,'N':722,'O':778,'P':667,'Q':778,'R':722,
		'S':667,'T':611,'U':722,'V':667,'W':944,'X':667,'Y':667,'Z':611,
		'[':278,'\\':278,']':278,'^':469,'_':556,'`':333,
		'a':556,'b':556,'c':500,'d':556,'e':556,'f':278,'g':556,'h':556,'i':222,
		'j':222,'k':500,'l':222,'m':833,'n':556,'o':556,'p':556,'q':556,'r':333,
		's':500,'t':278,'u':556,'v':500,'w':722,'x':500,'y':500,'z':500,
		'{':334,'|':260,'}':334,'~':584,
	}
	_CHAR_WIDTH_DEFAULT = 556

	@classmethod
	def _text_width(cls, text: str, font_size: float) -> float:
		units = sum(cls._CHAR_WIDTH.get(ch, cls._CHAR_WIDTH_DEFAULT) for ch in text)
		return units / 1000.0 * font_size

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

	def _style_for(self, accession: str, color: Dict[str, str], is_rna: bool = False):
		special = self._special_type(accession)
		if special == "pseudo":
			return self.PSEUDO
		if special == "rna":
			return self.RNA
		if special == "other":
			return self.OTHER
		if is_rna:
			return color.get(accession, self.RNA[0]), self.RNA[1]
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

	def _family_numbers(self, families: List[List[str]], rna_accessions=None):
		return family_numbers(families, rna_accessions)

	@staticmethod
	def _escape(text: str) -> str:
		return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class OperonView(_FlaGsBase):
	LABEL_STEP = 12   
	SECRETION_BAND_OPACITY = 0.28   
	SECRETION_BAND_PAD = 4          

	def __init__(self, mode: str = "families", row_h: int = 26, gene_h: int = 8,
				 bp_per_px: float = 10.4, pad: int = 16, font: int = 13,
				 domain_h: int = 6):
		self.mode = mode
		self.row_h = row_h
		self.gene_h = gene_h
		self.bp_per_px = bp_per_px
		self.pad = pad
		self.font = font
		self.domain_h = domain_h

	def render(self, neighborhoods, families, species=None, order=None,
			   domains=None, clans=None, features=None, labels=None, secretion=None):
		species = species or {}
		features = features or {}
		row_labels_map = labels or {}
		by_query = {}
		for g in neighborhoods:
			by_query.setdefault(g.query, []).append(g)
		for q in by_query:
			by_query[q].sort(key=lambda x: x.start)

		rows = [q for q in (order or list(by_query)) if q in by_query] or list(by_query)
		gene_color = self._family_colors(families)
		rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
		if self.mode == "domains":
			overlay = self._domain_overlay(rows, by_query, domains or {}, clans or {})
		elif self.mode == "secretion":
			overlay = self._secretion_overlay(rows, by_query, secretion or [])
		else:
			overlay = self._family_overlay(families, rna_accessions)
		def _row_label(q):
			if q in row_labels_map:
				return row_labels_map[q]
			return q if not species.get(q) else "{}  {}".format(q, species[q])
		labels_out = {q: _row_label(q) for q in rows}
		label_w = int(max((self._text_width(l, self.font) for l in labels_out.values()),
						   default=0)) + 2

		spans = {}
		row_reversed = {}
		max_left = max_right = 0
		for q in rows:
			genes = by_query[q]
			qg = next((g for g in genes if g.offset == 0), genes[len(genes) // 2])
			spans[q] = (qg.start + qg.end) / 2
			reversed_row = self._row_reversed(genes)
			row_reversed[q] = reversed_row
			if reversed_row:
				left_reach = (spans[q] - max(g.end for g in genes)) / self.bp_per_px
				right_reach = (spans[q] - min(g.start for g in genes)) / self.bp_per_px
			else:
				left_reach = (min(g.start for g in genes) - spans[q]) / self.bp_per_px
				right_reach = (max(g.end for g in genes) - spans[q]) / self.bp_per_px
			max_left = min(max_left, left_reach)
			max_right = max(max_right, right_reach)

		center = self.pad + label_w - max_left
		genes_right = center + max_right 
		W = int(genes_right + self.pad)

		legend_items, cols, col_w, legend_h = self._legend_layout(overlay, W)
		H = self.pad * 2 + len(rows) * self.row_h + legend_h

		svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
			   'font-family="{}" font-size="{}">'.format(W, H, self.FONT_FAMILY, self.font),
			   '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

		clip_n = 0
		for i, q in enumerate(rows):
			y = self.pad + i * self.row_h + self.row_h / 2
			q_mid = spans[q]
			reversed_row = row_reversed[q]

			svg.append('<text x="{}" y="{}">{}</text>'.format(
				self.pad, y + self.font / 3, self._escape(labels_out[q])))
			row_labels = []

			if self.mode == "secretion":
				gxs = []
				for g in by_query[q]:
					gxs.append(self._x_for(g.start, center, q_mid, reversed_row))
					gxs.append(self._x_for(g.end, center, q_mid, reversed_row))
				row_lo, row_hi = (min(gxs), max(gxs)) if gxs else (0, 0)
				for h in overlay["row_hits"].get(q, []):
					band, mid = self._secretion_band(
						h, overlay["color"][h.type], center, q_mid, reversed_row,
						y, row_lo, row_hi)
					if band:
						svg.append(band)
						row_labels.append((mid, h.type, "#000"))

			for g in by_query[q]:
				if reversed_row:
					gx0 = center + (q_mid - g.end) / self.bp_per_px
					gx1 = center + (q_mid - g.start) / self.bp_per_px
					drawn = "+" if g.strand == "-" else "-"
				else:
					gx0 = center + (g.start - q_mid) / self.bp_per_px
					gx1 = center + (g.end - q_mid) / self.bp_per_px
					drawn = g.strand
				drawn = g.strand if not reversed_row else drawn

				fill, outline = self._gene_style(g, gene_color)
				stroke = "#000" if g.offset == 0 else outline
				sw = 2 if g.offset == 0 else 1

				if self.mode == "domains":
					wedges, wlabels = self._domains_on_gene(
						g, gx0, gx1, y, g.strand, (domains or {}).get(g.accession, []),
						overlay)
					feats = self._features_on_gene(
						g, gx0, gx1, y, g.strand, features.get(g.accession, []))
					overlay_svg = wedges + feats
					if overlay_svg:
						clip_id = "clip{}".format(clip_n); clip_n += 1
						clip_pts = " ".join("{:.1f},{:.1f}".format(px, py)
											for px, py in self._arrow_points(gx0, gx1, y, g.strand))
						svg.append('<clipPath id="{}"><polygon points="{}"/></clipPath>'.format(
							clip_id, clip_pts))
						svg.append('<g clip-path="url(#{})">{}</g>'.format(clip_id, overlay_svg))
					svg.append(self._gene_outline(gx0, gx1, y, g.strand, stroke, sw))
					row_labels.extend(wlabels)
				elif self.mode == "secretion":
					svg.append(self._gene_fill(gx0, gx1, y, g.strand, fill, stroke, sw))
				else:
					svg.append(self._gene_fill(gx0, gx1, y, g.strand, fill, stroke, sw))
					num = overlay["number"].get(g.accession)
					if num is not None:
						row_labels.append(((gx0 + gx1) / 2, num, "#000"))

			svg.append(self._place_labels(row_labels, y))

		if legend_items:
			svg.append(self._legend(legend_items, self.pad,
									self.pad + len(rows) * self.row_h + 14, cols, col_w,
									title="Secretion systems" if self.mode == "secretion" else "Domains"))
		svg.append('</svg>')
		return "\n".join(svg)

	def _family_overlay(self, families, rna_accessions):
		return {
			"color": self._family_colors(families),
			"number": self._family_numbers(families, rna_accessions),
			"legend": None,
		}

	def _domain_overlay(self, rows, by_query, domains, clans):
		group_index, group_label, domain_group = {}, {}, {}
		for q in rows:
			for g in by_query[q]:
				for d in domains.get(g.accession, []):
					group = clans.get(d.name) or d.name
					if group not in group_index:
						group_index[group] = len(group_index)
						group_label[group] = group
					domain_group[d.name] = group
		return {
			"number": {name: group_index[grp] for name, grp in domain_group.items()},
			"color": {name: self._contrast_color(group_index[grp])
					  for name, grp in domain_group.items()},
			"legend": sorted(
				((group_index[grp], group_label[grp], self._contrast_color(group_index[grp]))
				 for grp in group_index), key=lambda t: t[0]),
		}

	def _secretion_overlay(self, rows, by_query, hits):
		by_loc: Dict[tuple, list] = {}
		for h in hits:
			by_loc.setdefault((h.assembly, h.contig), []).append(h)

		gene_hits: Dict[tuple, list] = {}   
		row_hits: Dict[str, list] = {}      
		seen_types = []
		for q in rows:
			assembly = q.rsplit("|", 1)[-1]
			for g in by_query[q]:
				matches = [h for h in by_loc.get((assembly, g.contig), [])
						   if h.start <= g.end and h.end >= g.start]
				if matches:
					gene_hits[(q, g.offset)] = matches
					for h in matches:
						if h not in row_hits.setdefault(q, []):
							row_hits[q].append(h)
						if h.type not in seen_types:
							seen_types.append(h.type)
		seen_types.sort()
		color = {t: self._contrast_color(i) for i, t in enumerate(seen_types)}
		return {
			"gene_hits": gene_hits,
			"row_hits": row_hits,
			"color": color,
			"legend": [(i, t, color[t]) for i, t in enumerate(seen_types)],
		}

	def _x_for(self, coord, center, q_mid, reversed_row):
		if reversed_row:
			return center + (q_mid - coord) / self.bp_per_px
		return center + (coord - q_mid) / self.bp_per_px

	def _secretion_band(self, hit, color, center, q_mid, reversed_row, y,
						row_lo, row_hi):
		xa = self._x_for(hit.start, center, q_mid, reversed_row)
		xb = self._x_for(hit.end, center, q_mid, reversed_row)
		x0, x1 = (xa, xb) if xa <= xb else (xb, xa)
		x0 = max(x0, row_lo)
		x1 = min(x1, row_hi)
		if x1 <= x0:
			return "", None
		h = self.gene_h + self.SECRETION_BAND_PAD * 2
		rect = ('<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" '
				'fill="{}" fill-opacity="{}" stroke="{}" stroke-opacity="{}" '
				'stroke-width="1" rx="2"/>').format(
					x0, y - h / 2, x1 - x0, h, color, self.SECRETION_BAND_OPACITY,
					color, min(1.0, self.SECRETION_BAND_OPACITY * 2))
		return rect, (x0 + x1) / 2

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

	def _gene_outline(self, x0, x1, cy, strand, stroke, sw):
		points = " ".join("{:.1f},{:.1f}".format(px, py)
						  for px, py in self._arrow_points(x0, x1, cy, strand))
		return ('<polygon points="{}" fill="rgba(0,0,0,0)" stroke="{}" '
				'stroke-width="{}"/>'.format(points, stroke, sw))

	def _gene_fill(self, x0, x1, cy, strand, fill, stroke, sw):
		points = " ".join("{:.1f},{:.1f}".format(px, py)
						  for px, py in self._arrow_points(x0, x1, cy, strand))
		return '<polygon points="{}" fill="{}" stroke="{}" stroke-width="{}"/>'.format(
			points, fill, stroke, sw)

	def _domains_on_gene(self, gene, gx0, gx1, cy, drawn_strand, hits, overlay):
		if not hits:
			return "", []
		minus = (drawn_strand == "-")
		u = self.domain_h / 4.0
		prot_len = max((gene.end - gene.start) // 3, 1)
		span = gx1 - gx0

		def res_to_x(res):
			frac = min(max(res / prot_len, 0.0), 1.0)
			return gx1 - frac * span if minus else gx0 + frac * span

		wedges, labels = [], []
		for d in hits:
			s, e = res_to_x(d.start), res_to_x(d.end)
			color = overlay["color"][d.name]
			pts = [(s, cy + 2 * u), (s, cy + 1 * u), (e, cy - 2 * u),
				   (e, cy + 2 * u), (s, cy + 2 * u)]
			points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
			wedges.append('<polygon points="{}" fill="{}"/>'.format(points, color))
			labels.append(((s + e) / 2, overlay["number"][d.name], color))
		return "".join(wedges), labels

	def _features_on_gene(self, gene, gx0, gx1, cy, drawn_strand, regions):
		if not regions:
			return ""
		minus = (drawn_strand == "-")
		prot_len = max((gene.end - gene.start) // 3, 1)
		span = gx1 - gx0

		def res_to_x(res):
			frac = min(max(res / prot_len, 0.0), 1.0)
			return gx1 - frac * span if minus else gx0 + frac * span

		h = self.gene_h
		out = []
		for kind, start, end in regions:
			x_s, x_e = res_to_x(start), res_to_x(end)
			lo, hi = (x_e, x_s) if x_e < x_s else (x_s, x_e)
			if kind == "tm":
				top, bot = cy - h / 2, cy + h / 2
				n = 4                                   # number of hatch lines
				for k in range(1, n + 1):
					hy = top + (bot - top) * k / (n + 1)
					out.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" '
							   'stroke="#d40000" stroke-width="0.8"/>'.format(lo, hy, hi, hy))
				mid = (lo + hi) / 2
				out.append('<line x1="{:.1f}" y1="{:.1f}" x2="{:.1f}" y2="{:.1f}" '
						   'stroke="#fff" stroke-width="1.2"/>'.format(mid, top, mid, bot))
			elif kind == "signal":
				tip = lo
				w = max(hi - lo, 3) * 0.6
				pts = [(tip, cy - h / 2), (tip - w / 2, cy + h / 2),
					   (tip + w / 2, cy + h / 2)]
				points = " ".join("{:.1f},{:.1f}".format(px, py) for px, py in pts)
				out.append('<polygon points="{}" fill="#000"/>'.format(points))
		return "".join(out)

	def _place_labels(self, row_labels, y):
		step = self.LABEL_STEP
		ly = y - self.gene_h / 2 - 4
		out, placed = [], []
		for lx, num, col in sorted(row_labels, key=lambda t: t[0]):
			if placed and lx - placed[-1] < step:
				lx = placed[-1] + step
			placed.append(lx)
			out.append('<text x="{:.1f}" y="{:.1f}" font-size="{}" fill="{}" '
					   'text-anchor="middle">{}</text>'.format(
						   lx, ly, self.font - 4, col, num))
		return "".join(out)

	def _gene_style(self, gene, gene_color):
		special = self._special_type(gene.accession)
		if special == "pseudo":
			return self.PSEUDO
		if special == "rna":
			return self.RNA
		if special == "other":
			return self.OTHER
		outline = self.RNA[1] if gene.is_rna else "#333"
		if self.mode == "domains":
			return "rgba(0,0,0,0)", outline
		if self.mode == "secretion":
			return "#ffffff", outline
		fill = gene_color.get(gene.accession, self.RNA[0] if gene.is_rna else self.GREY)
		return fill, outline

	def _legend_layout(self, overlay, W):
		items = overlay.get("legend")
		if not items:
			return None, 0, 0, 0
		longest = max(self._text_width("{}. {}".format(num, lbl), self.font)
					  for num, lbl, _ in items)
		col_w = int(longest) + 26
		avail = W - 2 * self.pad
		cols = max(1, min(len(items), avail // col_w))
		legend_rows = -(-len(items) // cols)
		legend_h = (legend_rows + 1) * 16 + 8
		return items, cols, col_w, legend_h

	def _legend(self, items, x, y, cols, col_w, title="Domains"):
		parts = ['<text x="{}" y="{}" font-weight="bold">{}</text>'.format(x, y, title)]
		for idx, (num, label, color) in enumerate(items):
			ex = x + (idx % cols) * col_w
			ey = y + 16 + (idx // cols) * 16
			parts.append('<rect x="{}" y="{}" width="12" height="10" fill="{}"/>'.format(
				ex + 2, ey - 8, color))
			parts.append('<text x="{}" y="{}">{}. {}</text>'.format(
				ex + 20, ey, num, self._escape(label)))
		return "".join(parts)

	@staticmethod
	def _contrast_color(n):
		hue = (n * 0.61803398875) % 1.0
		r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.75)
		return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))

	@staticmethod
	def _row_reversed(genes):
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