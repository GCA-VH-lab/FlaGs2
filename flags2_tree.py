import os
import shutil
import subprocess
import tempfile
from io import StringIO
from typing import Dict, List, Optional, Tuple

from Bio import Phylo, SeqIO

from flags2_view import _FlaGsBase
from FlaGs2 import FlankingGene


class TreeBuilder:

  GAP_CHARS = "-."

  def __init__(self, threads: int = 0, engine: str = "veryfasttree",
      gap_threshold: float = 0.1):
    self.threads = threads
    self.engine = engine
    self.gap_threshold = gap_threshold
    self.alignment: Dict[str, str] = {}

  def build(self, sequences: Dict[str, str]) -> Tuple[str, List[str]]:
    names = list(sequences)
    if len(names) < 3:
      return "", names

    try:
      with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "q.fasta")
        aln = os.path.join(tmp, "q.aln")
        with open(fasta, "w") as out:
          for name, seq in sequences.items():
            out.write(">{}\n{}\n".format(name, seq))

        with open(aln, "w") as out:
          subprocess.run(["mafft", "--auto", "--anysymbol", "--quiet",
          "--thread", str(self.threads), fasta],
          check=True, stdout=out, stderr=subprocess.DEVNULL)

        self.alignment = self._trim(self._read_alignment(aln), self.gap_threshold)
        trimmed = os.path.join(tmp, "q.trimmed.aln")
        with open(trimmed, "w") as out:
          for name, seq in self.alignment.items():
            out.write(">{}\n{}\n".format(name, seq))

        if self.engine == "iqtree":
          newick = self._run_iqtree(trimmed, tmp, len(names))
        else:
          newick = subprocess.run(["VeryFastTree", trimmed],
            check=True, capture_output=True, text=True).stdout.strip()
    except FileNotFoundError as e:
      print("Warning: tree building needs mafft and {} on PATH; skipping the tree "
            "({}).".format("iqtree" if self.engine == "iqtree" else "VeryFastTree", e))
      return "", names
    except subprocess.CalledProcessError as e:
      print("Warning: tree building failed, skipping the tree ({}).".format(e))
      return "", names

    leaf_order = [t.name for t in Phylo.read(StringIO(newick), "newick").get_terminals()]
    return newick, leaf_order

  def _run_iqtree(self, aln: str, tmp: str, n_taxa: int) -> str:
    binary = next((b for b in ("iqtree3", "iqtree2", "iqtree")
                   if shutil.which(b)), None)
    if binary is None:
      raise FileNotFoundError("no iqtree binary found")
    cmd = [binary, "-s", aln, "-m", "MFP", "--prefix", os.path.join(tmp, "iq"),
           "-T", str(self.threads) if self.threads else "AUTO", "--quiet"]
    if n_taxa >= 4:
      cmd += ["-B", "1000"]   # ultrafast bootstrap needs at least 4 taxa
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    with open(os.path.join(tmp, "iq.treefile")) as fh:
      return fh.read().strip()

  @staticmethod
  def _read_alignment(path: str) -> Dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}

  @classmethod
  def _trim(cls, alignment: Dict[str, str], gap_threshold: float) -> Dict[str, str]:
    """Drop columns where fewer than gap_threshold of sequences carry a residue.

    Equivalent to trimal -gt, which is what ete3's trimal01 ran in the old pipeline.
    """
    rows = list(alignment.values())
    if not rows:
      return alignment
    width = len(rows[0])
    need = gap_threshold * len(rows)
    keep = [i for i in range(width)
            if sum(1 for r in rows if r[i] not in cls.GAP_CHARS) >= need]
    if not keep or len(keep) == width:
      return alignment
    return {name: "".join(seq[i] for i in keep) for name, seq in alignment.items()}


class NeighborhoodVisualizer(_FlaGsBase): 

  def __init__(self, gene_h: int = 20, gene_gap: int = 1,
      row_h: int = 24, tree_w: int = 320, pad: int = 16, font: int = 13):
    self.gene_h = gene_h
    self.gene_w = gene_h * 0.95
    self.gene_gap = gene_gap
    self.row_h = row_h; self.tree_w = tree_w; self.pad = pad
    self.font = font

  def render(self, newick: str,
      neighborhoods: List["FlankingGene"],
      families: List[List[str]],
      species: Optional[Dict[str, str]] = None,
      labels: Optional[Dict[str, str]] = None) -> str:
    species = species or {}
    row_labels_map = labels or {}
    by_query: Dict[str, list] = {}
    for g in neighborhoods:
      by_query.setdefault(g.query, []).append(g)
    for q in by_query:
      by_query[q].sort(key=lambda x: x.offset)

    color = self._family_colors(families)
    rna_accessions = {g.accession for g in neighborhoods if g.is_rna}
    number = self._family_numbers(families, rna_accessions)
    tree = Phylo.read(StringIO(newick), "newick")
    try:
      tree.root_at_midpoint()
    except Exception:
      pass    
    tree.ladderize()
    tip_order = [t.name for t in tree.get_terminals()]
    rows = [q for q in tip_order if q in by_query] or list(by_query)
    def _lbl(q):
      if q in row_labels_map:
        return row_labels_map[q]
      return q if not species.get(q) else "{}  {}".format(q, species[q])
    labels = {q: _lbl(q) for q in rows}
    label_w = int(max((self._text_width(l, self.font) for l in labels.values()), default=0)) + 14

    cell = self.gene_w + self.gene_gap
    max_off = max((abs(g.offset) for gs in by_query.values() for g in gs), default=0)
    track_w = (2 * max_off + 1) * cell

    W = self.pad + self.tree_w + label_w + track_w + self.pad
    H = self.pad * 2 + len(rows) * self.row_h + 30

    y_of = {q: self.pad + i * self.row_h + self.row_h / 2 for i, q in enumerate(rows)}

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" '
      'font-family="{}" font-size="{}">'.format(W, H, self.FONT_FAMILY, self.font),
      '<rect width="{}" height="{}" fill="white"/>'.format(W, H)]

    tree_svg, xscale, maxd = self._draw_tree(tree, y_of)
    svg.append(tree_svg)

    track_x0 = self.pad + self.tree_w + label_w
    center = track_x0 + track_w / 2  

    for q in rows:
      y = y_of[q]
      svg.append('<text x="{}" y="{}">{}</text>'.format(
        self.pad + self.tree_w + 6, y + self.font / 3, self._escape(labels[q])))
      for g in by_query[q]:
        cx = center + g.offset * cell
        fill, outline = self._style_for(g.accession, color, g.is_rna)
        stroke = "#000" if g.offset == 0 else outline
        sw = 2 if g.offset == 0 else 1
        svg.append(self._arrow(cx, y, g.strand, fill, stroke, sw,
          number.get(g.accession)))

    base_y = self.pad + len(rows) * self.row_h + 16
    if xscale and maxd:
      svg.append(self._scale_bar(self.pad, base_y, xscale, maxd))
    svg.append('</svg>')
    return "\n".join(svg)

  def _scale_bar(self, x: float, y: float, xscale: float, maxd: float) -> str:
    import math
    target = maxd / 5 if maxd else 0.1
    if target <= 0:
      return ""
    mag = 10 ** math.floor(math.log10(target))
    nice = min((1, 2, 5, 10), key=lambda m: abs(m * mag - target)) * mag
    w = nice * xscale
    return ('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#333" stroke-width="1.5"/>'
        '<line x1="{0}" y1="{3}" x2="{0}" y2="{4}" stroke="#333"/>'
        '<line x1="{2}" y1="{3}" x2="{2}" y2="{4}" stroke="#333"/>'
        '<text x="{5}" y="{6}" text-anchor="middle" font-size="{7}">{8:g}</text>'
        ).format(x, y, x + w, y - 3, y + 3, x + w / 2, y + self.font + 2,
          self.font - 1, nice)

  def _arrow(self, cx, cy, strand, fill, stroke, sw, number=None) -> str:
    hw, hh = self.gene_w / 2, self.gene_h / 2
    if strand == "+":
      pts = [(cx - hw, cy - hh), (cx + hw, cy), (cx - hw, cy + hh)]
    else:
      pts = [(cx + hw, cy - hh), (cx - hw, cy), (cx + hw, cy + hh)]
    points = " ".join("{},{}".format(x, y) for x, y in pts)
    out = ['<polygon points="{}" fill="{}" stroke="{}" stroke-width="{}"/>'.format(
      points, fill, stroke, sw)]
    if number is not None:
      tx = cx - hw * 0.2 if strand == "+" else cx + hw * 0.2
      out.append('<text x="{}" y="{}" font-size="{}" fill="black" '
        'text-anchor="middle" dominant-baseline="central">{}</text>'.format(
          tx, cy, self.font - 5, number))
    return "".join(out)

  def _draw_tree(self, tree, y_of: Dict[str, float]) -> str:
    depths = tree.depths()
    if not any(depths.values()):
      depths = tree.depths(unit_branch_lengths=True)
    maxd = max(depths.values()) or 1
    xscale = (self.tree_w - 10) / maxd
    x0 = self.pad
    yc: Dict = {}

    def assign(clade):
      if clade.is_terminal():
        yc[clade] = y_of.get(clade.name, self.pad)
      else:
        ys = [assign(c) for c in clade.clades]
        yc[clade] = sum(ys) / len(ys)
      return yc[clade]
    assign(tree.root)

    seg = []
    def walk(clade, px):
      x = x0 + depths[clade] * xscale
      y = yc[clade]
      seg.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(px, y, x, y))
      if clade.is_terminal():
        seg.append('<circle cx="{}" cy="{}" r="2.2" fill="#555"/>'.format(x, y))
      else:
        cys = [yc[c] for c in clade.clades]
        seg.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(
          x, min(cys), x, max(cys)))
        if clade.confidence is not None:
          seg.append('<text x="{}" y="{}" font-size="8" fill="#8b0000" '
            'text-anchor="end">{:g}</text>'.format(x - 4, y + 9, clade.confidence))
        for c in clade.clades:
          walk(c, x)
    walk(tree.root, x0)
    return "\n".join(seg), xscale, maxd


def ladderized_leaf_order(newick):
  t = Phylo.read(StringIO(newick), "newick")
  try:
    t.root_at_midpoint()
  except Exception:
    pass
  t.ladderize()
  return [tip.name for tip in t.get_terminals()]