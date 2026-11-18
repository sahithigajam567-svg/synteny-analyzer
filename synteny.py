"""Core logic for the Comparative Genomics Synteny Analyzer (Biopython)."""
import random
from collections import Counter, defaultdict

from Bio import SeqIO


# ---------- Parsing ----------
def parse_genbank(handle, fmt="genbank"):
    """Read CDS features from a GenBank file handle -> list of gene dicts."""
    genes, offset = [], 0
    for rec in SeqIO.parse(handle, fmt):
        for f in rec.features:
            if f.type != "CDS":
                continue
            q = f.qualifiers
            protein = q.get("translation", [""])[0]
            if not protein:
                continue
            name = q.get("gene", q.get("locus_tag", [""]))[0]
            genes.append({
                "name": name,
                "start": int(f.location.start) + offset,
                "end": int(f.location.end) + offset,
                "strand": f.location.strand or 1,
                "product": q.get("product", [""])[0],
                "protein": protein,
            })
        offset += len(rec.seq)
    genes.sort(key=lambda g: g["start"])
    for i, g in enumerate(genes):
        g["idx"] = i
    return genes


# ---------- Ortholog detection ----------
def orthologs_by_name(a, b):
    """Pair genes that share the same gene name."""
    index = {}
    for g in b:
        if g["name"] and g["name"] not in index:
            index[g["name"]] = g["idx"]
    return [(g["idx"], index[g["name"]], 1.0)
            for g in a if g["name"] in index]


def _kmers(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def _best_hits(query, target, k):
    """Best k-mer Jaccard hit in `target` for each gene in `query`."""
    t_sets = [_kmers(g["protein"], k) for g in target]
    inv = defaultdict(list)
    for j, s in enumerate(t_sets):
        for km in s:
            inv[km].append(j)
    best = {}
    for g in query:
        qs = _kmers(g["protein"], k)
        counts = Counter(j for km in qs for j in inv.get(km, ()))
        if not counts:
            continue
        j, shared = max(counts.items(), key=lambda x: x[1])
        best[g["idx"]] = (j, shared / (len(qs) + len(t_sets[j]) - shared))
    return best


def orthologs_by_similarity(a, b, k=4, min_score=0.05):
    """Reciprocal best hits using protein k-mer similarity (pure Python)."""
    ab, ba = _best_hits(a, b, k), _best_hits(b, a, k)
    return [(i, j, s) for i, (j, s) in ab.items()
            if s >= min_score and ba.get(j, (None,))[0] == i]


# ---------- Synteny blocks ----------
def find_blocks(pairs, min_genes=3, max_gap=5):
    """Chain collinear ortholog pairs into synteny blocks.

    A block is a run of pairs (sorted by genome A order) whose positions in B
    move consistently forward (+) or backward (-, an inversion), with at most
    `max_gap` genes skipped between neighbours.
    """
    pairs = sorted(pairs)
    blocks, cur, direction = [], [], 0

    def close():
        if len(cur) >= min_genes:
            blocks.append({"pairs": list(cur), "direction": direction or 1})

    for p in pairs:
        if not cur:
            cur, direction = [p], 0
            continue
        da, db = p[0] - cur[-1][0], p[1] - cur[-1][1]
        d = 1 if db > 0 else -1
        if db != 0 and da <= max_gap and abs(db) <= max_gap and direction in (0, d):
            cur.append(p)
            direction = d
        else:
            close()
            cur, direction = [p], 0
    close()
    return blocks


# ---------- Demo data (no downloads needed) ----------
def make_demo(n=300, seed=7):
    """Two synthetic genomes with an inversion, a translocation, a deletion and insertions."""
    rnd = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    prot = lambda: "M" + "".join(rnd.choice(aa) for _ in range(rnd.randint(150, 400)))
    mutate = lambda p: "".join(rnd.choice(aa) if rnd.random() < 0.12 else c for c in p)

    base = [(f"gene{i:04d}", prot()) for i in range(n)]
    order = list(range(n))
    order[80:140] = order[80:140][::-1]              # inversion
    seg = order[200:230]; del order[200:230]; order[10:10] = seg  # translocation
    del order[250:270]                               # deletion
    genes_b = [(base[i][0], mutate(base[i][1])) for i in order]
    for pos in (40, 120, 200):                       # insertions of new genes
        genes_b[pos:pos] = [(f"new{pos}_{k}", prot()) for k in range(5)]

    def build(lst):
        out, pos = [], 0
        for i, (name, p) in enumerate(lst):
            L = len(p) * 3
            out.append({"idx": i, "name": name, "start": pos, "end": pos + L,
                        "strand": 1, "product": "hypothetical protein", "protein": p})
            pos += L + 100
        return out
    return build(base), build(genes_b)
