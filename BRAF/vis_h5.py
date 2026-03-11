"""
HDF5 → Interactive HTML Explorer
Generates a self-contained HTML file with collapsible tree + data preview.
Usage: python vis_h5.py [path.h5]
"""

import sys, os, json, base64
import numpy as np
import h5py

FILE = (
    "/Users/kateplas/Library/Mobile Documents/"
    "com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_H.dna+protein.h5"
)

# ── Helpers ───────────────────────────────────────────────────────────────────

def decode_val(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    if isinstance(v, np.generic):
        return v.item()
    return v

def dataset_preview(ds, max_rows=10, max_cols=12):
    """Return a dict describing the dataset and a small preview."""
    shape  = list(ds.shape)
    dtype  = str(ds.dtype)
    nbytes = int(np.prod(shape)) * ds.dtype.itemsize if shape else ds.dtype.itemsize

    if not shape:  # scalar
        raw = ds[()]
        return {"shape": [], "dtype": dtype, "bytes": nbytes,
                "kind": "scalar", "value": str(decode_val(raw))}

    ndim = len(shape)

    # String / bytes arrays
    if ds.dtype.kind in ("S", "O", "U"):
        flat = ds[()].ravel()[:max_rows * max_cols]
        vals = [str(decode_val(x)) for x in flat]
        return {"shape": shape, "dtype": dtype, "bytes": nbytes,
                "kind": "strings", "values": vals[:50]}

    # Numeric
    if ndim == 1:
        n   = min(max_rows * max_cols, shape[0])
        idx = np.sort(np.random.choice(shape[0], n, replace=False)) if shape[0] > n else np.arange(shape[0])
        vals = ds[idx].tolist()
        return {"shape": shape, "dtype": dtype, "bytes": nbytes,
                "kind": "1d", "values": vals}

    if ndim == 2:
        r = min(max_rows,  shape[0])
        c = min(max_cols,  shape[1])
        ridx = np.sort(np.random.choice(shape[0], r, replace=False)) if shape[0] > r else np.arange(shape[0])
        data = ds[ridx, :c].tolist()
        # stats
        sample = ds[ridx, :].astype(float)
        stats  = {
            "min":    float(np.nanmin(sample)),
            "max":    float(np.nanmax(sample)),
            "mean":   float(np.nanmean(sample)),
            "median": float(np.nanmedian(sample)),
        }
        return {"shape": shape, "dtype": dtype, "bytes": nbytes,
                "kind": "2d", "rows_shown": r, "cols_shown": c,
                "data": data, "stats": stats}

    # ≥3-D: show a 2-D slice
    mid = [s // 2 for s in shape[:-2]]
    idx_expr = tuple(mid) + (slice(0, min(max_rows, shape[-2])),
                              slice(0, min(max_cols, shape[-1])))
    sl = ds[idx_expr].tolist()
    return {"shape": shape, "dtype": dtype, "bytes": nbytes,
            "kind": "nd_slice", "slice_axes": mid, "data": sl}


def walk_h5(f):
    """Return nested dict tree."""
    def _node(name, obj):
        short = name.split("/")[-1] if name != "/" else "/"
        if isinstance(obj, h5py.Group):
            children = []
            for k in obj.keys():
                child_path = f"{name}/{k}".lstrip("/")
                children.append(_node(child_path, obj[k]))
            return {"type": "group", "name": short, "path": name,
                    "children": children}
        else:  # Dataset
            preview = dataset_preview(obj)
            return {"type": "dataset", "name": short, "path": name,
                    "preview": preview}

    root = {"type": "group", "name": os.path.basename(f.filename),
            "path": "/", "children": []}
    for k in f.keys():
        root["children"].append(_node(k, f[k]))
    return root

# ── HTML template ─────────────────────────────────────────────────────────────

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HDF5 Explorer</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Syne:wght@400;700;800&display=swap');

  :root {
    --bg:       #080c10;
    --panel:    #0d1117;
    --border:   #1e2633;
    --hover:    #131c26;
    --grp:      #4fa3e8;
    --ds:       #3fb950;
    --dim:      #c9921a;
    --txt:      #cdd9e5;
    --muted:    #545d68;
    --accent:   #58a6ff;
    --badge-bg: #182030;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background: var(--bg);
    color: var(--txt);
    font-family: 'JetBrains Mono', monospace;
    font-size: 13px;
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }

  /* ── Header ── */
  header {
    padding: 14px 24px;
    border-bottom: 1px solid var(--border);
    display: flex;
    align-items: center;
    gap: 16px;
    background: var(--panel);
    flex-shrink: 0;
  }
  header h1 {
    font-family: 'Syne', sans-serif;
    font-weight: 800;
    font-size: 17px;
    color: var(--accent);
    letter-spacing: -0.3px;
  }
  header .filename {
    color: var(--muted);
    font-size: 11px;
  }
  .badge {
    background: var(--badge-bg);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 10px;
    color: var(--muted);
  }
  .search-wrap {
    margin-left: auto;
    position: relative;
  }
  #search {
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--txt);
    font-family: 'JetBrains Mono', monospace;
    font-size: 12px;
    padding: 5px 10px 5px 28px;
    width: 220px;
    outline: none;
    transition: border-color .15s;
  }
  #search:focus { border-color: var(--accent); }
  .search-wrap::before {
    content: "⌕";
    position: absolute;
    left: 8px;
    top: 50%;
    transform: translateY(-50%);
    color: var(--muted);
    font-size: 14px;
    pointer-events: none;
  }

  /* ── Layout ── */
  .layout {
    display: flex;
    flex: 1;
    overflow: hidden;
  }

  /* ── Tree panel ── */
  .tree-panel {
    width: 380px;
    min-width: 240px;
    max-width: 560px;
    resize: horizontal;
    overflow: auto;
    border-right: 1px solid var(--border);
    padding: 10px 0;
    flex-shrink: 0;
  }

  ul { list-style: none; padding: 0; }
  ul ul { padding-left: 18px; border-left: 1px solid #1e2633; margin-left: 10px; }

  .node {
    cursor: pointer;
    padding: 3px 14px 3px 6px;
    border-radius: 4px;
    display: flex;
    align-items: center;
    gap: 6px;
    white-space: nowrap;
    transition: background .1s;
  }
  .node:hover { background: var(--hover); }
  .node.active { background: #132040; }
  .node.hidden { display: none; }

  .toggle {
    width: 14px;
    height: 14px;
    flex-shrink: 0;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 9px;
    color: var(--muted);
    transition: transform .15s;
  }
  .toggle.open { transform: rotate(90deg); }
  .toggle.leaf  { color: transparent; }

  .icon { font-size: 12px; flex-shrink: 0; }
  .icon.grp { color: var(--grp); }
  .icon.ds  { color: var(--ds); }

  .label { color: var(--txt); font-size: 12px; flex: 1; min-width: 0; overflow: hidden; text-overflow: ellipsis; }
  .shape-badge {
    font-size: 9px;
    color: var(--dim);
    background: #1a1200;
    border: 1px solid #3a2800;
    border-radius: 3px;
    padding: 1px 5px;
    flex-shrink: 0;
  }

  .children { overflow: hidden; }
  .children.collapsed { display: none; }

  /* ── Preview panel ── */
  .preview-panel {
    flex: 1;
    overflow: auto;
    padding: 24px 28px;
  }

  .placeholder {
    height: 100%;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    color: var(--muted);
    gap: 10px;
    font-size: 13px;
  }
  .placeholder .big { font-size: 36px; }

  .preview-header {
    margin-bottom: 20px;
    padding-bottom: 14px;
    border-bottom: 1px solid var(--border);
  }
  .preview-header h2 {
    font-family: 'Syne', sans-serif;
    font-weight: 700;
    font-size: 16px;
    color: var(--accent);
    margin-bottom: 4px;
  }
  .preview-header .path {
    color: var(--muted);
    font-size: 11px;
  }

  .meta-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(160px, 1fr));
    gap: 10px;
    margin-bottom: 22px;
  }
  .meta-card {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 10px 14px;
  }
  .meta-card .key { font-size: 10px; color: var(--muted); margin-bottom: 4px; text-transform: uppercase; letter-spacing: .5px; }
  .meta-card .val { font-size: 13px; color: var(--txt); font-weight: 600; }

  .section-title {
    font-family: 'Syne', sans-serif;
    font-size: 11px;
    font-weight: 700;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 10px;
  }

  /* ── Stats row ── */
  .stats-row {
    display: flex;
    gap: 10px;
    margin-bottom: 22px;
    flex-wrap: wrap;
  }
  .stat-pill {
    background: var(--badge-bg);
    border: 1px solid var(--border);
    border-radius: 20px;
    padding: 4px 14px;
    font-size: 11px;
    color: var(--txt);
  }
  .stat-pill span { color: var(--dim); margin-right: 4px; }

  /* ── Table ── */
  .table-wrap {
    overflow: auto;
    border: 1px solid var(--border);
    border-radius: 8px;
    margin-bottom: 10px;
  }
  table {
    border-collapse: collapse;
    width: max-content;
    min-width: 100%;
  }
  th {
    background: #0d1520;
    color: var(--dim);
    font-size: 10px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: .4px;
    padding: 7px 12px;
    border-bottom: 1px solid var(--border);
    text-align: right;
    white-space: nowrap;
    position: sticky;
    top: 0;
  }
  th:first-child { text-align: left; color: var(--muted); }
  td {
    padding: 5px 12px;
    border-bottom: 1px solid #111820;
    font-size: 11.5px;
    text-align: right;
    white-space: nowrap;
    color: var(--txt);
  }
  td:first-child { color: var(--muted); text-align: left; }
  tr:hover td { background: var(--hover); }
  .truncate-note { font-size: 10px; color: var(--muted); margin-top: 6px; }

  /* ── String list ── */
  .string-grid {
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    margin-bottom: 22px;
  }
  .str-chip {
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 4px;
    padding: 3px 10px;
    font-size: 11px;
    color: var(--ds);
  }

  /* ── 1-D sparkline bar ── */
  .bar-wrap {
    display: flex;
    flex-wrap: wrap;
    gap: 2px;
    margin-bottom: 22px;
    align-items: flex-end;
    height: 60px;
  }
  .bar-item {
    flex: 1;
    min-width: 3px;
    max-width: 18px;
    background: var(--accent);
    opacity: 0.7;
    border-radius: 2px 2px 0 0;
    transition: opacity .1s;
  }
  .bar-item:hover { opacity: 1; }

  /* scrollbar */
  ::-webkit-scrollbar { width: 6px; height: 6px; }
  ::-webkit-scrollbar-track { background: transparent; }
  ::-webkit-scrollbar-thumb { background: #2a3240; border-radius: 3px; }
</style>
</head>
<body>

<header>
  <div>
    <h1>HDF5 Explorer</h1>
    <div class="filename" id="fname"></div>
  </div>
  <div class="badge" id="node-count"></div>
  <div class="search-wrap">
    <input id="search" type="text" placeholder="Search nodes…" oninput="filterTree(this.value)">
  </div>
</header>

<div class="layout">
  <div class="tree-panel" id="tree"></div>
  <div class="preview-panel" id="preview">
    <div class="placeholder">
      <div class="big">⬡</div>
      <div>Click any node to preview its data</div>
    </div>
  </div>
</div>

<script>
const DATA = __DATA__;

// ── Count nodes ───────────────────────────────────────────────────────────────
let totalNodes = 0;
function countNodes(n) { totalNodes++; if (n.children) n.children.forEach(countNodes); }
countNodes(DATA);
document.getElementById('fname').textContent = DATA.name;
document.getElementById('node-count').textContent = totalNodes + ' nodes';

// ── Build tree ────────────────────────────────────────────────────────────────
function buildTree(node, container) {
  const li = document.createElement('li');

  const row = document.createElement('div');
  row.className = 'node';
  row.dataset.path = node.path;
  row.dataset.name = node.name.toLowerCase();

  const toggle = document.createElement('span');
  toggle.className = 'toggle';

  const icon = document.createElement('span');
  icon.className = 'icon';

  const label = document.createElement('span');
  label.className = 'label';
  label.textContent = node.name;

  if (node.type === 'group') {
    toggle.textContent = '▶';
    icon.textContent = '⬡';
    icon.classList.add('grp');

    const childUl = document.createElement('ul');
    childUl.className = 'children collapsed';

    row.addEventListener('click', e => {
      e.stopPropagation();
      const open = !childUl.classList.contains('collapsed');
      childUl.classList.toggle('collapsed', open);
      toggle.classList.toggle('open', !open);
      document.querySelectorAll('.node.active').forEach(n => n.classList.remove('active'));
      row.classList.add('active');
      showGroupPreview(node);
    });

    node.children.forEach(child => buildTree(child, childUl));
    li.appendChild(row);
    li.appendChild(childUl);

  } else {
    toggle.className = 'toggle leaf';
    icon.textContent = '▪';
    icon.classList.add('ds');

    const p = node.preview;
    if (p && p.shape && p.shape.length > 0) {
      const sb = document.createElement('span');
      sb.className = 'shape-badge';
      sb.textContent = p.shape.join('×');
      row.appendChild(toggle); row.appendChild(icon); row.appendChild(label); row.appendChild(sb);
    }

    row.addEventListener('click', e => {
      e.stopPropagation();
      document.querySelectorAll('.node.active').forEach(n => n.classList.remove('active'));
      row.classList.add('active');
      showDatasetPreview(node);
    });

    if (!row.children.length) {
      row.appendChild(toggle); row.appendChild(icon); row.appendChild(label);
    }
    li.appendChild(row);
    return container.appendChild(li);
  }

  row.insertBefore(toggle, row.firstChild);
  row.insertBefore(icon,   toggle.nextSibling);
  row.insertBefore(label,  icon.nextSibling);
  container.appendChild(li);
}

const rootUl = document.createElement('ul');
buildTree(DATA, rootUl);
document.getElementById('tree').appendChild(rootUl);

// auto-expand root
const firstToggle = rootUl.querySelector('.toggle');
const firstChildren = rootUl.querySelector('.children');
if (firstToggle && firstChildren) {
  firstChildren.classList.remove('collapsed');
  firstToggle.classList.add('open');
}

// ── Search ────────────────────────────────────────────────────────────────────
function filterTree(q) {
  q = q.toLowerCase().trim();
  document.querySelectorAll('.node').forEach(row => {
    if (!q) { row.classList.remove('hidden'); return; }
    const match = row.dataset.name && row.dataset.name.includes(q);
    row.classList.toggle('hidden', !match);
    if (match) {
      // reveal parents
      let p = row.parentElement;
      while (p) {
        if (p.classList && p.classList.contains('children')) p.classList.remove('collapsed');
        if (p.classList && p.classList.contains('node')) p.classList.remove('hidden');
        p = p.parentElement;
      }
    }
  });
}

// ── Format helpers ────────────────────────────────────────────────────────────
function fmtBytes(b) {
  if (b < 1024) return b + ' B';
  if (b < 1048576) return (b/1024).toFixed(1) + ' KB';
  if (b < 1073741824) return (b/1048576).toFixed(1) + ' MB';
  return (b/1073741824).toFixed(2) + ' GB';
}
function fmtNum(v) {
  if (v === null || v === undefined) return '—';
  if (typeof v === 'number') {
    if (!isFinite(v)) return String(v);
    return Math.abs(v) > 1e4 || (Math.abs(v) < 1e-3 && v !== 0)
      ? v.toExponential(3) : v.toPrecision(5).replace(/\.?0+$/, '');
  }
  return String(v);
}

// ── Group preview ─────────────────────────────────────────────────────────────
function showGroupPreview(node) {
  const groups = node.children.filter(c => c.type === 'group').length;
  const ds     = node.children.filter(c => c.type === 'dataset').length;
  document.getElementById('preview').innerHTML = `
    <div class="preview-header">
      <h2>${node.name}</h2>
      <div class="path">${node.path}</div>
    </div>
    <div class="meta-grid">
      <div class="meta-card"><div class="key">Type</div><div class="val">Group</div></div>
      <div class="meta-card"><div class="key">Subgroups</div><div class="val">${groups}</div></div>
      <div class="meta-card"><div class="key">Datasets</div><div class="val">${ds}</div></div>
      <div class="meta-card"><div class="key">Total children</div><div class="val">${node.children.length}</div></div>
    </div>
    <div class="section-title">Children</div>
    <div class="string-grid">${node.children.map(c =>
      `<span class="str-chip" style="color:${c.type==='group'?'var(--grp)':'var(--ds)'}">
        ${c.type==='group'?'⬡':'▪'} ${c.name}
       </span>`).join('')}
    </div>`;
}

// ── Dataset preview ───────────────────────────────────────────────────────────
function showDatasetPreview(node) {
  const p = node.preview;
  let html = `
    <div class="preview-header">
      <h2>${node.name}</h2>
      <div class="path">${node.path}</div>
    </div>
    <div class="meta-grid">
      <div class="meta-card"><div class="key">Shape</div><div class="val">${p.shape && p.shape.length ? p.shape.join(' × ') : 'scalar'}</div></div>
      <div class="meta-card"><div class="key">dtype</div><div class="val">${p.dtype}</div></div>
      <div class="meta-card"><div class="key">Size</div><div class="val">${fmtBytes(p.bytes)}</div></div>
      <div class="meta-card"><div class="key">Kind</div><div class="val">${p.kind}</div></div>
    </div>`;

  if (p.kind === 'scalar') {
    html += `<div class="section-title">Value</div>
      <div class="str-chip" style="font-size:14px;color:var(--txt)">${p.value}</div>`;

  } else if (p.kind === 'strings') {
    html += `<div class="section-title">Values (first ${p.values.length})</div>
      <div class="string-grid">${p.values.map(v =>
        `<span class="str-chip">${v}</span>`).join('')}
      </div>`;

  } else if (p.kind === '1d') {
    const nums = p.values.filter(v => typeof v === 'number' && isFinite(v));
    const mn = nums.length ? Math.min(...nums) : 0;
    const mx = nums.length ? Math.max(...nums) : 1;
    const range = mx - mn || 1;
    html += `<div class="section-title">Distribution preview (${p.values.length} samples)</div>
      <div class="bar-wrap">${p.values.map(v => {
        const h = Math.max(4, Math.round(((v - mn) / range) * 55));
        return `<div class="bar-item" style="height:${h}px" title="${fmtNum(v)}"></div>`;
      }).join('')}</div>`;
    if (nums.length) {
      const mean = nums.reduce((a,b)=>a+b,0)/nums.length;
      html += `<div class="stats-row">
        <div class="stat-pill"><span>min</span>${fmtNum(mn)}</div>
        <div class="stat-pill"><span>max</span>${fmtNum(mx)}</div>
        <div class="stat-pill"><span>mean</span>${fmtNum(mean)}</div>
      </div>`;
    }
    html += `<div class="section-title">Raw values</div>
      <div class="string-grid">${p.values.map(v =>
        `<span class="str-chip" style="color:var(--txt)">${fmtNum(v)}</span>`).join('')}
      </div>`;

  } else if (p.kind === '2d' || p.kind === 'nd_slice') {
    if (p.stats) {
      html += `<div class="stats-row">
        <div class="stat-pill"><span>min</span>${fmtNum(p.stats.min)}</div>
        <div class="stat-pill"><span>max</span>${fmtNum(p.stats.max)}</div>
        <div class="stat-pill"><span>mean</span>${fmtNum(p.stats.mean)}</div>
        <div class="stat-pill"><span>median</span>${fmtNum(p.stats.median)}</div>
      </div>`;
    }
    if (p.kind === 'nd_slice') {
      html += `<div class="section-title">Slice at axes ${JSON.stringify(p.slice_axes)}</div>`;
    } else {
      html += `<div class="section-title">Preview — ${p.rows_shown} rows × ${p.cols_shown} cols shown</div>`;
    }
    const nCols = p.data[0].length;
    html += `<div class="table-wrap"><table>
      <thead><tr><th>row</th>${Array.from({length:nCols},(_,i)=>`<th>col ${i}</th>`).join('')}</tr></thead>
      <tbody>${p.data.map((row,ri)=>
        `<tr><td>${ri}</td>${row.map(v=>`<td>${fmtNum(v)}</td>`).join('')}</tr>`
      ).join('')}</tbody>
    </table></div>`;
    if (p.rows_shown) {
      html += `<div class="truncate-note">Showing ${p.rows_shown} of ${p.shape[0]} rows, ${p.cols_shown} of ${p.shape[1]} cols — random sample</div>`;
    }
  }

  document.getElementById('preview').innerHTML = html;
}
</script>
</body>
</html>"""

# ── Generate ───────────────────────────────────────────────────────────────────

fp = sys.argv[1] if len(sys.argv) > 1 else FILE
if not os.path.exists(fp):
    sys.exit(f"File not found: {fp}")

print(f"Reading {fp} …")
np.random.seed(42)

with h5py.File(fp, "r") as f:
    tree = walk_h5(f)

data_json = json.dumps(tree, allow_nan=False)
out_html  = HTML.replace("__DATA__", data_json)

out_path = os.path.join(os.path.dirname(fp), "h5_explorer.html")
with open(out_path, "w", encoding="utf-8") as fh:
    fh.write(out_html)

size_kb = os.path.getsize(out_path) / 1024
print(f"✅  Saved: {out_path}  ({size_kb:.0f} KB)")
print("    Open in any browser — no server needed.")