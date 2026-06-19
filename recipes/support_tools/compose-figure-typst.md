---
name: compose-figure-typst
description: Compose publication-ready figures and short documents (figure+caption, multi-panel with letters, manuscript-style pages) by wrapping bare data plots in a Typst layout
when_to_use: User asks for a "Nature-style figure with caption / legend", "multi-panel figure", "panels A/B/C", "Figure 1 / Figure 2 with bold label", "manuscript figure", "compose A and B side by side", or any request that mixes a data plot with typeset prose/labels. Reach for this BEFORE attempting grid/strwrap/gridspec text overlays — those are fragile.
requires_tools: [make_revision, view_artifact]
capabilities_needed: [typst]
keywords: [Nature-style, figure caption, figure legend, multi-panel, panel letters, Figure 1, manuscript, composed PDF, typst, publication-ready, subcaption, bold label]
produces: [composed PDF figure (revision of the current figure in the chain)]
domain: support_tools
---

# Compose figures (and short documents) with Typst

The job here is **typesetting**, not graphics. A bare data plot (PNG / PDF) is
the figure's *content*; arranging it with a caption, panel letters, or
sibling panels is *layout*. Doing both in R/Python's plotting stack means
fighting `grid`/`gridspec` for things they were never designed to do
(wrap-and-measure, full-width text with bold runs, multi-paragraph captions).
**Typst handles all of that in one line.**

Use this recipe whenever the user wants a figure that combines a plot with
typeset prose or with sibling panels. Keep the data-plotting code stable;
the layout step is a separate, pure-text composition.

## Why Typst (and not LaTeX or matplotlib subplots)

- Single Rust binary, shipped inside the `typst` PyPI package — no
  `tlmgr` swamp, no system packages.
- Fast compile (sub-second on small figures).
- Readable syntax + readable error messages.
- Vector output (`.pdf`); no rasterization of text.

## Preconditions

```python
# Once per environment — the propose call is auto-approved (PyPI default).
propose_capability(name="typst", archetype="library",
                   summary="Typst (Rust binary + Python wrapper) for typeset figure composition.",
                   tags=["typesetting", "pdf", "figure-composition"])
ensure_capability(name="typst")
```

After `ensure_capability`, the **Python API** (`import typst; typst.compile(...)`)
is available. The PyPI package does NOT install a standalone `typst` CLI on
PATH — every invocation goes through the Python module. From R, you bridge
via `system2()` (see the R-track template below).

You also need a **bare data figure already on disk** — PDF preferred (vector;
scales without resampling), PNG acceptable for raster plots.

## How this plays with the revision chain — IMPORTANT

The composed PDF is a **revision of the entity the user is currently
viewing** (the chain head). NOT of some "bare-figure anchor" upstream.

```
v1 (bare)  ←  v2 (with grid, R)  ←  v3 (composed PDF with caption, this recipe)
                                     wasRevisionOf points HERE  ─ at v2.
```

Pass `entity_id = <id-of-currently-displayed-figure>` to `make_revision`.
This is usually what `read_entity(<result_id>)` reports as the chain head
for the figure member, or what the user is looking at in the panel.

Why this rule: `make_revision` refuses to branch off a non-latest revision
(it would create a fork with `supersede_newer` required). Pinning to the
anchor when the chain has newer entries always fails. The "wasRevisionOf
points at the bare figure" rule from earlier drafts of this recipe is
**wrong in practice**.

If the user wants to iterate the caption later, each new composition
revises the **previous composition** — the chain grows linearly. If they
want to roll a version back, the per-panel `⋯ → Remove this version` menu
hard-deletes a specific entry and re-links the chain (the entity-mgmt
refactor 2026-06-08).

## Language plumbing — Python-chain vs R-chain

The composed PDF must be produced in the language of the **parent
entity's exec record** (i.e. the currently-displayed figure). `make_revision`
runs `modified_code` in that language's kernel. Two tracks:

### Python-track (parent figure was produced by run_python)

Direct `import typst; typst.compile(...)`. The simple case.

### R-track (parent figure was produced by run_r)

R must shell out to Python because typst is a Python module. From R:

```R
# Find the Python interpreter that has typst installed. ABA's R kernel
# is launched on a separate conda tools env whose `python3` does NOT see
# the typst module — `Sys.which("python3")` is the wrong target. The
# kernel boot sets `ABA_PYTHON` to the backend's own venv python (the
# one with all ABA's deps); read THAT.
.py_for_typst <- function() {
  p <- Sys.getenv("ABA_PYTHON")
  if (nzchar(p) && file.exists(p)) return(p)
  stop("ABA_PYTHON env var not set; this should be set by the kernel boot. ",
       "If you see this in a stale session, restart the kernel.")
}

.compile_typst <- function(typ_path, out_pdf, root = NULL) {
  py <- .py_for_typst()
  if (!nzchar(py)) stop("no python3 binary found; ensure_capability(typst) first")
  if (is.null(root)) {
    code <- sprintf("import typst; typst.compile(%s, %s)",
                    shQuote(typ_path), shQuote(out_pdf))
  } else {
    code <- sprintf("import typst; typst.compile(%s, %s, root=%s)",
                    shQuote(typ_path), shQuote(out_pdf), shQuote(root))
  }
  # shQuote(code) is REQUIRED — system2 won't auto-quote args containing
  # parentheses/semicolons (sh would otherwise split on the ";" and parse
  # the "(" as a subshell-open, yielding a confusing "Syntax error: word
  # unexpected" failure).
  out <- system2(py, c("-c", shQuote(code)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0)
    stop("typst compile failed: ", paste(out, collapse = "\n"))
  invisible(out_pdf)
}
```

Keep these two helpers at the top of every R-track composition.

## Path resolution — read this before writing `image("...")`

Typst sandboxes file reads to a **root directory** (default = the cwd
when compile was invoked). Inside the `.typ` source, an `image("...")`
path is resolved relative to that root. **Absolute filesystem paths
will fail** with `file not found (searched at /tmp<your-path>)` unless
you pass `root=` explicitly.

Two working patterns:

```python
# Pattern A (recommended): write the .typ next to the figure, use the
# bare filename. Simplest — works for single-panel L1.
fig_dir = os.path.dirname(BARE_FIG)
typ_path = os.path.join(fig_dir, "compose.typ")
open(typ_path, "w").write(f'#image("{os.path.basename(BARE_FIG)}", ...)')
typst.compile(typ_path, "/tmp/out.pdf", root=fig_dir)

# Pattern B: set root to a directory ABOVE both the .typ and any images
# (use this when panels live in different subdirectories). Reference
# images with paths relative to root, starting with "/".
typst.compile("/tmp/compose.typ", "/tmp/out.pdf", root="/tmp")
# inside compose.typ:  #image("/figures/panel_a.pdf")
#                      #image("/figures/panel_b.pdf")
```

Pick Pattern A unless you have panels in multiple directories.

## Verify the output — IMPORTANT

After `make_revision` succeeds, **call `view_artifact(entity_id=<new
revision id>)` to LOOK at what was actually produced** before reporting
the change to the user. The tool returns the rendered image as a vision
content block, so you SEE the composed PDF (page 1 by default; pass
`page=N` for multi-page manuscripts).

This catches:
- Typst rendering that doesn't match your intent (e.g. ggplot themes
  silently overridden by Seurat's patchwork wrapper, font fallback to
  the wrong family, image clipped because page size doesn't match the
  bare-figure aspect).
- Caption-vs-image misalignment.
- Layout drift when you re-rendered the bare panel and a setting was
  lost.

DO NOT skip this. The user has spent revision cycles asking for the
same fix because the agent reasoned from its own source code instead of
verifying. If `view_artifact` shows the change didn't land, fix the
code and re-revise — don't claim "done" until what you see matches the
user's request.

## Intermediate-artifact handling — IMPORTANT

`make_revision` runs your code and harvests every figure/table artifact it
produces. If you render a bare panel inside the same call (e.g. R-track
re-renders the plot, then composes), **the bare panel will be auto-harvested
as a sibling revision unless you hide it**. To avoid producing two
revisions per call:

- **Render the bare panel to a temp path** outside the working directory:
  Python `tempfile.NamedTemporaryFile(suffix=".pdf").name`; R `tempfile(fileext=".pdf")`.
- **Copy** it next to `compose.typ` only for the Typst step (Typst needs
  the file under its `root`). Use a stable filename so the .typ references
  resolve.
- **`unlink()` / `os.remove()` the staged copy** after compile. The harvester
  scans the working dir at end-of-call; if the bare copy is gone, only the
  composed PDF (which YOU wrote to the working dir) is picked up.

## Layout templates

Four graded levels. Pick the lowest level that fits the user's request;
don't reach for L3 when L1 will do. Each template is a complete,
working Typst document — substitute the placeholder strings, drop into a
`compose.typ`, compile.

### L1 — figure + caption ("Nature-style with caption")

Single panel, multi-line caption with bold `Figure N.` label, full-width
wrap. This replaces the entire grid/strwrap codepath that fails on
caption-bold-overlay edge cases.

Substitute `FIGURE_FILENAME` with the bare filename of the figure (e.g.
`marker_dotplot.pdf`). Place `compose.typ` in the same directory as the
figure, and pass that directory as `root=` to `typst.compile`.

```typst
#set page(width: 7in, height: 5.5in, margin: (x: 0.3in, y: 0.3in))
#set text(font: "Liberation Sans", size: 10pt)
#set par(justify: false, leading: 0.7em)

#figure(
  image("FIGURE_FILENAME", width: 100%),
  caption: [
    *Figure 1.* CAPTION_BODY_HERE — Typst wraps to the full page
    width automatically; bold/italic/`code` work inline; multi-paragraph
    captions are supported (start a new paragraph with a blank line).
  ],
  supplement: none,
  numbering: none,
)
```

Walkthrough:
- `#set page(...)` — total composed-PDF size. Match `width` to the bare
  figure's aspect (e.g. a 7×4 plot → 7in × ~5in to leave room for caption).
- `image(..., width: 100%)` — scales the bare figure to the content area;
  Typst handles both PDF and PNG embeds, vector when possible.
- The caption is a Typst content block (`[...]`); `*bold*`, `_italic_`,
  `` `code` `` work inline.
- `supplement: none, numbering: none` — keeps "Figure 1." as caption text
  rather than auto-numbering it as part of a larger document.

### L2 — multi-panel with letters ("compose A and B side by side")

Two or more panels arranged in a grid, with bold panel letters in the
top-left of each. Each panel can have its own subcaption.

```typst
#set page(width: 7.5in, height: 4.5in, margin: 0.3in)
#set text(font: "Liberation Sans", size: 9pt)

#grid(
  columns: (1fr, 1fr),
  gutter: 12pt,
  // ── A ──
  block[
    #place(top + left, dx: -4pt, dy: -4pt, [*A*])
    #image("PANEL_A_FILENAME", width: 100%)
  ],
  // ── B ──
  block[
    #place(top + left, dx: -4pt, dy: -4pt, [*B*])
    #image("PANEL_B_FILENAME", width: 100%)
  ],
)

#v(8pt)
#par(leading: 0.7em)[
  *Figure 1.* (*A*) DESCRIPTION_OF_PANEL_A. (*B*) DESCRIPTION_OF_PANEL_B.
]
```

If panels live in different directories, switch to Pattern B in the
"Path resolution" section: set `root=` to a common parent and reference
panels with absolute-style paths starting with `/`.

Walkthrough:
- `grid(columns: (1fr, 1fr), ...)` — two equal-width columns. For 2×2,
  use `columns: (1fr, 1fr)` and pass four panels (rows auto-flow).
- `place(top + left, dx, dy, [*A*])` — absolute-position the panel letter
  in the corner of its containing block. Negative `dx`/`dy` nudges it
  outside the image, hovering over the corner.
- Single caption below the grid carrying (*A*)/(*B*) callouts — this is
  the standard scientific convention.

For 3 panels (A | B / C spanning bottom), nest grids:

```typst
#grid(
  rows: (auto, auto),
  row-gutter: 10pt,
  grid(columns: (1fr, 1fr), gutter: 12pt,
       block[...A...], block[...B...]),
  block[...C spans here...],
)
```

### L3 — multi-section / subcaptions

When each panel deserves its own short caption beneath itself, not pooled
into a single bottom caption.

```typst
#set page(width: 7.5in, height: 7in, margin: 0.3in)
#set text(font: "Liberation Sans", size: 9pt)

#grid(
  columns: (1fr, 1fr),
  gutter: 14pt,
  // Each cell is a vertical stack: panel letter, image, subcaption.
  stack(dir: ttb, spacing: 6pt,
    [*A* — Pre-integration UMAP],
    image("PANEL_A_FILENAME", width: 100%),
    text(size: 8pt)[Day-0 (blue) and day-7 (red) form offset islands per
      cell type — a clear batch effect.],
  ),
  stack(dir: ttb, spacing: 6pt,
    [*B* — Post-Harmony UMAP],
    image("PANEL_B_FILENAME", width: 100%),
    text(size: 8pt)[Cells now interleave within shared clusters; the
      sample-skewed islands that remain (clusters 3, 8) are real.],
  ),
)

#v(10pt)
*Figure 1.* OVERALL_CAPTION_TYING_PANELS_TOGETHER.
```

### L4 — manuscript-style multi-figure document

The user asks for "the full story", "manuscript figure" with abstract
+ figures + methods, or wants several composed figures stitched into one
multi-page PDF. Use Typst's document structure.

```typst
#set page(paper: "us-letter", margin: 1in)
#set text(font: "Liberation Serif", size: 11pt)
#set par(justify: true, leading: 0.65em)

#align(center)[
  #text(size: 16pt, weight: "bold")[Title of the analysis]
  #v(4pt)
  Author Names · Date
]

#v(20pt)

= Abstract
ONE_PARAGRAPH_ABSTRACT.

= Results

== Subsection heading
Body text introducing Figure 1.

#figure(
  image("FIG1_FILENAME", width: 90%),
  caption: [*Figure 1.* SHORT_CAPTION.],
  numbering: none,
)

Body text introducing Figure 2.

#figure(
  image("FIG2_FILENAME", width: 90%),
  caption: [*Figure 2.* SHORT_CAPTION.],
  numbering: none,
)

= Methods
METHODS_PROSE.
```

For L4, the composed PDF is conceptually a **new artifact** (a document,
not a figure-with-caption), and you should pin it as its own Result rather
than as a revision of any single figure. The Result's members can include
the constituent bare figures as siblings if useful for navigation.

## End-to-end pattern — Python track (parent figure produced by run_python)

```python
import os, tempfile, typst
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Step 1: render the bare panel to a TEMP path so make_revision's
# auto-harvest doesn't pick it up as a sibling figure. (The harvester
# scans the working dir at end-of-call; tempfile paths live outside it.)
bare_pdf_tmp = tempfile.NamedTemporaryFile(suffix=".pdf", delete=False).name
plt.figure(figsize=(7, 4.5))
# ... your plotting calls ...
plt.savefig(bare_pdf_tmp, bbox_inches="tight")
plt.close("all")

# Step 2: copy the bare PDF into the working dir under a stable filename
# (Typst needs it under root). Stage compose.typ next to it.
bare_local = "bare_panel.pdf"
import shutil; shutil.copyfile(bare_pdf_tmp, bare_local)

caption_body = "YOUR CAPTION TEXT — Typst wraps it to full width."
typst_src = f'''
#set page(width: 7in, height: 5.5in, margin: 0.3in)
#set text(font: "Liberation Sans", size: 10pt)
#set par(leading: 0.7em)

#figure(
  image("{bare_local}", width: 100%),
  caption: [*Figure 1.* {caption_body}],
  supplement: none, numbering: none,
)
'''
with open("compose.typ", "w") as f:
    f.write(typst_src)

# Step 3: compile (root=cwd because compose.typ + bare_local are both there).
out_pdf = "figure_composed.pdf"
typst.compile("compose.typ", out_pdf, root=os.getcwd())

# Step 4: REMOVE the staged bare copy. The harvester will pick up only
# out_pdf as the artifact for this make_revision call.
os.remove(bare_local)
os.remove("compose.typ")            # optional — keeps the wd clean
os.remove(bare_pdf_tmp)             # cleanup tempfile

print("composed:", os.path.getsize(out_pdf), "bytes")
```

Pin via `make_revision(entity_id=<chain-head-figure-id>, modified_code=<the
code above>)`. The chain head is normally the entity the user is currently
viewing (the displayed revision of the figure member in the focused Result).

## End-to-end pattern — R track (parent figure produced by run_r)

```R
suppressPackageStartupMessages({
  library(Seurat)   # or ggplot2 / pagoda2 / ... — whatever your bare plot needs
})

# Helpers — keep these at the top of every R-track composition. See the
# "Language plumbing" section above for why R must shell out to Python.
.py_for_typst <- function() {
  p <- Sys.getenv("ABA_PYTHON")
  if (nzchar(p) && file.exists(p)) return(p)
  stop("ABA_PYTHON env var not set; restart the kernel.")
}
.compile_typst <- function(typ_path, out_pdf, root = NULL) {
  py <- .py_for_typst()
  if (!nzchar(py)) stop("no python3 binary found; ensure_capability(typst) first")
  if (is.null(root)) {
    code <- sprintf("import typst; typst.compile(%s, %s)",
                    shQuote(typ_path), shQuote(out_pdf))
  } else {
    code <- sprintf("import typst; typst.compile(%s, %s, root=%s)",
                    shQuote(typ_path), shQuote(out_pdf), shQuote(root))
  }
  # shQuote(code) is REQUIRED — system2 won't auto-quote args containing
  # parentheses/semicolons (sh would otherwise split on the ";" and parse
  # the "(" as a subshell-open, yielding a confusing "Syntax error: word
  # unexpected" failure).
  out <- system2(py, c("-c", shQuote(code)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0)
    stop("typst compile failed: ", paste(out, collapse = "\n"))
  invisible(out_pdf)
}

# Step 1: render the bare panel to a TEMP path so make_revision's
# auto-harvest doesn't pick it up as a sibling figure.
bare_tmp <- tempfile(fileext = ".pdf")
ggsave(bare_tmp, plot = p, width = 7, height = 4.5, device = cairo_pdf)

# Step 2: copy bare into wd under a stable filename + write compose.typ.
bare_local <- "bare_panel.pdf"
file.copy(bare_tmp, bare_local, overwrite = TRUE)

caption_body <- "YOUR CAPTION TEXT — Typst wraps it to full width."
typst_src <- sprintf('
#set page(width: 7in, height: 5.5in, margin: 0.3in)
#set text(font: "Liberation Sans", size: 10pt)
#set par(leading: 0.7em)

#figure(
  image("%s", width: 100%%),
  caption: [*Figure 1.* %s],
  supplement: none, numbering: none,
)
', bare_local, caption_body)
writeLines(typst_src, "compose.typ")

# Step 3: compile via the helper (shells to python3 -c "import typst; ...").
out_pdf <- "figure_composed.pdf"
.compile_typst("compose.typ", out_pdf, root = getwd())

# Step 4: REMOVE the staged bare copy. Only out_pdf gets harvested.
unlink(c(bare_local, "compose.typ", bare_tmp))

cat("composed:", file.info(out_pdf)$size, "bytes\n")
```

Pin via `make_revision(entity_id=<chain-head-figure-id>, modified_code=<the
R code above>)`. The chain head is normally the entity the user is currently
viewing.

## Common gotchas

- **CLI is not on PATH.** The PyPI `typst` package only exposes the Python
  module (`typst.compile`). There is no standalone `typst` binary you can
  `system2()` directly — always shell to `python3 -c "import typst; ..."`.
- **Image paths are sandboxed by Typst's root.** Absolute filesystem paths
  in `image("...")` fail unless they fall under the `root=` you pass to
  `typst.compile`. The simplest pattern is: write `compose.typ` in the
  same directory as the bare figure, reference the figure by bare
  filename, and pass that directory as `root=`.
- **Fonts**: Typst falls back to a system sans-serif if your named font
  isn't available. `Liberation Sans` / `Liberation Serif` ship with most
  Linux containers including ours.
- **No autonumbering by default.** The templates above use
  `numbering: none` so the caption reads literally "Figure 1." If you
  want Typst's auto-numbering (e.g. for L4 manuscripts), drop the
  `numbering: none` line.
- **PDF inputs work the same as PNG**. Prefer bare figures saved as PDF
  (vector) so the composed PDF stays infinitely zoomable.
- **Width vs page width**. `image(width: 100%)` fills the *content area*
  (page minus margins). To match a specific bare-figure aspect ratio,
  set `#set page(width, height)` first.
- **Long caption breaks**: Typst wraps automatically; don't pre-wrap with
  `\n`s in your source — they become literal line breaks.

## Quick Typst syntax reference

```typst
// Inline formatting
*bold*  _italic_  `code`

// Set page geometry
#set page(width: 7in, height: 5in, margin: 0.3in)
// or paper presets: paper: "a4" | "us-letter"

// Set text style
#set text(font: "Liberation Sans", size: 10pt, weight: "regular")

// Set paragraph style
#set par(justify: false, leading: 0.7em)

// Embed an image
#image("path.pdf", width: 100%)

// Figure with caption (Typst's first-class figure block)
#figure(image("p.pdf"), caption: [*Figure 1.* ...])

// Grid layout (rows auto-flow)
#grid(columns: (1fr, 1fr), gutter: 12pt, [A...], [B...])

// Stack (vertical box)
#stack(dir: ttb, spacing: 6pt, [top], [bottom])

// Absolute positioning inside a containing block
#place(top + left, dx: -4pt, dy: -4pt, [*A*])

// Vertical / horizontal spacers
#v(8pt)   // vertical
#h(8pt)   // horizontal
```

## When NOT to use this recipe

- The user wants a plain data plot (no caption, no panel composition).
  Stay in matplotlib / ggplot. Don't compose.
- The user wants an interactive figure (HTML, Plotly). Typst outputs
  static PDFs only.
- The user wants editable presentation slides. Use a different tool.
- The user is iterating on the *data* (changing what's plotted, not the
  layout). Iterate the bare figure first; compose once the data is stable.
