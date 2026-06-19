# aba-recipe-pack — standard ABA bioinformatics cookbook

This repo is a **content overlay** for [ABA](https://github.com/kharchenkolab/aba).
It carries the bioinformatics recipes, capability catalog, and reference
knowhow that most installations layer on top of the platform's system
bundle. The platform repo (`aba`) keeps the runtime + a tiny set of
operating skills; the cookbook lives here so it can evolve independently.

## Layout

```
.
├── LICENSE                  Apache License 2.0
├── NOTICE                   Third-party attributions (see "License" below)
├── CONTEXT.md               Prepended to the system prompt as the cookbook header
├── policies.yaml            Sandbox-preinstall lists, approval gates, runtime knobs
├── catalog/
│   ├── python_bio.yaml          Python capabilities (clean, hand-curated)
│   ├── r_bioconductor.yaml      R/Bioconductor capabilities (clean, hand-curated)
│   └── biomni-derived/          Reference catalogue mined from upstream Biomni
├── knowhow/                 Reference markdown that recipes cross-link to
└── recipes/
    ├── <domain>/                Hand-curated recipes (each domain dir)
    │   ├── *.md                     Flat recipes
    │   └── <folder-skill>/SKILL.md  Folder skills (with optional references/)
    └── biomni-derived/<domain>/ Distilled-from-Biomni recipes, by domain
```

**Why the split.** Hand-curated recipes get rigorous review and live
under `recipes/<domain>/`. Recipes mined from upstream
[Biomni](https://github.com/snap-stanford/Biomni) and reformatted for ABA
live under `recipes/biomni-derived/<domain>/`. Each derived recipe carries
`source: biomni:tool/...` provenance frontmatter. Both branches are
discoverable by the same loader (the agent ranks them together); the
folder split is so humans can tell hand-curated content from automated
derivation at a glance.

Each recipe is a `.md` file with YAML frontmatter (`name`, `description`,
`when_to_use`, `capabilities_needed`, `domain`, `keywords`, `aliases`, `source`,
…) and a procedural body the agent reads via `Skill(skill=…)`. Folder skills
add a `SKILL.md` + sibling `references/` and `scripts/` directories.

## How an institution uses this

Clone alongside the platform repo (recommended):

```bash
cd /path/to/aba/..
git clone https://github.com/kharchenkolab/aba-recipe-pack
```

…or to a known absolute path:

```bash
git clone https://github.com/kharchenkolab/aba-recipe-pack /srv/aba/content/aba-recipe-pack
```

Declare it as a content layer in the ABA deployment config:

```yaml
# /etc/aba/deployment.yaml  (or ~/.aba/deployment.yaml for dev)
layers:
  - name: aba-recipe-pack
    path: /srv/aba/content/aba-recipe-pack
```

ABA's skill loader walks layers in order. An optional *institution* layer
can sit on top of this one to override recipes or add lab-specific ones —
see `misc/content_layers.md` in the platform repo for the layering rules.

## Adding a recipe

Drop a `.md` file in the appropriate domain folder
(`recipes/<domain>/<name>.md`, or `recipes/biomni-derived/<domain>/<name>.md`
if it's mined from Biomni). Frontmatter must include at least `name` and
`description`. After commit + push:

```bash
# on the deployment host:
git -C /srv/aba/content/aba-recipe-pack pull
curl -X POST http://localhost:8000/api/skills/reload
```

## Overriding a recipe at the institution layer

Two patterns, both supported by the same loader:

- **Alias-style (recommended).** Author a NEW recipe with a fresh canonical
  `name:` and list the base name under `aliases:`. The agent's references
  to the base name resolve to your new recipe; the base stays in
  `aba-recipe-pack` unchanged, so `git pull` doesn't fight your overlay.

- **Same-name override (heavy hand).** Author a recipe with the SAME
  canonical `name:` as the base. Loader's last-write-wins (overlay loaded
  after `aba-recipe-pack`) makes your version authoritative. You own
  update tracking against the base — `aba-overlay-status` flags drift.

## License

Apache License, Version 2.0 — see `LICENSE`. Third-party attributions
(notably the upstream Biomni-derived material under
`recipes/biomni-derived/` and `catalog/biomni-derived/`) are listed in
`NOTICE`. Procedural bodies in this repo are written to use ABA's own
libraries and conventions; upstream Biomni's code is NOT imported or
executed at runtime — `source:` frontmatter records provenance only.
