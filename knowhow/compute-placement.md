---
name: compute-placement
description: Decision guide for WHERE a compute step runs — the interactive kernel (local), a background job, or a Slurm job with a resource request. Covers the local-vs-Slurm tradeoff (cores / mem / GPU / walltime vs queue wait), how to size a request, and the fresh-process caveat. Use when planning a heavy / parallel / GPU / long-running step, especially on a cluster.
when_to_use: You are about to run a compute-heavy / parallel / GPU / long step and must decide interactive vs background vs Slurm; or the per-turn "Compute environment:" line / describe_compute shows a Slurm cluster and you're weighing local vs Slurm; or a step needs more cores/mem/GPU than this node has, or might exceed the session's remaining walltime.
avoid_when: A short interactive step on a local (non-cluster) deployment — just run it with run_python/run_r, no placement decision needed. Not for the mechanics of a specific analysis (those are the recipes).
invocation: interactive
kind: knowhow_draft
requires_tools: [run_python, run_r]
keywords: [compute, Slurm, HPC, background job, sbatch, cores, memory, GPU, walltime, partition, queue wait, QOS, parallel, estimated_runtime_min, est_cores, est_mem_gb, est_gpu, describe_compute, interactive kernel, where to run, placement]
domain: meta
source: "ABA execution router (core/exec/router.py) + compute_env/describe_compute; the environment-aware routing redesign (local = explicit-background only; Slurm = resource/walltime driven)."
audience: agent
---

# Compute placement — interactive, background, or Slurm?

The question: a step is about to run — should it run in the **interactive kernel**
(state persists; you see the result now), as a **background job** (async; you get
a `job_id` and poll), or — on a cluster — as a **Slurm job** with a resource
request? The wrong call either wastes a node, gets killed by the allocation
walltime, or (the classic trap) silently relocates a state-dependent cell into a
fresh process where `panel`/`con`/your loaded matrices no longer exist.

First read the standing **"Compute environment:"** line in your context (it rides
every turn) — it tells you the mode, this node's cores/mem/GPU, the remaining
walltime, and on a cluster the live partitions. For full detail (all partitions,
your access, current load) call **`describe_compute`**.

## The one rule that prevents the worst bug

**A background or Slurm job is a FRESH PROCESS.** It has NONE of your interactive
kernel's in-memory objects. Always: load inputs from disk inside the job, write
outputs to disk; never reference a variable created in an earlier interactive
cell. (The failure: a backgrounded preprocessing step → `object 'panel' not
found`.)

## Local deployment (mode = local)

There is no scheduler to escape to, so:

- **Run interactively.** Kernel state persists across cells — that's the point.
- A **long** cell is fine: raise `timeout_s` (e.g. `run_python(code, timeout_s=1800)`).
  Do NOT background a cell just to "avoid a timeout."
- Use **`background=True`** ONLY when: the **user asks** for it, or you want to
  **fan out several independent jobs in parallel** (each self-contained). Nothing
  is auto-backgrounded on a runtime guess.

## Slurm deployment (mode = slurm)

Most steps still run interactively on the node. Send a step to Slurm
(`run_python`/`run_r` with `background=True` + the `est_*` sizing args) when:

| Trigger | Why |
|---|---|
| Needs **more cores** than this node has, AND a partition offers bigger nodes | parallel speedup you can't get locally |
| Needs a **GPU** this node lacks (and a GPU partition exists) | only Slurm has the device |
| Needs **more memory** than this node has | won't fit locally |
| Might **approach/exceed the session's remaining walltime** | a long interactive cell dies when the allocation ends; a Slurm job gets its own walltime |
| The **user asks**, or you're fanning out independent jobs | explicit |

If a step fits comfortably on this node and is short, just run it interactively —
even on a cluster. Don't pay a queue wait for a 2-minute job.

### Weighing speed vs queue wait (the 64-vs-4 case)

The router auto-routes the *hard* cases (won't fit / would be killed). The
**speed** case is YOUR judgment, because it depends on the computation:

1. Estimate the time locally (this node's cores) vs on Slurm (a partition's
   cores/GPU), and read the partition's **`wait`** label from `describe_compute`.
2. Go to Slurm when `local_time − slurm_time` comfortably beats the wait. Example:
   **64 cores → ~30 min + ~10 min queue wait** beats **4 cores → ~8 h locally**.
   For a 2-minute job, the wait never pays off → run local.
3. Express the choice by setting `est_cores` (and `est_gpu`/`est_mem_gb`) to what
   the job *wants*. Requesting more cores than this node has routes it to Slurm
   and sizes the `sbatch`.

## Sizing a Slurm request

Pass the estimate on `run_python`/`run_r`:

- `est_cores` — cores the job parallelizes over (set thread counts / `n.cores` to match).
- `est_mem_gb` — peak memory (be generous; an OOM kill wastes the whole queue wait).
- `est_gpu` — `true` if it needs a GPU.
- `estimated_runtime_min` — wall-clock estimate (drives the walltime request + the
  walltime guardrail; it is NOT a local-background trigger anymore).

The deployment maps these to a concrete partition/QOS/walltime (clamped to the
partition limits shown by `describe_compute`). Keep the request within a
partition's `cpus_per_node` / `mem_gb_per_node` / `max_walltime`.

## Checklist before you background / submit to Slurm

- [ ] The job loads its inputs **from disk** and writes outputs **to disk** (fresh process).
- [ ] On Slurm: `est_cores`/`est_mem_gb`/`est_gpu`/`estimated_runtime_min` set, within a partition's limits.
- [ ] The local node genuinely can't (or shouldn't) run it, OR the user asked, OR you're parallelizing.
- [ ] You'll `get_job_status(job_id)` to report progress, not deflect to the UI.

## See also

- `describe_compute` — the live landscape (call it when planning a heavy step).
- Executable recipes for the analysis itself (this guide is only about *where* it runs).
