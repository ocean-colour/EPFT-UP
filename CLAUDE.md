# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EPFT-UP (Empirical Phytoplankton Functional Types with Uncertainty and Provenance)
derives empirical PFTs from ocean optics and related observations, with rigorous
treatment of uncertainty and full provenance for every product. We will generate the
algorithms, metrics, and diagnostics to share with the community.

Repository layout:

- `claude_prompts/` — prompts and task definitions that drive this work. Read the
  relevant prompt doc (starting with `start_up.md`) before acting, and execute only
  the numbered task you were pointed at.

## Working Conventions

- **Git:** The user (J. Xavier Prochaska) will perform all git commands (add, commit,
  push, etc.). Do not run git commands that change repository state unless explicitly
  asked. Read-only git commands (e.g. `git status`, `git diff`, `git log`) are fine.
- **Calculations:** If you do any calculation, generate it as a Python script and write
  it to disk so that it can be added to the repository. Do not perform one-off
  calculations only in memory or in the chat.
- **Python environment:** If you need to run Python, use the `ocean14` conda environment
  (e.g. `conda run -n ocean14 python script.py`).
- **Logging:** Record completed work under the `## Logs` section of the prompt doc that
  drove it, dated, including what was learned about the repository.

## Related Repositories

- **IOPtics:** Related work lives on this computer at
  `/Users/xavier/Oceanography/python/IOPtics`. It is also the upstream source for this
  project's `.claude/skills/` and `.claude/settings.json`.
