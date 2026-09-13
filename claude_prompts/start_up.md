# Getting started 

## Goals

This repository will be used to derive Empirical Phytoplankton Functional Types (PFTs)
from ocean optics and related observations, with rigorous treatment of uncertainty and
full provenance for every product.  We will generate the algorithms, metrics, and
diagnostics to share with the community.

## Prompts

1. Read this file.  Execute the 1st task under "Claude/CLAUDE.md file"
2. Read this file.  Execute the 1st task under "Claude/Skills"
3. Read this file.  Execute the 1st task under "Claude/Settings"
4. Read this file.  Execute the 1st task under "Basic start up"


## Claude

### CLAUDE.md file

1. Please generate a basic CLAUDE.md file for this project.  Have it indicate:

- I will perform git commands
- Add to the CLAUDE.md file:  If you do any calculation, generate it as a python script and write it to disk so that I can add it to the Repository.
- Add to the CLAUDE.md file:  If you need to run Python, use the "ocean14" conda environment.
- Add to the CLAUDE.md file:  Related work lives in /Users/xavier/Oceanography/python/IOPtics

### Skills

1. Copy over the skills/ files from the `IOPtics` repository.

### Settings

1. Copy over the settings.json file from the `IOPtics` repository.  Copy the policy
   (the allow/deny/ask structure), not the accumulated path-specific entries that
   reference IOPtics' own scratchpad and test files.

## Basic start up

1. Generate the basic files that one needs for a Python GitHub repository, e.g. a file for dependencies.  Examine the other Repositories in Oceanography/python to see how I tend to organize things.  Also, make a suggestion of which of their settings.json files to copy into this one.


## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-12 (Generate the root CLAUDE.md)

Executed prompt #1: generated `CLAUDE.md` at the repository root.

What was done:

- Wrote `CLAUDE.md` following the house form used by the sibling repos, closest to
  `cugn-climatology/CLAUDE.md` (which is the only other one carrying all four of the
  requested directives, including a related-repository pointer).
- Sections: Project Overview (with repository layout), Working Conventions
  (git / calculations / python environment / logging), Related Repositories.
- All four directives from the task are present: the user performs git commands
  (read-only git permitted), calculations are written to disk as Python scripts,
  Python runs in the `ocean14` conda environment, and IOPtics at
  `/Users/xavier/Oceanography/python/IOPtics` is named as related work.

What I learned about this repository:

- EPFT-UP is essentially empty at this point: `.gitignore`, `LICENSE`, a two-line
  `README.md`, and `claude_prompts/`. There is no Python package directory, no
  `requirements.txt`/`setup.py`, and no `.claude/` yet — those arrive with prompts
  #2-#4. The layout section of `CLAUDE.md` therefore lists only `claude_prompts/`;
  it should be extended once the package and supporting files exist.
- Sibling `CLAUDE.md` files vary in wording but agree on substance. IOPtics states the
  git rule most explicitly (naming which commands are forbidden vs. permitted), and
  that phrasing was adopted here.
- IOPtics is upstream for `.claude/skills/` and `.claude/settings.json`; that was
  recorded in the Related Repositories section so a future session does not have to
  rediscover it.

No git commands that change repository state were run.
