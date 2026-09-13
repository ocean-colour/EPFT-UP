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

### 2026-09-12 (Copy the skills from IOPtics)

Executed prompt #2: copied `.claude/skills/` from the IOPtics repository.

What was done:

- Created `.claude/` in this repository and copied
  `/Users/xavier/Oceanography/python/IOPtics/.claude/skills/` into it verbatim.
- Two skills came across, one `SKILL.md` each:
  - `critical-partner` — constructive disagreement: surface an untested assumption
    before agreeing, argue the strongest opposing case, retreat only on new evidence,
    review weaknesses first, name emotional investment, do not manufacture flaws, end
    with a question rather than a summary.
  - `grill-me` — interview the user one question at a time down each branch of the
    design tree, always offering a recommended answer, and explore the codebase rather
    than asking what the codebase can answer.
- Verified with `diff -r` that the copies are byte-identical to the IOPtics originals.

What I learned about this repository:

- `.claude/` did not exist before this task; it now holds only `skills/`.
  `settings.json` arrives with prompt #3.
- The two skills were copied unmodified on purpose: unlike `settings.json` (whose
  allow-list accretes repo-specific paths), the skill files contain nothing
  IOPtics-specific, so there is nothing to prune. Both are repo-agnostic behavior
  specifications.
- IOPtics' skills are authored as a YAML frontmatter block (`name`, `description`)
  followed by the instruction body. The `description` field is what a session matches
  against when deciding whether to invoke the skill, which is why it is long and
  trigger-phrase heavy in `critical-partner`.

No git commands that change repository state were run.

### 2026-09-12 (Copy the settings.json policy from IOPtics)

Executed prompt #3: created `.claude/settings.json` from the IOPtics policy.

What was done:

- Copied the permission policy from
  `/Users/xavier/Oceanography/python/IOPtics/.claude/settings.json`, keeping the
  allow/deny/ask structure and dropping the accumulated path-specific entries.
- **Kept (37 allow entries):** read-only shell (`ls`, `cat`, `cd`, `pwd`, `echo`,
  `find`, `grep`, `rg`, `head`, `tail`, `wc`, `sort`, `which`, `env`), benign file ops
  (`mkdir`, `touch`, `cp`, `mv`), read-only git (`status`, `diff`, `log`, `show`,
  `branch`), conda/python/pip/pytest/jupyter including `conda run -n ocean14`, and the
  four publisher `WebFetch` domains (Optica, Wiley/ASLO, doi.org, Taylor & Francis).
- **Pruned (7 entries):** the two one-off Crossref `curl` invocations with a literal DOI
  and a literal bibliographic query baked in; the two `python3 -c` JSON-parsing
  one-liners that went with them; two absolute-path `pytest` invocations naming
  machine-specific interpreters (`/home/xavier/miniconda3/...` and
  `/Users/xavier/miniforge3/...`); an import check for the `bing` package; and a `grep`
  of an IOPtics report file. None of these would ever match in this repository.
- **deny** (unchanged): `sudo`, `rm -rf /`, `rm -rf ~`, and `git push` / `git commit` /
  `git reset` / `git rebase`. **ask** (unchanged): `rm:*`.
- Validated the result parses as JSON.

What I learned about this repository:

- `.claude/` now holds both `settings.json` and `skills/`, completing the Claude
  configuration block; prompt #4 (`Basic start up`) is the remaining task.
- The git rule is enforced in two places, not one: prose in `CLAUDE.md` and the deny
  list here. The deny list is the binding one — it blocks `git commit`/`push`/`reset`/
  `rebase` at the permission layer regardless of what a session is asked to do.
- The pruned entries are a good illustration of how these allow-lists decay: every one
  of them was a session-specific permission grant that got written back to the
  committed policy file. Worth re-pruning this file periodically rather than letting it
  accrete, and worth preferring generic patterns (`Bash(pytest:*)`) over the literal
  command a session happened to run.
- The Crossref entries that were pruned are evidence that literature checking is a
  routine activity in these repos. If it becomes routine here, add a *generic* Crossref
  allow pattern rather than re-adding per-query ones.

No git commands that change repository state were run.
