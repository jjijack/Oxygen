---
description: Draft a commit message following the project's style, print for review, but never commit automatically.
---

## draft-commit

Draft a commit message for the current changes. The user commits **one big
feature at a time**, usually scoped to `track.py`.

### Mental model

One commit = one headline feature. The message *introduces that feature* — it
is **not** a changelog of the diff. The diff almost always also contains small
fixes and refactors made along the way; those do **not** go in the message.
Find the single headline, describe it, stop.

### Style (from project history)

**Subject** — English, imperative, short. Name the main function or concept;
backtick function names.
- `Introduce the \`func_name\` function`
- `Update the \`func_name\` function`
- `Rework the \`func_name\` X selection`
- `Add X and Y aggregated T-S diagram functions`
- `Introduce the <concept> and <concept>` — when the feature spans several functions

**Body** — 2–4 lines of flowing sentences (present tense). Briefly describe
what the feature does and its key mechanism, as if introducing it to a reader.
No bullets, no blank-line sections, no headers. Backtick function/param names.

### Omit — never mention
Incidental bug fixes, dead-code removal, docstring/comment edits, config
(`processing.yml`) tweaks, notebook changes — unless one of them *is* the
primary change.

### Good vs. too-complex

✅ House style (real commit):
> Introduce the `compute_spiciness_anomaly` function
>
> Compute the kernel-weighted spiciness percentile in σ₀ space, returning
> signed δπ and percentile for water-mass deviation. Integrate into T-S
> diagram titles and the hotspot overview pipeline via `annotate_spice`.

❌ Avoid: multi-section or bulleted bodies, `Changes:` / `Details:` headers,
enumerating every new function and parameter, restating the diff file-by-file,
or noting the small bug you fixed in passing. If it reads like a release-notes
page, it's wrong — cut it back to one feature in 2–4 lines.

### Procedure
1. `git diff --cached --stat` + `git diff --stat` — see scope.
2. `git diff --cached` + `git diff` — read the actual change (focus `track.py`).
3. `git log --format="%s%n%b" -5` — match current phrasing.
4. Pick the ONE headline feature; draft a short subject + 2–4 line body.
5. **Print** the message for review. **Never** run `git commit` — the user commits manually.
