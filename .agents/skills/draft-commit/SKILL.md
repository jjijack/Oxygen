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

**Body** — 1–2 sentences (present tense) saying what the feature does. Stop
there. The mechanism, the function names involved, the per-mode plumbing — all
of that is in the diff; only add a clause if the subject alone is genuinely
unclear. When in doubt, cut. No bullets, no blank-line sections, no headers.
Backtick function/param names. A one-sentence body is the target, not the floor.

### Omit — never mention
Incidental bug fixes, dead-code removal, docstring/comment edits, config
(`processing.yml`) tweaks, notebook changes — unless one of them *is* the
primary change.

### Good vs. too-complex

✅ House style (real commit):
> Introduce the local winter-MLD ventilation threshold
>
> Replace the fixed 300 m heave cutoff for the ventilated/isolated split with a
> per-point local winter mixed-layer depth, selectable via `ventilation_mode`.

❌ Avoid: multi-section or bulleted bodies, `Changes:` / `Details:` headers,
enumerating every new function and parameter, restating the diff file-by-file,
spelling out the mechanism or downstream wiring, or noting the small bug you
fixed in passing. If the body runs past two sentences, it's almost certainly
too long — cut it back to one feature in one sentence.

### Procedure
1. `git diff --cached --stat` + `git diff --stat` — see scope.
2. `git diff --cached` + `git diff` — read the actual change (focus `track.py`).
3. `git log --format="%s%n%b" -5` — match current phrasing.
4. Pick the ONE headline feature; draft a short subject + 1–2 sentence body.
5. **Print** the message for review. **Never** run `git commit` — the user commits manually.
