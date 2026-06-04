---
description: After editing track.py (and any accompanying config/processing.yml or notebook consumer cell), review the changes for: section-label comments, stale docstrings, mismatched return types, dead parameters, style consistency, a config version bump when YAML changed, and notebook cell consistency with sibling cells.
---

## Self-review checklist

After making changes to `track.py`, run through:

1. **Comments**: Remove any newly added `# --- section label ---` or `# step name` comments that label *what* the code does. Only keep comments that explain *why* (non-obvious constraints, algorithmic choices, edge cases).

2. **Return types**: Verify `->` annotations match actual return behavior. If the function only does `plt.show(); plt.close(fig)`, the return type should be `-> None`, not `-> plt.Figure | None`. If early exits use `return None`, that's fine.

3. **Docstrings**: Public / externally-called entry points (e.g. the `plot_*` functions) use the full **说明 → 参数 → 返回** format — a description paragraph, then a `参数:` block documenting *every* parameter, then a `返回:` block if it returns a value. Internal / underscore-prefixed helpers may use a brief one-paragraph docstring. **Don't leave a public entry point with only brief prose** — match its siblings (e.g. `plot_regional_vertical_argo_overview`). Also verify docstrings track the current signature: no removed params still documented, no added params undocumented, no stale default values (e.g. a `默认 80` left after the default changed to 60).

4. **Dead code**: Remove leftover parameters, variables, or branches from earlier iterations of the change.

5. **Style consistency**: Match the patterns of surrounding code — docstring tier per #3, Chinese period `。` in Chinese text, multi-line `float(_CFG.get(...))` for config-driven module constants, etc.

6. **Config version bump**: If this change touched `config/processing.yml` (new key, new sub-section, changed default), record it in the header changelog matching the existing entries. **But first check `git diff`/`git log` whether the top `# version X.Y` entry is already committed**: if that latest version is still uncommitted (only staged/working-tree), fold the new change INTO that pending entry — extend its description — rather than minting a fresh version. Only bump to a new `# version` number when the current top entry is already committed (a released version is frozen). A YAML change with no changelog touch at all is incomplete. Don't put version numbers or line counts in `CLAUDE.md` (they go stale); point to the YAML header changelog instead.

7. **Notebook cell consistency**: If a new function got a consumer cell in a notebook (`GLORYS*.ipynb`), check its 话风/style against the sibling cells that call related functions (e.g. the other `plot_hotspot_anomaly_*` cells). Conventions to match:
   - **Markdown headers are pure English** and bare — no Chinese, no `（parenthetical）`, no crammed-in description. Title cell holds *only* the title.
   - **Descriptions go in their own separate markdown cell** (the prose may be Chinese), not in the header cell — see the `Argo Data Flaws` section for the pattern.
   - **Code cells stay thin**: comment density 2–3 lines (not a wall); usually just call the function with `return_details=False` — don't assign `res =` or add a trailing display expression unless the siblings do; inline `#` hints for key optional args are fine.
   The notebook is a thin consumer — keep cells uniform with their neighbors, no logic.

## Verification

```bash
python3 -c "import py_compile; py_compile.compile('track.py', doraise=True); print('Syntax OK')"
```

```bash
git diff --stat
```

Report issues found with specific line references.
