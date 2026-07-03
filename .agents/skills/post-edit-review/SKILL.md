---
description: After editing track.py (and any accompanying config/processing.yml or notebook consumer cell), review the changes for: section-label comments, stale docstrings, mismatched return types, dead parameters, style consistency, a config version bump when YAML changed, and notebook cell consistency with sibling cells.
---

## Self-review checklist

After making changes to `track.py`, run through:

1. **Comments**: Remove any newly added `# --- section label ---` or `# step name` comments that label *what* the code does. Only keep comments that explain *why* (non-obvious constraints, algorithmic choices, edge cases).

2. **Return types**: Verify `->` annotations match actual return behavior. If the function only does `plt.show(); plt.close(fig)`, the return type should be `-> None`, not `-> plt.Figure | None`. If early exits use `return None`, that's fine.

3. **Docstrings**: All module-level public functions (no `_` prefix) must follow the **canonical template** documented in `docstring-convention.md` (memory file). The template is authoritative; this checklist is a summary. Key checks:

   - **Structure**: 摘要行 → 散文描述（吸收旧 功能: 内容）→ `参数:` → `返回:`（有返回值时）→ `输出:`（写文件时）→ `说明:`（可选，重要注意事项）。节顺序不可颠倒。
   - **`输出:` vs `返回:` for file-writers**: a function that writes files must document the written paths *somewhere*. If it returns `None` and mainly draws/saves → give it an `输出:` section listing each `` `path/pattern` `` (note the `save_fig`-style gate). If it returns a `dict` that already carries the paths (e.g. `output_dir`/`*_path`/`figures` keys) → those live in `返回:` and a separate `输出:` is redundant — do NOT flag such a function for "missing 输出:". Both is allowed when the dict is rich AND explicit path patterns add value (`plot_argo_hotspots`). Some plotters use `'''` not `"""` — scan via AST/`get_docstring`, not `grep '"""'`.
   - **Param entries**: MUST be bulleted `- name (type): 描述。` — no plain paragraphs after `参数:`.
   - **Section indent**: Section headers (`参数:`/`返回:`/`输出:`/`说明:`) at **4 spaces**; list items at **8 spaces**. Sub-headings under `说明:` at 8+ spaces.
   - **Legacy headers** (`功能:`/`模式差异:`/`显示模式:`/`统计口径:`/`流程:`/`数据源:`/`返回契约:`) are BANNED as top-level sections; fold into `说明:` sub-headings.
   - **No `•`**: Use `-` for all list markers.
   - **Bare `<...>`**: Any path/placeholder pattern with angle brackets (e.g. `<kind>`, `<region>`, `<method>`) MUST be wrapped in backticks — `` `<kind>` `` — or they get parsed as HTML and swallow following text.
   - **Signature sync**: Every parameter in the signature must appear in `参数:`; no removed params left behind; no stale default values (e.g. `默认 80` left after the default changed to 60).
   - **Return type**: `->` annotation must match actual behavior. Functions that only `plt.show(); plt.close(fig)` should be `-> None`, not `-> plt.Figure | None`. Conversely, a function that returns a dict but has an early `return None` (e.g. empty-data abort) is `-> dict | None`, not `-> dict` — grep the body for ALL `return` statements, not just the last one.

   Internal / underscore-prefixed helpers may use a brief one-paragraph docstring; the full template only applies to public functions.

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
