---
name: claude-memory
description: Use when a task is about the Oxygen repo and needs durable project memory from prior Claude work, or when Codex should read, update, add, or reconcile notes in the Claude project memory directory for Oxygen.
---

# Claude Memory

Use this skill for Oxygen tasks when durable project memory matters: research background, prior experiment history, literature judgments, workflow preferences, environment caveats, or explicit requests to update project memory.

The live memory source for Oxygen is Claude project memory:

- `/home/user3/.claude/projects/-mnt-w2-scratch-user3-Oxygen/memory/`

## Workflow

1. Read `MEMORY.md` in that directory first.
2. Open only the note files relevant to the current task.
3. If the task needs a durable update, edit the relevant Claude memory file directly.
4. If a new durable topic does not fit an existing note, create one new note in that directory and add one index line to `MEMORY.md`.

## Update Rules

- Treat the Claude memory directory as the hot source of truth for Oxygen.
- Prefer updating an existing note over creating a new near-duplicate note.
- Keep note style compatible with the existing memory set: concise title, short summary line in `MEMORY.md`, and Chinese-friendly prose when the surrounding files use it.
- Update only durable project knowledge. Do not store one-off scratch work, temporary execution logs, or ephemeral command output.

## Do Not Do

- Do not maintain a separate live Codex Oxygen memory tree.
- Do not recreate `docs/ai-memory/` as a second hot memory source unless the user explicitly asks for that structure again.
- Do not update Codex top-level memory for Oxygen except to preserve the routing pointer to this Claude memory directory.

## Scope Notes

- `/home/user3/scratch/Oxygen` and `/mnt/w2/scratch/user3/Oxygen` refer to the same repo.
- This skill is Oxygen-specific. Do not reuse this Claude memory directory for other projects.
