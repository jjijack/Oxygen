---
name: oxygen-background-runner
description: Use when launching, handing off, or debugging long-running background Python jobs in the Oxygen repo, especially Dask/plot workflows that must survive Codex tool sessions, write real-time logs, avoid chatty polling, and optionally finish under a separate background agent.
---

# Oxygen Background Runner

Use this skill for Oxygen background jobs such as hotspot sweeps, detector sensitivity runs, and other long Dask/plot workflows. The goal is to start exactly one durable job, with real-time logs, using the `plot` environment correctly, while avoiding wasteful polling in the main rollout.

## Default Policy

Do not keep the main rollout alive just to watch a long-running job.

Choose between three modes:

1. **Short foreground run**: for quick checks, bounded validations, or tasks where the next decision depends immediately on the result.
2. **Detached launch and exit**: default for long Oxygen jobs that mainly write outputs/logs and do not need interactive steering.
3. **Detached launch plus background reminder**: only when the user explicitly wants a completion notification after the current turn ends.

Use these heuristics:

- Prefer **short foreground run** when the command is expected to finish quickly and the result is needed right now.
- Prefer **detached launch and exit** when the task is expected to run for minutes to hours, especially Dask sweeps, multi-year exports, and plot pipelines.
- Prefer **detached launch plus background reminder** only when the user explicitly asks for a callback-style workflow. In Codex, this should be handed to a separate background agent when that capability is available.

## No-Polling Rule

Avoid tight polling loops such as repeated `ps`, `sed`, `tail`, or short `write_stdin` waits.

- If a foreground run is appropriate, use one blocking command with a generous wait rather than many short polls.
- If a job should outlive the current turn, launch it detached, verify startup once, report the log path, and stop.
- Do not keep checking the same log unless there is a concrete debugging reason.
- When a completion reminder is required, move the waiting responsibility to a background agent instead of the main rollout.

## Default Pattern

Prefer direct execution with the `plot` environment Python plus explicit environment variables. Avoid `conda run` for detached background jobs because it can buffer output or fail to leave a durable child process.

Use this shape, replacing the script and arguments:

```bash
zsh -lc 'mkdir -p plot_outputs/test/logs; log="plot_outputs/test/logs/<job>_$(date +%Y%m%d_%H%M%S).log"; launch="plot_outputs/test/logs/<job>_launch.out"; env CONDA_PREFIX=/home/user3/.conda/envs/plot PATH=/home/user3/.conda/envs/plot/bin:$PATH PROJ_LIB=/home/user3/.conda/envs/plot/share/proj GDAL_DATA=/home/user3/.conda/envs/plot/share/gdal MPLCONFIGDIR=/tmp/matplotlib-user3 PYTHONUNBUFFERED=1 setsid -f /home/user3/.conda/envs/plot/bin/python -u <script.py> <args> --log-file "$log" > "$launch" 2>&1 < /dev/null; printf "%s\n%s\n" "$log" "$launch"'
```

Run detached starts with `sandbox_permissions="require_escalated"` because Codex sandboxed processes can be killed with the tool session.

## Preferred Long-Task Workflow

For most Oxygen production runs, the default should be:

1. Write or update a small runner script if needed.
2. Launch the job detached.
3. Verify startup once with `ps` and the first part of the logs.
4. Report the launch details to the user.
5. End the turn.

Do not keep the conversation active only to monitor progress. The user can come back later and ask for a status check, or the task can be delegated to a background agent if they asked for a completion reminder.

## Background Reminder Mode

When the user explicitly wants a Claude-style "submit and forget" workflow with a later reminder:

- Prefer spawning a **separate worker/background agent** to own the long task, rather than keeping the main rollout active.
- That background agent should either:
  - launch the detached OS job, then wait for completion with a sparse, blocking strategy, or
  - own the whole task in its own worktree/context if code edits and the long run belong together.
- The main rollout should hand off the task, confirm the handoff, and then end the turn.

Important:

- The background agent should not busy-poll either. Use long waits and one-shot log reads.
- The main rollout should not also monitor the same job.
- Use this mode only when the user asked for it; otherwise detached launch and exit is the default.

## Script Requirements

For repeatable jobs, write a small runner script under `plot_outputs/test/` unless the workflow is formal enough to belong in `track.py`.

The runner should:

- Import the repo root before importing `track`.
- Call `track.switch_region(...)` explicitly.
- Support `--workers`, `--memory-limit`, and `--log-file`.
- Skip existing outputs by default, with `--overwrite` only when intentional.
- Use `print(..., flush=True)` for progress messages.
- If `--log-file` is provided, tee both `stdout` and `stderr` to that file with line buffering.
- Pass `dask_scheduler="processes"` and a deliberate `dask_workers` value to `track` workflows.

## Validation After Launch

Immediately after launching, verify both process state and logs:

```bash
ps -eo pid,ppid,pcpu,pmem,stat,etime,cmd | rg '<script-name>|<job-log-stem>'
sed -n '1,120p' plot_outputs/test/logs/<job>_YYYYMMDD_HHMMSS.log
sed -n '1,80p' plot_outputs/test/logs/<job>_launch.out
```

Expected signs:

- Main Python process has `PPID 1` after `setsid -f`.
- Dask worker child processes appear under that main process after the workflow starts.
- The task log contains the region switch, settings, and current `[run]` line.
- Tqdm/progress output is visible in the log.
- A Dask dashboard port warning is acceptable if another job already uses port `8787`.

After this one-time validation, stop watching unless you are actively debugging startup failure.

Problem signs:

- No task log file appears: the process likely never reached Python.
- Only an empty launch log appears: the detached command failed before script execution.
- `PROJ: proj_create_from_database` appears: restart with `PROJ_LIB` and `GDAL_DATA` set as above.
- The process is visible only inside sandboxed `ps`: recheck with escalated `ps`.

## Worker Policy

Treat the Oxygen server as a work machine. When the job is mostly CPU-bound and reads local cached/parquet inputs, it is acceptable to use most or all available cores. Check the machine size with `nproc`, then choose `--workers` close to the available core count unless another heavy job is already running.

Be conservative when the job is I/O-bound, especially workflows reading GLORYS NetCDF over sshfs or repeatedly loading remote-mounted data. In those cases, too many workers can overload I/O and make the job slower or less stable.

For mixed workflows such as Argo hotspot sweeps, scale up when CPU is idle and memory is comfortable, but watch process CPU and output progress after launch. Multiple parameter combinations inside one runner usually execute sequentially unless the script explicitly parallelizes combinations, so tune workers for the active combination, not for the total number of queued combinations.

## Monitoring Policy

Use the lightest monitoring that still answers the question:

- **Launch verification**: one `ps` check and one log read.
- **Later status check on request**: one fresh `ps` check and one fresh `tail` or `sed`.
- **Completion reminder mode**: background agent owns the waiting and sends one completion update.

Avoid "heartbeat" monitoring from the main rollout.

## Codex-Specific Notes

This skill is Codex-specific. It does not need a mirrored `.claude` peer.

When deciding whether to use this skill by default:

- Do **not** use it for every task.
- Do use it by default for clearly long, non-interactive Oxygen jobs.
- Before launching, explicitly assess whether the task is:
  - short and decision-blocking,
  - long but fire-and-forget,
  - long and requiring a completion reminder.

That assessment should happen up front so Codex does not drift into token-heavy monitoring by habit.

## Cleanup

If a bad duplicate job is launched, stop only the bad job's main process after checking its command and parent/child tree:

```bash
kill <main-pid>
```

Then recheck `ps` for children. Never kill unrelated worker processes by broad patterns; identify the parent process first.
