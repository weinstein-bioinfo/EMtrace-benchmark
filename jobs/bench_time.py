#!/usr/bin/env python3
"""bench_time.py — run a command, append one timing row to a CSV.

Shell- and OS-agnostic replacement for the shell `time` keyword. Child
stdout/stderr are inherited, so the caller's `2>&1 | tee` and the tool's own
--logs keep working unchanged.

Usage:
    python3 bench_time.py --csv FILE --label NAME -- CMD [ARGS...]

CSV columns: sample,user,system,cpu,total,max_rss_mb,exit_status
  user/system = CPU-seconds (incl. subprocesses), cpu = % (may exceed 100% with
  multiple threads), total = wall-clock seconds, max_rss_mb = peak RSS of the
  largest process in the tree, exit_status = child return code.
"""
import argparse, csv, os, resource, subprocess, sys, time

COLUMNS = ["sample", "user", "system", "cpu", "total", "max_rss_mb", "exit_status"]


def main():
    argv = sys.argv[1:]
    if "--" not in argv:
        sys.exit("usage: bench_time.py --csv FILE --label NAME -- CMD [ARGS...]")
    i = argv.index("--")
    opts, cmd = argv[:i], argv[i + 1:]
    if not cmd:
        sys.exit("bench_time.py: no command after --")

    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True)
    p.add_argument("--label", required=True)
    a = p.parse_args(opts)

    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    rc = subprocess.run(cmd).returncode          # stdout/stderr inherited
    wall = time.perf_counter() - t0
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)

    user = r1.ru_utime - r0.ru_utime
    system = r1.ru_stime - r0.ru_stime
    cpu = round((user + system) / wall * 100) if wall > 0 else 0
    maxrss = r1.ru_maxrss                          # bytes (macOS) / KB (Linux)
    max_rss_mb = maxrss / 1024 / 1024 if sys.platform == "darwin" else maxrss / 1024

    d = os.path.dirname(a.csv)
    if d:
        os.makedirs(d, exist_ok=True)
    new = (not os.path.exists(a.csv)) or os.path.getsize(a.csv) == 0
    with open(a.csv, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(COLUMNS)
        w.writerow([a.label, f"{user:.3f}", f"{system:.3f}", cpu,
                    f"{wall:.3f}", f"{max_rss_mb:.1f}", rc])

    sys.exit(rc)


if __name__ == "__main__":
    main()
