# Working on the HPC

These scripts live on the `tools` branch only; the default branch
(`devel`) holds just the package, as Bioconductor requires. The folder
`hpc/` is ignored by git on `devel`, so the scripts can sit in a clone of
`devel` without being committed.

One-time setup (pulls the Bioconductor devel image, R 4.6.x):

```bash
git clone https://github.com/sinar9776/scdRs.git
cd scdRs
git fetch origin tools
git archive origin/tools hpc | tar -x     # puts hpc/ into the clone
bash hpc/setup_container.sh
```

Each time the package changes:

```bash
git pull
git fetch origin tools && git archive origin/tools hpc | tar -x
bash hpc/check.sh test     # fast: unit tests only
bash hpc/check.sh          # full: R CMD check + BiocCheck
```

Paste the summary at the end of the log (`hpc/logs/*.log`) back to Claude.

Override paths with `IMG=...` and `SCDRS_RLIB=...` if needed.
