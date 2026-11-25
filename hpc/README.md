# Working on the HPC

One-time setup (pulls the Bioconductor devel image, R 4.6.x):

```bash
git clone git@github.com:sinar9776/scdRs.git   # or https://...
cd scdRs
bash hpc/setup_container.sh
```

Each time Claude pushes changes:

```bash
git pull
bash hpc/check.sh test     # fast: unit tests only
bash hpc/check.sh          # full: R CMD check + BiocCheck
```

Paste the summary at the end of the log (or the whole `hpc/logs/*.log`)
back to Claude.

Override paths with `IMG=...` and `SCDRS_RLIB=...` if needed.
