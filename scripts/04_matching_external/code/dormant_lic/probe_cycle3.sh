#!/usr/bin/env bash
grep -R -l "MCM5" /opt/bioinfo/envs/singlecell-py/lib/python3.11/site-packages/scanpy | head
grep -n "def score_genes_cell_cycle" -A 40 /opt/bioinfo/envs/singlecell-py/lib/python3.11/site-packages/scanpy/tools/_score_genes.py | head -50
