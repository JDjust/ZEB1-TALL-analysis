#!/usr/bin/env bash
rg -l "MCM5" /opt/bioinfo/envs/singlecell-py/lib/python3.11/site-packages/scanpy | head
rg -n "s_genes" /opt/bioinfo/envs/singlecell-py/lib/python3.11/site-packages/scanpy/tools/_score_genes.py | head
