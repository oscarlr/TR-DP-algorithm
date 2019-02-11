#!/bin/bash
set -e -x

python pacMonStr_V1.py \
    test_query.fasta \
    test_prefix.fasta \
    test_suffix.fasta \
    test_tr.fasta
