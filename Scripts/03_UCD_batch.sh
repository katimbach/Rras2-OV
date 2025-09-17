#!/bin/bash

eval "$(conda shell.bash hook)"

conda activate ucdenv

python ./Scripts/03_UCD_batch.py OVA522

python ./Scripts/03_UCD_batch.py OVA547

python ./Scripts/03_UCD_batch.py OVA708

python ./Scripts/03_UCD_batch.py OVA818
