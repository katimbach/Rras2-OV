#!/bin/bash

eval "$(conda shell.bash hook)"

conda activate ucdenv

python ./OV_UCD_batch.py OVA522

python ./OV_UCD_batch.py OVA547

python ./OV_UCD_batch.py OVA708

python ./OV_UCD_batch.py OVA818
