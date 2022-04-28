#!/bin/bash

# script to run all the EBV analyses






#---------------------- 0 get sequence data and align it
# download formatted EBV sequences from genbank
python3 scripts/02_accession_list_rename_download.py > data/ebv_genbank_raw.fa

