#!/usr/bin/env bash

# Download pretrained models and associated data required for predictions from Zenodo

curl -L -o zenodo_archive.zip "https://zenodo.org/records/5155330/files/Nardus/zoonotic_rank-v1.2.1.zip?download=1"

unzip -jo zenodo_archive.zip "Nardus-zoonotic_rank-25ef305/CalculatedData/*" -d "CalculatedData/"
unzip -jo zenodo_archive.zip "Nardus-zoonotic_rank-25ef305/RunData/PN_LongRun/*" -d "RunData/PN_LongRun/"
unzip -jo zenodo_archive.zip "Nardus-zoonotic_rank-25ef305/RunData/AllGenomeFeatures_LongRun/*" -d "RunData/AllGenomeFeatures_LongRun/"

rm zenodo_archive.zip
