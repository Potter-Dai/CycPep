# Accurate Structure Prediction for Cyclic Peptides Containing Proline Residues with High-Temperature Molecular Dynamics

This repository describes a simple protocol which can be used to run High-T MD simulations, perform probability density reweighting, and do clustering analysis in **Accurate Structure Prediction for Cyclic Peptides Containing Proline Residues with High-Temperature Molecular Dynamics**.

Complete documentation for using DPA to cluster all conformations is available on the [Read The Docs](https://dpaclustering.readthedocs.io/en/latest/) website.

## Examples
#### Run MD simulations at 700 K with RSFF2C force field

```bash
sh run_700K.sh
```

#### Calculate the probabilities of grouped clusters before and after reiweighting

```
python Proba_now_w.py
```

#### Calculate the top RMSD of grouped clusters' centroids

```
python Toprmsd_now_w.py
```
