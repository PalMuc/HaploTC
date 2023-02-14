# PHYLUCE Bait Design: parameters and settings

Settings used for the programs and Python script for our bait set design following the online tutorial of PHYLUCE IV:
https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html

#### References
Faircloth BC. 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32:786-788. doi:10.1093/bioinformatics/btv646.

Faircloth BC, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi:10.1093/sysbio/SYS004.

#### Transcriptome data used:
Transcriptome data used for bait design can be found back in the folder 'transcriptomes_baitdesign'.

This folder contains two subfolders:
* transcriptomes_haplosclerids
* transcriptomes_outgroups

A list with the abbreviations can be found in a separate text document.

## Tutorial IV: Identifying UCE loci and designing baits to target them

#### Starting directory structure

```python
> mkdir uce-haplosclerida
> cd uce-haplosclerida
```

### DATA AND PREPARATION

```python
> mkdir transcriptomes
> cd transcriptomes
```

The 'transcriptomes' folder contains all the transcriptomes necessary for a particular bait set (see Table S3 for an overview). For simplicity, we outline bait set version 50 (haplo-only, with Amphimedon queenslandica (AMQU) as base transcriptome) here.

#### Haplosclerid transcriptomes:

AMQU = Amphimedon queenslandica
HCIN = Haliclona (Reniera) cinerea
HIND = Haliclona (Rhizoniera) indistincta
HOCU = Haliclona (Haliclona) oculata
HSIM = Haliclona (Haliclona) simulans
HTUB = Haliclona (Reniera) tubifera
HVIS = Haliclona (Rhizoniera) viscosa
NCOM = Neopetrosia compacta