# PHYLUCE Bait Design: parameters and settings

Settings used for the programs and Python script for our bait set design following the online tutorial of PHYLUCE IV:
https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html (Faircloth, 2012; Faircloth et al., 2016)

#### Transcriptome data used:
Transcriptome data used for bait design can be found back in the folder 'transcriptomes_baitdesign'.

This folder contains two subfolders:
* transcriptomes_haplosclerids
* transcriptomes_outgroups

A list with the abbreviations can be found in a separate text document.

## Tutorial IV: Identifying UCE loci and designing baits to target them

### STARTING DIRECTORY STRUCTURE

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
* AMQU = Amphimedon queenslandica
* HCIN = Haliclona (Reniera) cinerea
* HIND = Haliclona (Rhizoniera) indistincta
* HOCU = Haliclona (Haliclona) oculata
* HSIM = Haliclona (Haliclona) simulans
* HTUB = Haliclona (Reniera) tubifera
* HVIS = Haliclona (Rhizoniera) viscosa
* NCOM = Neopetrosia compacta

### CLEANUP THE TRANSCRIPTOMES
Raw data was collected of the species used during bait design, and cleaned, assembled and checked for quality using TransPi: https://github.com/PalMuc/TransPi (Rivera-Vicéns et al. 2021)

#### Put transcriptomes in their own directories

```python
> cd uce-haplosclerida/transcriptomes
> for critter in *; do mkdir ${critter%.*}; mv $critter ${critter%.*}; done
```

#### Convert transcriptome to 2bit format

```python
> for critter in *; do faToTwoBit $critter/$critter.fasta $critter/${critter%.*}.2bit; done
```

#### Simulate reads from transcriptomes

Install 'Art' Package (Huang et al. 2012)

```python
> conda install -c bioconda art

> cd uce-haplosclerida
> mkdir reads
> cd reads
```

```python
> art_illumina --paired --in ../transcriptomes/AMQU/AMQU.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HCIN/HCIN.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HIND/HIND.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HOCU/HOCU.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HSIM/HSIM.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HTUB/HTUB.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/HVIS/HVIS.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
> art_illumina --paired --in ../transcriptomes/NCOM/NCOM.fasta --out AMQU-pe100-reads --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na
```

#### Merging reads together
```python
> for critter in AMQU HCIN HIND HOCU HSIM HTUB HVIS NCOM;
do
	echo "working on $critter";
	touch $critter-pe100-reads.fq;
	cat $critter-pe100-reads1.fq > $critter-pe100-reads.fq;
	cat $critter-pe100-reads2.fq >> $critter-pe100-reads.fq;
	rm $critter-pe100-reads1.fq;
	rm $critter-pe100-reads2.fq;
	gzip $critter-pe100-reads.fq;
done;
```

### READ ALIGNMENT TO THE BASE TRANSCRIPTOME

Prepare the base transcriptome

```python
> cd uce-haplosclerida
> mkdir base
> cd base
```

Copy base transcriptome to base-directory

```python
> cp -p ../transcriptomes/AMQU/AMQU.fasta ./
```

## References
Faircloth BC. 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32:786-788. doi:10.1093/bioinformatics/btv646.

Faircloth BC, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi:10.1093/sysbio/SYS004.

Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics, 28(4), 593–594.

Rivera-Vicéns RE, García-Escudero CA, Conci N, Eitel M, Wörheide G. 2021. TransPi–a comprehensive TRanscriptome ANalysiS PIpeline for de novo transcriptome assembly. bioRxiv 2021.02.18.431773; doi:10.1101/2021.02.18.431773








