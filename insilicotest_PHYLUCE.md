# *In-Silico* Test Settings

Settings used for the programs and python script used for the MSDA in-silico test following the online tutorial (I) of PHYLUCE: https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html#tutorial-i-uce-phylogenomics (Faircloth, 2012; Faircloth et al., 2016)

Transcriptome data used for the *in-silico* test can be found back in the folder 'data_insilico_test'.

## Tutorial I: UCE Phylogenomics

### STARTING DIRECTORY STRUCTURE

Open conda environment:

```python
> conda activate phyluce-1.7.1
```

```python
> mkdir in-silico-test
> cd in-silico-test
```

* In this in-silico-test folder we also load the bait-set (Supplementary Information III): 'uce-20k-probes.fasta'
* We create a separate folder called 'fastas'. This is where we store all the assemblies. The assemblies are stored in the data_insilico_test/fastas.

```python
> mkdir fastas
```´

### FINDING UCE LOCI

```python
> phyluce_assembly_match_contigs_to_probes \
  --contigs fastas \
  --probes uce-20k-probes.fasta \
  --output uce-search-results \
	--min-identity 85 \
	--min-coverage 85
```´

The output of this file can be found back under 'data_insilico_test/phyluce_assembly_match_contigs_to_probes.log'

### EXTRACTING UCE LOCI

Create a configuration file that lists all the transcriptomes used for the in-silico test: 'taxon-set.conf'

```python
[all]
AAPS
AMQU
CELE
CVAR
etc.
```´

```python
> cd in-silico-test
> mkdir -p taxon-sets/all

> phyluce_assembly_get_match_counts \
  --locus-db uce-search-results/probe.matches.sqlite \
  --taxon-list-config taxon-set.conf \
  --taxon-group 'all' \
  --incomplete-matrix \
  --output taxon-sets/all/all-taxa-incomplete.conf

> cd taxon-sets/all
> mkdir log

> phyluce_assembly_get_fastas_from_match_counts \
  --contigs ../../fastas \
  --locus-db ../../uce-search-results/probe.matches.sqlite \
  --match-count-output all-taxa-incomplete.conf \
  --output all-taxa-incomplete.fasta \
  --incomplete-matrix all-taxa-incomplete.incomplete \
  --log-path log
```´











