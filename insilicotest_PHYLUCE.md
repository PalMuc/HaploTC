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
```

### FINDING UCE LOCI

```python
> phyluce_assembly_match_contigs_to_probes \
  --contigs fastas \
  --probes uce-20k-probes.fasta \
  --output uce-search-results \
	--min-identity 85 \
	--min-coverage 85
```

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
```

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
```

### GET STATS ON UCE ASSEMBLIES 

```python
> phyluce_assembly_explode_get_fastas_file \
  --input all-taxa-incomplete.fasta \
  --output exploded-fastas \
  --by-taxon

> for i in exploded-fastas/*.fasta;
do
 phyluce_assembly_get_fasta_lengths --input $i --csv;
done
```

### ALIGNING UCE LOCI

Edge-trimming

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_seqcap_align \
  --input all-taxa-incomplete.fasta \
  --output mafft-nexus-edge-trimmed \
  --taxa 39 \
  --aligner mafft \
  --cores 4 \
  --incomplete-matrix \
  --log-path log
```

Summary stats

```python
> phyluce_align_get_align_summary_data \
  --alignments mafft-nexus-edge-trimmed \
  --cores 4 \
  --log-path log
```

Internal-trimming

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_seqcap_align \
  --input all-taxa-incomplete.fasta \
  --output mafft-nexus-internal-trimmed \
  --taxa 39 \
  --aligner mafft \
  --cores 4 \
  --incomplete-matrix \
  --output-format fasta \
  --no-trim \
  --log-path log
```

GBLOCKS to trimthe loci

```python
> phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
  --alignments mafft-nexus-internal-trimmed \
  --output mafft-nexus-internal-trimmed-gblocks \
  --cores 4 \
  --log log
```

Summary stats

```python
> phyluce_align_get_align_summary_data \
  --alignments mafft-nexus-internal-trimmed-gblocks \
  --cores 4 \
  --log-path log
```

### ALIGNMENT TRIMMING

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_remove_locus_name_from_files \
  --alignments mafft-nexus-internal-trimmed-gblocks \
  --output mafft-nexus-internal-trimmed-gblocks-clean \
  --cores 4 \
  --log-path log
```

### FINAL DATA MATRICES

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_get_only_loci_with_min_taxa \
  --alignments mafft-nexus-internal-trimmed-gblocks-clean \
  --taxa 39 \
  --percent 0.35 \
  --output mafft-nexus-internal-trimmed-gblocks-clean-35p \
  --cores 4 \
  --log-path log
```

### PREPARATION DATA FOR DOWNSTREAM ANALYSIS

Nexus file

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_concatenate_alignments \
  --alignments mafft-nexus-internal-trimmed-gblocks-clean-35p \
  --output mafft-nexus-internal-trimmed-gblocks-clean-35p-raxml \
  --nexus \
  --log-path log
```

Phylip file

```python
> cd in-silico-test/taxon-sets/all

> phyluce_align_concatenate_alignments \
  --alignments mafft-nexus-internal-trimmed-gblocks-clean-35p \
  --output mafft-phylip-internal-trimmed-gblocks-clean-35p-raxml \
  --phylip \
  --log-path log
```















