# PHYLUCE Bait Design: parameters and settings

Settings used for the programs and Python script for our bait set design following the online tutorial of PHYLUCE IV:
https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html (Faircloth, 2012; Faircloth et al., 2016)

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

Put transcriptomes in their own directories

```python
> cd uce-haplosclerida/transcriptomes
> for critter in *; do mkdir ${critter%.*}; mv $critter ${critter%.*}; done
```

Convert transcriptome to 2bit format

```python
> for critter in *; do faToTwoBit $critter/$critter.fasta $critter/${critter%.*}.2bit; done
```

Simulate reads from transcriptomes

* Install 'Art' Package (Huang et al. 2012)

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

Merging reads together

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

Prepare AMQU transcriptome for use, with stampy (Lunter & Goodson, 2011)

* Install package 'Stampy v.1.0.32'

```python
> make python=python2.7

> cd uce-haplosclerida
> cd base

> python2 ../stampy-1.0.32/stampy.py --species="Amphimedon queenslandica" --assembly="AMQU" -G AMQU AMQU.fasta
> python2 ../stampy-1.0.32/stampy.py -g AMQU -H AMQU
```

Align reads to base transcriptome

```python
> cd uce-haplosclerida
> mkdir alignments
```

Perform alignments (taxon-by-taxon)

* Substitutionrate = 0.05 (sequence divergence <5%)

```python
> export cores=4
> export base=AMQU
> export base_dir=../uce-haplosclerida/alignments
for critter in HCIN HIND HOCU HSIM HTUB HVIS NCOM;     
do         
 export reads=$critter-pe100-reads.fq.gz;         
 mkdir -p $base_dir/$critter;         
 cd $base_dir/$critter;
 python2 ../stampy-1.0.32/stampy.py --maxbasequal 93 -g ../../base/$base -h ../../base/$base \
 --substitutionrate=0.05 -t$cores --insertsize=400 -M \
 ../../reads/$reads | samtools view -Sb - > $critter-to-$base.bam;     
done;
```

Remove unmapped reads

```python
> cd uce-haplosclerida
> cd alignments
> mkdir all

> for critter in NCOM HCIN HIND HOCU HSIM HVIS HTUB;
do
 samtools view -h -F 4 -b $critter/$critter-to-AMQU.bam > $critter/$critter-to-AMQU-MAPPING.bam;
 rm $critter/$critter-to-AMQU.bam;
 ln -s ../$critter/$critter-to-AMQU-MAPPING.bam all/$critter-to-AMQU-MAPPING.bam;
done;
```


### CONSERVED LOCI IDENTIFICATION

Convert BAMS to BEDS

```python
> cd uce-haplosclerida
> mkdir bed
> cd bed

> for i in ../alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done
```

Sort the converted BEDS

```python
> for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done
```

Merge (nearly) overlapping intervals

```python
> for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done
```

Total N of merged, putatively conserved regions in each exemplar taxon shared with the base transcriptome

```python
> for i in *.bam.sort.merge.bed; do wc -l $i; done
```

Remove repetitive intervals

* base genome shorter than 80-bp
* where more than 25% of the base transcriptome is masked

```python
> for i in *.sort.merge.bed; do phyluce_probe_strip_masked_loci_from_set --bed $i --twobit ../transcriptomes/AMQU/AMQU.2bit --output ${i%.*}.strip.bed --filter-mask 0.25 --min-length 80; done;
```

Determining locus presence in multiple transcriptomes
* create a configuration file > bed-files.conf

```python
[beds]
HIND:HIND-to-AMQU-MAPPING.bam.sort.merge.strip.bed
HOCU:HOCU-to-AMQU-MAPPING.bam.sort.merge.strip.bed
HCIN:HCIN-to-AMQU-MAPPING.bam.sort.merge.strip.bed
HTUB:HTUB-to-AMQU-MAPPING.bam.sort.merge.strip.bed
HVIS:HVIS-to-AMQU-MAPPING.bam.sort.merge.strip.bed
HSIM:HSIM-to-AMQU-MAPPING.bam.sort.merge.strip.bed
NCOM:NCOM-to-AMQU-MAPPING.bam.sort.merge.strip.bed
```

Create record of alignment intervals shared among taxa

```python
> phyluce_probe_get_multi_merge_table \
  --conf bed-files.conf \
  --base-taxon AMQU \
  --output Haplosclerida-to-AMQU.sqlite
```

Table output

```python
> sqlite3 Haplosclerida-to-AMQU.sqlite
> .tables
> select * from AMQU limit 10;
```

Determining shared, conserved loci

* query the number of shared loci by the base transcriptome and the other exemplar taxa
* the results are documented in Table S3

```python
> phyluce_probe_query_multi_merge_table \
  --db Haplosclerida-to-AMQU.sqlite \
  --base-taxon AMQU
```

Here we selected those loci being shared between base transcriptome and all other exemplar taxa (n = 7)

```python
> phyluce_probe_query_multi_merge_table \
  --db Haplosclerida-to-AMQU.sqlite \
  --base-taxon AMQU \
  --output AMQU+7.bed \
  --specific-counts 7
```

### CONSERVED LOCUS VALIDATION

Extract fasta sequence from base transcriptome for temporary bait design

```python
> phyluce_probe_get_genome_sequences_from_bed \
	--bed AMQU+7.bed \
	--twobit ../transcriptomes/AMQU/AMQU.2bit \
	--buffer-to 120 \
	--output AMQU+7.fasta
```

Design duplicated temporary bait set

```python
> phyluce_probe_get_tiled_probes \
	--input AMQU+7.fasta \
	--probe-length 80 \
	--probe-prefix "uce-" \
	--design Haplosclerida-v50 \
	--designer vandersprong \
	--tiling-density 2 \
	--two-probes \
	--overlap middle \
	--masking 0.25 \
	--remove-gc \
	--output AMQU+7.temp.probes
```

Remove duplicated temporary bait set

```python
> phyluce_probe_easy_lastz \
  --target AMQU+7.temp.probes \
  --query AMQU+7.temp.probes \
  --identity 50 \
  --coverage 50 \
  --output AMQU+7.temp.probes-TO-SELF-PROBES.lastz 
```

Screen alignments + remove duplicate baits

```python
> phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
	--fasta AMQU+7.temp.probes \
	--lastz AMQU+7.temp.probes-TO-SELF-PROBES.lastz \
	--probe-prefix=uce- 
```

Align baits against exemplar transcriptomes

* using a duplicate-free (or putatively duplicate free) set of temporary baits designed from conserved loci in the base transcriptome
* how 'sticky' are the designed baits when they are used to enrich loci across divergent groups?

```python
> cd uce-haplosclerida
> mkdir probe-design
> cd probe-design 
```

Align temporary sequences to each transcriptome

```python
> mkdir haplosclerida-transcriptome-lastz

> phyluce_probe_run_multiple_lastzs_sqlite \
	--probefile ../bed/AMQU+7.temp-DUPE-SCREENED.probes \
	--scaffoldlist NCOM HCIN HIND HOCU HSIM HTUB HVIS \
	--genome-base-path ../transcriptomes \
	--identity 70 \
	--cores 4 \
	--db AMQU+7.sqlite \
	--output haplosclerida-transcriptome-lastz
```

The temporary bait set ('x' number of baits, targeting a 'x' number of loci) is here aligned back to AMQU and the 7 exemplar taxa, with an identity value of 70%. This is the minimum sequence identity for which a bait could be an accepted match to the transcriptome + a minimum coverage of 83% (Quattrini et al. 2018).


Extract- sequence around conserved loci from exemplar transcriptomes

* create configuration file 'haplosclerida-transcriptome.conf'

```python
[scaffolds]
NCOM:../transcriptomes/NCOM/NCOM.2bit
HIND:../transcriptomes/HIND/HIND.2bit
HCIN:../transcriptomes/HCIN/HCIN.2bit
HOCU:../transcriptomes/HOCU/HOCU.2bit
HVIS:../transcriptomes/HVIS/HVIS.2bit
HSIM:../transcriptomes/HSIM/HSIM.2bit
HTUB:../transcriptomes/HTUB/HTUB.2bit
AMQU:../transcriptomes/AMQU/AMQU.2bit
```

```python
> phyluce_probe_slice_sequence_from_genomes \
  --conf haplosclerida-transcriptome.conf \
  --lastz haplosclerida-transcriptome-lastz \
  --probes 140 \
  --name-pattern "AMQU+7.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
  --output haplosclerida-transcriptome-fasta
```

Exploring fasta files

```python
> less probe-design/haplosclerida-transcriptome-fasta/AMQU.fasta
```

Find consistently detected loci

```python
> phyluce_probe_get_multi_fasta_table \
	--fastas ../uce-haplosclerida/probe-design/haplosclerida-transcriptome-fasta \
	--output multifastas.sqlite \
	--base-taxon AMQU
```

Visualize table: shows detection of conserved loci in each of the exemplar taxa

```python
> sqlite3 multifastas.sqlite
> sqlite > select * from AMQU limit 10;
```

Detection shared loci (after cleaning) > results can be found in Table S3

```python
> phyluce_probe_query_multi_fasta_table \
  --db multifastas.sqlite \
  --base-taxon AMQU
```

Extract conserved loci

```python
> phyluce_probe_query_multi_fasta_table \
	--db multifastas.sqlite \
	--base-taxon AMQU \
	--output AMQU+7-back-to-7.conf \
	--specific-counts 7
```


### FINAL BAIT SET DESIGN

Bait set design + all exemplar transcriptomes (and base transcriptome)

```python
> phyluce_probe_get_tiled_probe_from_multiple_inputs \
	--fastas haplosclerida-transcriptome-fasta \
	--multi-fasta-output AMQU+7-back-to-7.conf \
	--probe-prefix "uce-" \
	--probe-length 80 \
	--designer vandersprong \
	--design haplosclerida-v50 \
	--tiling-density 2 \
	--overlap middle \
	--masking 0.25 \
	--remove-gc \
	--two-probes \
	--output haplosclerida-v50-amqu-master-probe-list.fasta
```

Remove duplicates from bait set

```python
> phyluce_probe_easy_lastz \
	--target haplosclerida-v50-amqu-master-probe-list.fasta \
	--query haplosclerida-v50-amqu-master-probe-list.fasta \
	--identity 50 \
	--coverage 50 \
	--output haplosclerida-v50-amqu-master-probe-list-TO-SELF-PROBES.lastz
```

Filtering for duplicate loci

```python
> phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
	--fasta haplosclerida-v50-amqu-master-probe-list.fasta \
	--lastz haplosclerida-v50-amqu-master-probe-list-TO-SELF-PROBES.lastz \
	--probe-prefix=uce-
```

After this step, a new file has been generated: 'haplosclerida-v50-amqu-master-probe-list-DUPE-SCREENED.fasta'

## Notes

The different steps as described above were repeated for every bait set (n = 137).

* using different species as base transcriptome
* using different levels of stringency
* using the haplosclerids + one demosponge as outgroup


To keep track of the probes designed by a specific bait set (and avoid different probes having similar name-codes), we gave every UCE a unique code:
* starting with the version number (in this case, the UCEs designed with v50 will start with '50')
* followed by three '000'
* followed by the UCE number resulting from the PHYLUCE workflow

```python
> sed -e 's/uce-/uce-[insert-unique-nr]/g' <input.fasta > output.fasta
```

Once the probes/baits have their unique codes, we compiled the different sets together in a (step-wise).

* we started with concatenating the bait sets
* then we re-ran the program for removing duplicates

```python
> phyluce_probe_easy_lastz \
	--target haplosclerida-master-probe-list-concat-1.fasta \
	--query haplosclerida-master-probe-list-concat-1.fasta \
	--identity 50 \
	--coverage 50 \
	--output haplosclerida-master-probe-list-concat-1-TO-SELF-PROBES.lastz
```

```python
> phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
	--fasta haplosclerida-master-probe-list-concat-1.fasta \
	--lastz haplosclerida-master-probe-list-concat-1-TO-SELF-PROBES.lastz \
	--probe-prefix=uce-
```

After this step, we again had a new master probe list filtered of putatively duplicate loci: 'haplosclerida-master-probe-list-concat-1-DUPE-SCREENED.fasta'



## References
Faircloth BC. 2016. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics 32:786-788. doi:10.1093/bioinformatics/btv646.

Faircloth BC, McCormack JE, Crawford NG, Harvey MG, Brumfield RT, Glenn TC. 2012. Ultraconserved elements anchor thousands of genetic markers spanning multiple evolutionary timescales. Systematic Biology 61: 717–726. doi:10.1093/sysbio/SYS004.

Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics, 28(4), 593–594.

Lunter, G., & Goodson, M. (2011). Stampy: a statistical algorithm for sensitive and fast mapping of Illumina sequence reads. Genome Research, 21(6), 936–939. doi:10.1101/gr.111120.110

Quattrini, A. M., Faircloth, B. C., Dueñas, L. F., Bridge, T. C. L., Brugler, M. R., Calixto-Botía, I. F., … McFadden, C. S. (2018). Universal target-enrichment baits for anthozoan (Cnidaria) phylogenomics: new approaches to long-standing problems. Molecular Ecology Resources, 18(2), 281–295. doi: 10.1111/1755-0998.12736

Rivera-Vicéns RE, García-Escudero CA, Conci N, Eitel M, Wörheide G. 2021. TransPi–a comprehensive TRanscriptome ANalysiS PIpeline for de novo transcriptome assembly. bioRxiv 2021.02.18.431773; doi:10.1101/2021.02.18.431773








