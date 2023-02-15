# ML phylogeny: settings and parameters

For the Maximum Likelihood phylogeny, we used RAxML (Stamatakis, 2014)

```python
./raxmlHPC-PTHREADS-AVX2 -m GTRGAMMAX -p 1505 -s 'redsea-complete'INPUT.phylip -n OUTPUT -f a -x 3012 -T 16 -# 1000
```

phylip files can be found in the *data_insilico_test*  and *data_invitro_test* folders.

## Reference
Stamatakis, A. (2014). RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312–1313. https://doi.org/10.1093/bioinformatics/btu033
