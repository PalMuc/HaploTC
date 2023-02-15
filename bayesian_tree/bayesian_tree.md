# RevBayes Script

The bayesian trees were computed using RevBayes (Hönha et al. 2016)
The file should be saved accordingly (.Rev). See below for all the details.

## Reading in the data

Add the path to your file

```python
data <- readDiscreteCharacterData("INPUT.nexus")
```

Get some useful variables from the data. We need these later on

```python
n_species <- data.ntaxa()
n_sites <- data.nchar()
names <- data.names()
n_branches <- 2 * n_species - 3
```

Create some vector for the moves and monitors of this analysis

```python
moves    = VectorMoves()
monitors = VectorMonitors()
```

## Substitution Model 


Specify the stationary frequency parameters

```python
pi_prior <- v(1,1,1,1) 
pi ~ dnDirichlet(pi_prior)
moves.append( mvBetaSimplex(pi, weight=2.0) ) 
moves.append( mvDirichletSimplex(pi, weight=2.0) ) 
```

Specify the exchangeability rate parameters

```python
er_prior <- v(1,1,1,1,1,1)
er ~ dnDirichlet(er_prior)
moves.append( mvBetaSimplex(er, weight=3.0) ) 
moves.append( mvDirichletSimplex(er, weight=3.0) )
```

Create a deterministic variable for the rate matrix, GTR

```python
Q := fnGTR(er,pi)
```

Among Site Rate Variation

* Among site rate variation, +Gamma4

```python
alpha ~ dnUniform( 0, 1E6 )
alpha.setValue( 10 )
gamma_rates := fnDiscretizeGamma( alpha, alpha, 4, false )

sr := fnDiscretizeGamma( alpha, alpha, 4, false )
moves.append( mvScale(alpha, lambda=1.0, weight=3.0) ) 
```


## Tree model 

Specify a uniform prior on the tree topology #### 

```python
phylogeny ~ dnUniformTopologyBranchLength(branchLengthDistribution=dnExponential(10.0), names)
```

Moves on the tree

```python
moves.append( mvNNI(phylogeny, weight=n_species*5) )
moves.append( mvSPR(phylogeny, weight=n_species*5/10) )
moves.append( mvBranchLengthScale(phylogeny, weight=n_species*5) )

TL := phylogeny.treeLength()
```


## PhyloCTMC Model 

The sequence evolution model

```python
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, siteRates=gamma_rates, type="DNA")
```

Attach the data

```python
seq.clamp(data)
```


## THE Model 

We define our model.
##### We can use any node of our model as a handle, here we chose to use the rate matrix.

```python
mymodel = model(Q)

monitors = VectorMonitors()

monitors.append( mnModel(filename = "output_bayesian_tree_posterior.log",printgen=1, separator = TAB) )
monitors.append( mnFile(filename = "output_bayesian_tree_posterior.trees",printgen=1, separator = TAB, phylogeny) )
monitors.append( mnScreen(printgen=1000, TL) )


mymcmc = mcmc(mymodel, monitors, moves, nruns=4, combine="mixed")
mymcmc.burnin(generations=20000,tuningInterval=200)
mymcmc.run(generations=200000,tuningInterval=100)

q()
```

## Check for convergence

The outcome of this script will result in 4 different trees, which need to be checked for convergence using the following Rscript:

```python
library(convenience)

check_haplos <- checkConvergence(list_files = c("output_bayesian_tree_posterior_run_1.log", "output_bayesian_tree_posterior_run_1.trees",
						"output_bayesian_tree_posterior_run_2.log", "output_bayesian_tree_posterior_run_2.trees",
						"output_bayesian_tree_posterior_run_3.log", "output_bayesian_tree_posterior_run_3.trees",
						"output_bayesian_tree_posterior_run_4.log", "output_bayesian_tree_posterior_run_4.trees"))

check_haplos


output$message_complete
```

* More details can be found: https://github.com/lfabreti/convenience
* For the demosponge tree, we increased the moves on the tree to 9 reach convergence

## Final tree

To create the final bayesian tree we used the following RevBayes script:

```python
trace = readTreeTrace("output_bayesian_tree_posterior.trees", treetype="non-clock")
trace.setBurnin(0.25)

# Summarize tree trace and the consensus tree to file
map_tree = mapTree(trace, file="redsea_bayesian.map.tre")
con_tree = consensusTree(trace, file="redsea_bayesian.majrule.tre")
```



## References
Höhna, S., Landis, M. J., Heath, T. A., Boussau, B., Lartillot, N., Moore, B. R., … Ronquist, F. (2016). RevBayes: Bayesian phylogenetic inference using graphical models and an interactive model-specification language. Systematic Biology, 65(4), 726–736. https://doi.org/10.1093/sysbio/syw021