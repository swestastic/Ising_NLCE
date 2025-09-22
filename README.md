# 2D Ising Model NLCE

This repository implements the Numerical Linked Cluster Expansion (NLCE) method for the 2D classical Ising model on the square lattice with nearest neighbor interactions. NLCE allows us to calculate extensive observables in the thermodynamic limit $L \rightarrow \infty$.

This method is great for generating results that don't suffer from finite size effects, although it struggles at low temperatures because the clusters are not large enough to correctly capture the system's long correlation lengths. Better results can be had with large orders, although that leads to high memory usage and long runtimes.

## Background

The 2D Ising model is the simplest spin model which shows spontenaeous symmetry-breaking phase transition. This transition occurs at $T_c = 2.269185...$
For simplicity, we assume no external magnetic field, $h=0$.

This gives us the following Hamiltonian to describe the system:

```math
H = -J \sum_i\sigma_i\sigma_j
```

## NLCE Implementation

### Key Terms

- **Order:** Number of sites in a cluster
- **Symmetrically Distinct Clusters:** Unique "base" clusters which cannot be defined in terms of other symmetrically distinct clusters in the same order.
- **Topologically Distinct Clusters:** Clusters which share the same Hamiltonian. Although there are symmetrically distinct, they have 
- **Multiplicity:** Number of possible embeddings of a specific symmetrically distinct cluster in the lattice.
- **Weight**: Number of possible embeddings of a specific smaller order cluster in a larger order cluster. 

### Cluster Creation

To begin, we want to calculate all of the possible valid clusters for each order we want to use. The order refers to the number of sites in each cluster, so an order 1 cluster will have 1 site, order 2 will have 2 sites, etc.

Something to consider is the difference between all clusters, symmetrically distinct clusters, and topologically distinct clusters; This is discussed later. A table containing the number of clusters is below for the square lattice. 

| Order | Total Clusters (naive) | Symmetrically Distinct | Topologically Distinct |
| :---: | :---: | :---: | :---: |
| 1 | 1 | 1 | 1 |
| 2 | 4 | 1 | 1 |
| 3 | 6 | 2 | 1 |
| 4 | 15 | 5 | 3 |
| 5 |43 | 12 | 4 |
| 6 | 119 | 35 | 10 |
| 7 | 389 | 108 | 19 |
| 8 | 1343 | 369 | 51 |
| 9 | 5029 | 1285 | 112 |
| 10 | 19126 | 4655 | 300 |

We start with the order 1 cluster of 1 site at $(0,0)$. The order 2 clusters are then created by adding all possible neighboring sites to each of the previous order's clusters. So a naive generation of order 2 clusters will give 4 clusters as follows:

$[(0,0), (0,1)], [(0,0), (1,0)], [(0,1),(0,0)], [(1,0),(0,0)]$

We do not care about specific locations of sites in the cluster, only all of the combined locations of sites. For example looking at an order 2 cluster located at $[(0,0), (1,0)]$ with 2 sites $S_1, S_2$, It does not matter if $S_1$ is located at $(0,0)$ or $(1,0)$, these are considered to be the same cluster. For order 2, this reduces our 4 generated clusters to only 2 clusters, as follows:

$[(0,0), (1,0)], [(0,0), (0,1)]$

After calculating all possible clusters for an order (the number in "Total Clusters (naive)"), we then remove redundant clusters to reduce us to the number of clusters in "Symmetrically Distinct", which are then used to calculate clusters for the next order.

### Multiplicity Calculation

The multiplicity of a cluster, $L(c)$, is given by the number of different ways a certain cluster can be orientated in the lattice, ignoring specific site locations. These are called **embeddings**. As discussed previously, we have 2 distinct order 2 clusters:

$[(0,0), (0,1)], [(0,0), (1,0)]$.

If we examine the clusters, we see that they are effectively the same cluster, just rotated $90\degree$. This means they are **not symmetrically distinct**. Symmetrically distinct clusters are clusters which cannot be created by transformations and translations of other symmetrically distinct clusters on the lattice through rotations and/or mirroring. The multiplicity for each cluster is given by its number of unique embeddings. For order 2, these two clusters yield 1 symmetrically distinct cluster with a multiplicity of 2.

### Bond Calculation

After each cluster is generated, we iterate through all sites in the cluster and look at each site's neighboring sites. If a neighboring site is in the cluster, then a bond is added between the two sites. We then repeat this for every site in the cluster until all possible bonds are calculated, while being careful to not double count any bonds.

### Cluster Energy Calculation

Each cluster has $2^{\text{order}}$ possible configurations of spins, because each site has 2 possible configurations ($\pm 1$). We first calculate all of these configurations.

We then calculate the partition function $Z$ for each cluster by summing over each configuration of $S$ values: $\{S\} = S_1, S_2, ..., S_\text{order}$ where $S_1, S_2, ... , S_\text{order} = \pm 1$:

```math
Z_\text{cluster} = \sum_{\{S\}} e^{-\beta H \{S\}} = \sum_{\{S\}}e^{\beta J \sum_{\langle i, j \rangle} S_i S_j}
```

For example, the partition function for the order 2 cluster could be calculated as follows:

$[S_1,S_2] =
[[1,1],
[1,-1],
[-1,1],
[-1,-1]]$

($2^2 = 4$ configurations).

```math
Z_\text{cluster} = e^{\beta J (1 \cdot 1)} + e^{\beta J (1 \cdot -1)} + e^{\beta J (-1 \cdot 1)} + e^{\beta J (-1 \cdot -1)} 
```

```math
Z_\text{cluster} = e^{\beta J} + e^{-\beta J} + e^{-\beta J} + e^{\beta J} 
```

```math
Z_\text{cluster} = 2e^{\beta J} + 2e^{-\beta J}
```

```math
Z_\text{cluster} = 4\text{cosh}(\beta J)
```

Now that we have the partition function, we can calculate the energy of each cluster by summing over all of the bonds in the cluster for each configuration:

```math
E_\text{cluster} = \frac{\sum_{\{S\}}\sum_\text{bonds} H_\text{bond in \{S\}} e^{-\beta H_\text{bond in \{S\}}}}
{Z_\text{cluster}}
 = \frac{\sum_{\{S\}}\sum_\text{bonds} -J S_i S_j e^{\beta J S_i S_j}}
 {Z_\text{cluster}}
```

For order 2 there is only one bond between $(0,0)$ and $(1,0)$, so the calculation is as follows:

```math
E_\text{cluster} = \frac{1}{Z_\text{cluster}}((1)(1)e^{\beta J (1)(1)}
 + (1)(-1)e^{\beta J (1)(-1)}
 + (-1)(1)e^{\beta J (-1)(1)}
 + (-1)(-1)e^{\beta J (-1)(-1)})
```

```math
E_\text{cluster} = \frac{-J}{4\text{cosh}(\beta J)}(e^{\beta J} - e^{-\beta J} - e^{-\beta J} + e^{\beta J})
```

```math
E_\text{cluster} = \frac{-J}{4\text{cosh}(\beta J)}(2e^{\beta J} - 2e^{-\beta J})
```

```math
E_\text{cluster} = \frac{-4J\text{sinh}(\beta J)}{4\text{cosh}(\beta J)}
```

```math
E_\text{cluster} = -J\text{tanh}(\beta J)
```

### Weight Calculation

The weight for each cluster is given by its number of embeddings in a cluster of a higher order. For example, the order 1 cluster can be embedded in the order 2 cluster twice, since it is a single site and the order 2 cluster contains 2 sites. Similarly, the order 2 cluster, which is 2 sites connected by a single bond, can be embedded twice in each order 3 cluster. The first order 3 cluster is a straight line, and the second forms a right angle. In both cases there are two places were an order 2 cluster can be placed.

The weights are calculated differently for each extensive property that we want to use, as follows:

```math
W_p(c) = P(c) - \sum_{s \subset c} W_p(s)
```

For the energy specificially:

```math
W_E(\text{cluster}) = E(\text{cluster}) - \sum_\text{embedded clusters} W_E(\text{embedded cluster }) 
```

### Total Energy Calculation

Once we have calculated the weights and multiplicities, we can calculate the property per site across the whole lattice as follows:

```math
P(\mathcal{L})/N = \sum_\text{clusters} L(\text{cluster}) \times W_p(c)
```

For the energy specifically:

```math
E/N = \sum_\text{clusters} L_\text{cluster} \times W_E(\text{cluster})
```

### Resummation Techniques

There are two resummation techniques we can look at which accelerate the convergence of NLCE, allowing us to reach convergence at lower temperatures per order: Wynn and Euler.

These are not yet implemented, but are discussed in reference [1].

### Runtimes

Calculating everything (clusters, bonds, energies, etc.) for up to order 9 takes about 167 minutes on my personal machine (old, less efficient cluster generation)

For up to order 8, it takes about 24 minutes to calculate everything from scratch.

### TODO

- The Hamiltonian is currently recomputed for all clusters, although it only *needs* to be recomputed for topologically distinct clusters.

- Allow for loading of clusters, energies, bonds, etc. This way if we have previously computed up to order 8, then to compute order 9 we do not need to recalculate the lower orders.

- Implement resummation techniques.

- Possible bug? Orders 2&3, 4&5, 6&7, ... return near identical results. Even orders create more bonds, so it may make sense that going from even $\rightarrow$ odd order would give better results than going from odd $\rightarrow$ even order, although I'm not sure.

- Update lists to NumPy arrays.

### Acknowledgements

Many thanks to Francisco, who works with much more complicated implementations of NLCE and gave me a lot of good information about the method and its implementation.

### References

[1] A Short Introduction to Numerical Linked-Cluster Expansions - Baoming Tang, Ehsan Khatami, Marcos Rigol
[arXiv:1207.3366 [cond-mat.stat-mech]](https://arxiv.org/abs/1207.3366).