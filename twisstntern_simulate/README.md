# Running `twisstntern_simulate`

## Why simulate data?
The `twisstntern` module enables the implementation of a wide range of neutral demographic simulations. Simulating data allows researchers to explore how different demographic scenarios—such as population divergence, migration, and changes in population size—affect the distribution of genealogical topologies. These simulated patterns can be visualized using ternary plots, helping to build intuition about expected signals under neutrality and to compare against empirical genomic data.

## Configuration file example


```
###########################################
# CONFIGURATION FILE FOR MSPRIME SIMULATION
# Four-population demographic model
# You can modify the values below or leave them as defaults
###########################################

#=== SPLIT TIMES (in generations before present) ===#
# These define the timing of historical splits between populations.
# Larger values indicate older divergence events.

t1 = 1000         # Time of split between p1 and p2 (forms ancestor p12)
t2 = 5000         # Time of split between p12 and p3 (forms ancestor p123)
t3 = 10000        # Time of split between p123 and outgroup O

#=== EFFECTIVE POPULATION SIZES (diploid individuals) ===#
# Ne values for each extant and ancestral population.
# Affects genetic drift and coalescent times.

ne_p1 = 10000     # Population size of extant population 1
ne_p2 = 10000     # Population size of extant population 2
ne_p3 = 10000     # Population size of extant population 3
ne_p12 = 15000    # Size of ancestral population of p1 and p2
ne_p123 = 20000   # Size of ancestral population of p1, p2, and p3
ne_O = 30000      # Size of outgroup population

#=== MIGRATION RATES (per generation) ===#
# These rates are proportions: e.g. 0.001 means 0.1% of individuals replaced by migrants per generation.
# Define direction as "source > destination"
# Default m = 0

# Migration between extant populations
m_p1>p2 = 0.001
m_p2>p1 = 0.001
m_p1>p3 = 0.0005
m_p3>p1 = 0.0005
m_p2>p3 = 0.0005
m_p3>p2 = 0.0005

# Migration between outgroup and extant populations
m_O>p1 = 0.0001
m_p1>O = 0.0001
m_O>p2 = 0.0001
m_p2>O = 0.0001
m_O>p3 = 0.0001
m_p3>O = 0.0001

# Migration involving ancestral populations
m_p12>p3 = 0.0002   # Migration from ancestor of p1/p2 to p3 (before t2)
m_p3>p12 = 0.0002
m_O>p123 = 0.00005  # Migration from outgroup to ancestral p123 (before t3)
m_p123>O = 0.00005

###########################################
# END OF CONFIGURATION
###########################################

```


> 💡 **Note for developers**  
> Implemented in ...